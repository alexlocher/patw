/**
 * This file is part of PATW (Progressive All The Way).
 *
 * Copyright (C) 2016 Alex Locher <alocher at ethz dot ch> (ETH Zuerich)
 * For more information see <https://github.com/alexlocher/patw>
 *
 * PATW is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PATW is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PATW. If not, see <http://www.gnu.org/licenses/>.
 */
#include <iostream>
#include <string>
#include <functional>
#include <sstream>

#include <glog/logging.h>
#include <gflags/gflags.h>
#include <map>

#include <stlplus3/file_system.hpp>

#include <hpmvs/NVMReader.h>
#include <hpmvs/CellProcessor.h>
#include <hpmvs/Scene.h>

#include <omp.h>

#include <CImg.h>
#include <theia/sfm/camera/camera.h>
#include <theia/sfm/triangulation/triangulation.h>

using namespace std;


// ===============================================================

void toPly(const char* name, const std::vector<mo3d::Ppatch3d>& patches, bool binary) {
	int32_t n;
	bool bigEndian = *(char *) &n == 1;

	// output the header
	std::ofstream pFile(name, std::ofstream::out);
	pFile << "ply" << std::endl;
	if (binary && bigEndian)
		pFile << "format binary_big_endian 1.0" << std::endl;
	else if (binary)
		pFile << "format binary_little_endian 1.0" << std::endl;
	else
		pFile << "format ascii 1.0" << std::endl;
	pFile << "element vertex " << (int) patches.size() << std::endl;
	pFile << "property float x" << std::endl;
	pFile << "property float y" << std::endl;
	pFile << "property float z" << std::endl;
	pFile << "property float nx" << std::endl;
	pFile << "property float ny" << std::endl;
	pFile << "property float nz" << std::endl;
	pFile << "property uchar red" << std::endl;
	pFile << "property uchar green" << std::endl;
	pFile << "property uchar blue" << std::endl;
	pFile << "property float scalar_scale" << std::endl;
	pFile << "element point_visibility " << (int) patches.size() << std::endl;
	pFile << "property list uint uint visible_cameras" << std::endl;
	pFile << "end_header" << std::endl;
	pFile.close();

	// output the points
	std::ofstream pData(name,
			(binary ? std::ofstream::binary | std::ofstream::app : std::ofstream::app));

	for (const auto& p : patches) {
		if (binary) {
			Eigen::Vector3f v(p->x(), p->y(), p->z());
			pData.write((char*) v.data(), 3 * sizeof(float));

			v = p->normal_.head(3);
			pData.write((char*) v.data(), 3 * sizeof(float));

			Eigen::Matrix<unsigned char, 3, 1> c(p->color_[0], p->color_[1], p->color_[2]);
			pData.write((char*) c.data(), 3 * sizeof(unsigned char));

			pData.write((char*) &p->scale_3dx_, sizeof(p->scale_3dx_));

		} else {
			pData << p->x() << " " << p->y() << " " << p->z() << " ";
			pData << p->normal_[0] << " " << p->normal_[1] << " " << p->normal_[2] << " ";
			pData << (int) p->color_[0] << " " << (int) p->color_[1] << " " << (int) p->color_[2]
					<< " ";
			pData << p->scale_3dx_ << " ";
			pData << std::endl;
		}
	}

	// now output the visibility
	for (const auto& p : patches) {
		if (binary) {
			uint32_t nrImgs = p->images_.size();
			pData.write((char*) &nrImgs, sizeof(uint32_t));
			for (uint32_t imgId : p->images_)
				pData.write((char*) &imgId, sizeof(uint32_t));
		} else {
			pData << (int) p->images_.size() << " ";
			for (int imgId : p->images_)
				pData << (uint32_t) imgId << " ";
			pData << std::endl;
		}
	}
	pData.flush();
	pData.close();
}

// ==============================================================================

int main(int argc, char* argv[]) {
	if (argc < 2)
		exit(EXIT_FAILURE);

	int minPoints =  (argc > 2) ? std::atoi(argv[2]) : 2;

	// load the nvm file
	std::vector<mo3d::NVM_Model> nvm;
	mo3d::NVMReader::readFile(argv[1], nvm, true);
	CHECK(nvm.size() > 0) << " no models found in <" << argv[1] << ">";

	//
	std::vector<mo3d::Ppatch3d> patches;

	patches.clear();
	for (size_t ii = 0; ii < nvm[0].points.size(); ii++) {
		const mo3d::NVM_Point& pt=nvm[0].points[ii];

		if (pt.measurements.size() < minPoints)
			continue;

		patches.emplace_back(new mo3d::Patch3d);
		patches.back()->center_ << pt.xyz.cast<float>(), 1.0f;
		patches.back()->color_ << pt.rgb.cast<float>();

		// visible images
		for (auto imgIdx : pt.measurements)
			patches.back()->images_.emplace_back(imgIdx.imgIndex);
	}

// --------------------------------------------------
	std::string outfolder = stlplus::folder_part(argv[1]);
	std::string basename = stlplus::basename_part(argv[1]);


	toPly(stlplus::create_filespec(outfolder,basename, "ply").c_str(), patches, false);

}

