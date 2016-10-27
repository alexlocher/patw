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
#include <hpmvs/doctree.h>
#include <patw/EstimateSE3.hpp>

#include <omp.h>

#include <CImg.h>
#include <theia/sfm/camera/camera.h>
#include <theia/sfm/triangulation/triangulation.h>

#include <yaply/YaPLY.hpp>

#include <pcl/point_cloud.h>
#include <pcl/octree/octree.h>

using namespace std;

// FLAGS
DEFINE_string(nvm1, "", "old camera calibration file (as nvm)");
DEFINE_string(nvm2, "", "new camera calibration file (as nvm)");
DEFINE_string(patches, "", "original patches file");
DEFINE_bool(ismesh, false, "the patches file is actually a mesh file");
DEFINE_double(maxdist, 1e-3, "maximum mean distance for the recursive transformation estimation");
DEFINE_double(maxtriang, 10, "maximum triangulation error in order vertex is not deleted");

DEFINE_string(outdir, "/tmp/hpmvs", "output directory");
DEFINE_int32(nn, 4, "number nearest neighbors for consistency estimation");
DEFINE_bool(skip_estimation, false, "if set to true, points are only regtriangulated");
DEFINE_bool(skip_clean, true,
		"skip patches which are clean in HPMVS (good for finalized pointclouds)");
DEFINE_bool(skip_dirtyfromvisible, false , "if true, dirty flag is not set on patches visible in new Cameras");

DEFINE_int32(subtrees, 100, "min number subtrees the model is split into");
DEFINE_int32(maxprio, 1000, "maximum priority until exit of algo");

// ===============================================================

template<class Element>
void getSubTrees(DynOctTree<Element>& tree,
		std::vector<std::shared_ptr<DynOctTree<Element> > >& subTrees, const int minTrees = 2) {

	if (minTrees < 2) {
		subTrees.emplace_back(std::shared_ptr<DynOctTree<Element> >(&tree));
		return;
	}

	// do a first split
	tree.getSubTrees(subTrees);

	while (subTrees.size() < minTrees) {

		// get the subtree with the most leafs
		int maxIndex = -1;
		int maxLeafs = -1;
		for (int ii = 0; ii < subTrees.size(); ii++) {
			int nrLeafs = subTrees[ii]->getRoot()->nrLeafs(); // recursive and slow!
			if (nrLeafs > maxLeafs) {
				maxLeafs = nrLeafs;
				maxIndex = ii;
			}
		}

		// doesn't make sense to split, if we have to few points
		if (maxLeafs < 100)
			return;

		std::shared_ptr<DynOctTree<Element> > maxTree = subTrees[maxIndex];
		std::vector<std::shared_ptr<DynOctTree<Element> >> newSubTrees;
		maxTree->getSubTrees(newSubTrees); // works because we are working with pointers
		for (int ii = 0; ii < subTrees.size(); ii++)
			if (ii != maxIndex)
				newSubTrees.push_back(subTrees[ii]);

		subTrees.swap(newSubTrees);

	}

	LOG(INFO)<< "Split to " << subTrees.size() << " subtrees";
	for (int ii = 0; ii < subTrees.size(); ii++) {
		LOG(INFO)<< "Nr " << ii << " => " << subTrees[ii]->getRoot()->nrLeafs() << " Leafs";
	}
}

// ===============================================================

theia::Camera nvmToTheia(const mo3d::NVM_Camera& nvmCam) {
	theia::Camera c;

	// set extrinsics
	Eigen::Quaterniond rq(nvmCam.rq[0], nvmCam.rq[1], nvmCam.rq[2], nvmCam.rq[3]);
	c.SetOrientationFromRotationMatrix(rq.matrix());
	c.SetPosition(nvmCam.c);

	// set the intrinsics
//	c.SetCameraIntrinsicsModelType(theia::CameraIntrinsicsModelType::PINHOLE_RADIAL_TANGENTIAL);
	c.SetFocalLength(nvmCam.f);
//	c.SetRadialDistortion(nvmCam.r, 0);

	// load the image to get the width / height
	cimg_library::CImg<float> img;
	img.load(nvmCam.filename.c_str());
	c.SetImageSize(img.width(), img.height());
	c.SetPrincipalPoint(img.width() / 2.0, img.height() / 2.0);

//	LOG(INFO)<< "image size: " << img.width() << " " << img.height();
//	LOG(INFO)<< "focal: " << nvmCam.f;

	return c;

}

// ==============================================================================


bool saveColoredMesh(const std::string& plyFile_in, const std::vector<bool>& dirtyVertices) {
	// load the ply file
	if (!stlplus::file_readable(plyFile_in))
		return false;
	yaply::PlyFile ply(plyFile_in.c_str());

	auto ply_r = ply["vertex"].getProperty<yaply::PLY_PROPERTY_SCALAR<unsigned char>>("red");
	auto ply_g = ply["vertex"].getProperty<yaply::PLY_PROPERTY_SCALAR<unsigned char>>("green");
	auto ply_b = ply["vertex"].getProperty<yaply::PLY_PROPERTY_SCALAR<unsigned char>>("blue");
	auto ply_a = ply["vertex"].getProperty<yaply::PLY_PROPERTY_SCALAR<unsigned char>>("alpha");

	CHECK(ply_r && ply_g && ply_b && ply_a) << "vertex colors missing";

	CHECK_EQ(ply_r->data.size(), dirtyVertices.size()) << "number of flags in dirtyVertices does not match ply file";

	// mark dirty vertices as red
	for (int ii = 0; ii < dirtyVertices.size(); ii++){
		if (dirtyVertices[ii]){
			ply_r->data[ii] = 255;
			ply_g->data[ii] = 0;
			ply_b->data[ii] = 0;
			ply_a->data[ii] = 255;
		} else
			ply_a->data[ii] = 128;
	}

	// save it
	LOG(INFO) << "save recolored ply file";
	ply.save(stlplus::create_filespec(FLAGS_outdir, stlplus::basename_part(plyFile_in)+"-colored", "ply").c_str(), false);

	return true;
}

bool loadPatchesMeshPly(const std::vector<theia::Camera>& oldCameras, const std::string& plyFile,
		std::vector<mo3d::Ppatch3d>& patches) {
	if (!stlplus::file_readable(plyFile))
		return false;
	yaply::PlyFile ply(plyFile.c_str());

	// get pointers to the elements
	auto ply_x = ply["vertex"].getProperty<yaply::PLY_PROPERTY_SCALAR<float>>("x");
	auto ply_y = ply["vertex"].getProperty<yaply::PLY_PROPERTY_SCALAR<float>>("y");
	auto ply_z = ply["vertex"].getProperty<yaply::PLY_PROPERTY_SCALAR<float>>("z");

	auto ply_nx = ply["vertex"].getProperty<yaply::PLY_PROPERTY_SCALAR<float>>("nx");
	auto ply_ny = ply["vertex"].getProperty<yaply::PLY_PROPERTY_SCALAR<float>>("ny");
	auto ply_nz = ply["vertex"].getProperty<yaply::PLY_PROPERTY_SCALAR<float>>("nz");

	auto ply_r = ply["vertex"].getProperty<yaply::PLY_PROPERTY_SCALAR<unsigned char>>("red");
	auto ply_g = ply["vertex"].getProperty<yaply::PLY_PROPERTY_SCALAR<unsigned char>>("green");
	auto ply_b = ply["vertex"].getProperty<yaply::PLY_PROPERTY_SCALAR<unsigned char>>("blue");

	auto ply_s = ply["vertex"].getProperty<yaply::PLY_PROPERTY_SCALAR<float>>("radius");

	if (!(ply_x && ply_y && ply_z && ply_r && ply_g && ply_b && ply_s)) {
		LOG(WARNING)<< "Missing property";
		return false;
	}

	// ok, now we can start copy the ply data int the patch data
	for (size_t ii = 0; ii < ply_x->data.size(); ii++) {

		patches.emplace_back(new mo3d::Patch3d);
		patches.back()->center_ << ply_x->value(ii), ply_y->value(ii), ply_z->value(ii), 1.0f;
		patches.back()->color_ << ply_r->value(ii), ply_g->value(ii), ply_b->value(ii);
		patches.back()->normal_ << ply_nx->value(ii), ply_ny->value(ii), ply_nz->value(ii), 1.0f;

		patches.back()->scale_3dx_ = ply_s->value(ii);

		// assume that the camera is visible if we can project the point into
		// the camera
		// TODO => here we could also check the normal
		for (int imgIdx = 0; imgIdx < oldCameras.size(); imgIdx++) {
			Eigen::Vector2d pixel;
			if (oldCameras[imgIdx].ProjectPoint(patches.back()->center_.cast<double>(), &pixel) < 0)
				continue;
			if (pixel[0] < 0 || pixel[0] >= oldCameras[imgIdx].ImageWidth())
				continue;
			if (pixel[1] < 0 || pixel[1] >= oldCameras[imgIdx].ImageHeight())
				continue;
			patches.back()->images_.emplace_back(imgIdx);
		}
	}

	return patches.size() > 0;
}

bool loadPatchesPly(const char* plyFile, std::vector<mo3d::Ppatch3d>& patches) {
	if (!stlplus::file_readable(plyFile))
		return false;
	yaply::PlyFile ply(plyFile);

	// get pointers to the elements
	auto ply_x = ply["vertex"].getProperty<yaply::PLY_PROPERTY_SCALAR<float>>("x");
	auto ply_y = ply["vertex"].getProperty<yaply::PLY_PROPERTY_SCALAR<float>>("y");
	auto ply_z = ply["vertex"].getProperty<yaply::PLY_PROPERTY_SCALAR<float>>("z");

	auto ply_nx = ply["vertex"].getProperty<yaply::PLY_PROPERTY_SCALAR<float>>("nx");
	auto ply_ny = ply["vertex"].getProperty<yaply::PLY_PROPERTY_SCALAR<float>>("ny");
	auto ply_nz = ply["vertex"].getProperty<yaply::PLY_PROPERTY_SCALAR<float>>("nz");

	auto ply_r = ply["vertex"].getProperty<yaply::PLY_PROPERTY_SCALAR<unsigned char>>("red");
	auto ply_g = ply["vertex"].getProperty<yaply::PLY_PROPERTY_SCALAR<unsigned char>>("green");
	auto ply_b = ply["vertex"].getProperty<yaply::PLY_PROPERTY_SCALAR<unsigned char>>("blue");

	auto ply_s = ply["vertex"].getProperty<yaply::PLY_PROPERTY_SCALAR<float>>("scalar_scale");

	auto ply_vis =
			ply["point_visibility"].getProperty<yaply::PLY_PROPERTY_LIST<uint32_t, uint32_t>>(
					"visible_cameras");

	if (!(ply_x && ply_y && ply_z && ply_r && ply_g && ply_b && ply_s && ply_vis)) {
		LOG(WARNING)<< "Missing property";
		return false;
	}

	// ok, now we can start copy the ply data int the patch data
	for (size_t ii = 0; ii < ply_x->data.size(); ii++) {

		patches.emplace_back(new mo3d::Patch3d);
		patches.back()->center_ << ply_x->value(ii), ply_y->value(ii), ply_z->value(ii), 1.0f;
		patches.back()->color_ << ply_r->value(ii), ply_g->value(ii), ply_b->value(ii);
		patches.back()->normal_ << ply_nx->value(ii), ply_ny->value(ii), ply_nz->value(ii), 1.0f;

		patches.back()->scale_3dx_ = ply_s->value(ii);

		// visible images
		for (int imgIdx : ply_vis->data[ii])
			patches.back()->images_.emplace_back(imgIdx);
	}

	return patches.size() > 0;
}


// ###################################################################

// struct containing a 3d point plus an idx

struct PointWithIdx {
	PointWithIdx(Eigen::Vector3f pt, unsigned int idx) :
			pt_(pt), idx_(idx) {
	}
	inline float x() {
		return pt_[0];
	}
	inline float y() {
		return pt_[1];
	}
	inline float z() {
		return pt_[2];
	}
	Eigen::Vector3f pt_;
	unsigned int idx_;
};

float estimateCell(Leaf<shared_ptr<PointWithIdx>>* leaf, std::vector<Eigen::Vector3f>& oldPoints,
		std::vector<Eigen::Vector3f>& newPoints, std::vector<Sophus::Sim3f>& transformations,
		float maxError) {

	// estimate an SE3 transformation using all elements in this leaf
	const int nrPoints = leaf->data.size();
	if (nrPoints < 3) {
		return 0.0;
	}

	std::vector<Eigen::Vector3f> pl(nrPoints), pr(nrPoints);
	for (int ii = 0; ii < nrPoints; ii++) {
		pl[ii] = oldPoints[leaf->data[ii]->idx_];
		pr[ii] = newPoints[leaf->data[ii]->idx_];
	}

	Sophus::Sim3f sim3;
	if (!mo3d::Sim3fromCorrspondence(pl, pr, sim3)) {
		LOG(WARNING)<< "SE3 estimation failed";
	}

//	// estimate the transformation
	std::vector<unsigned int> idxs;
	for (auto pt : leaf->data)
		idxs.emplace_back(pt->idx_);
//	theia::RansacSummary summary;
//	Sophus::Sim3f sim3;
//	float error = mo3d::Sim3_ransac(oldPoints, newPoints, idxs, sim3, maxError, summary);
//
//	LOG(INFO) << "Support is " << (summary.inliers.size() * 1.0 / nrPoints);

// save the transformations
	for (int ii = 0; ii < nrPoints; ii++) {
		transformations[idxs[ii]] = sim3;
	}

	// ok, we have a transformation for this quad now => test its support
	Eigen::Matrix<float, Eigen::Dynamic, 1> errors(nrPoints, 1);
	for (int ii = 0; ii < nrPoints; ii++) {
		// get the error
		errors[ii] = (newPoints[idxs[ii]] - (sim3 * oldPoints[idxs[ii]])).norm();
	}

	double mean = errors.mean();
	errors = errors.array() - errors.mean();
	double stdev = std::sqrt((errors.cwiseProduct(errors)).mean());
//
//	LOG(INFO)<< "calculated SE3 from " << nrPoints << " points";
//	LOG(INFO)<< "mean error: " << mean << "(normalized: "<< mean/maxError << ")";
//	LOG(INFO)<< "stddev: " << stdev << "(normalized: "<< stdev/maxError << ")";

	if (mean <= maxError /* &&  (summary.inliers.size() * 1.0 / nrPoints) > 0.95 */) {
		// this cell is finished => add the estimated transformation

		Eigen::Vector3f cl(Eigen::Vector3f::Zero());
		for (int ii = 0; ii < nrPoints; ii++)
			cl += oldPoints[idxs[ii]];
		cl /= nrPoints;

//		tCenters.push_back(cl);
//		tSims.push_back(sim3);

		return mean * nrPoints;
	}

	// otherwise we have to split this leaf and process all of them
	// now change the cell
	std::vector<std::shared_ptr<PointWithIdx>> oldCellPoints;
	Branch<std::shared_ptr<PointWithIdx>>* newBranch = leaf->split(oldCellPoints); // from here on, the cell* is not valid

	// => leaf is not valid anymore
	delete leaf;
	leaf = nullptr;

	// add the points to the new branch
	for (const auto& p : oldCellPoints)
		newBranch->at(p->x(), p->y(), p->z())->data.emplace_back(p);

//	LOG(INFO) << "branching";

	// process those cells
	mean = 0;
	for (int ii = 0; ii < 8; ii++) {
		Leaf<shared_ptr<PointWithIdx>>* l =
				reinterpret_cast<Leaf<shared_ptr<PointWithIdx>>*>(newBranch->children[ii]);
		mean += estimateCell(l, oldPoints, newPoints, transformations, maxError);
	}

	return mean;
}

// ###################################################################

#define PRINT(x) \
{ \
Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "[", "]", "np.array([", "])"); \
LOG(INFO) << #x << " = " << x.format(CommaInitFmt); \
}

int main(int argc, char* argv[]) {
	google::InitGoogleLogging(argv[0]);
	FLAGS_colorlogtostderr = true;
	FLAGS_logtostderr = true;
	GFLAGS_NAMESPACE::ParseCommandLineFlags(&argc, &argv, true);

	LOG(INFO)<< " =============================================================== ";
	LOG(INFO)<< " ======== welcome to the Progressive All The Way         ======= ";
	LOG(INFO)<< " =============================================================== ";


	// check the arguments
	CHECK(stlplus::file_readable(FLAGS_nvm1)) << "input file <" << FLAGS_nvm1 << "> not readable";
	CHECK(stlplus::file_readable(FLAGS_nvm2)) << "input file <" << FLAGS_nvm2 << "> not readable";
	CHECK(stlplus::file_readable(FLAGS_patches)) << "input file <" << FLAGS_patches
															<< "> not readable";
	CHECK(stlplus::folder_exists(FLAGS_outdir) || stlplus::folder_create(FLAGS_outdir))
																									<< "Unable to create output folder <"
																									<< FLAGS_outdir

																									<< ">";

	// create a timing file
	std::ofstream f_timing(stlplus::create_filespec(FLAGS_outdir, "timing", "txt"));
	f_timing << "# timing of the progressive 3d modelling all the way" << std::endl;
	f_timing << "# --nvm1=" << FLAGS_nvm1 << std::endl;
	f_timing << "# --nvm2=" << FLAGS_nvm2 << std::endl;
	f_timing << "# --patches=" << FLAGS_patches << std::endl;
	f_timing << "# --subtees=" << FLAGS_subtrees << std::endl;
	f_timing << "# cores: " << omp_get_max_threads() << std::endl;
	f_timing << "# ############################## " << std::endl;

	auto t0 = std::chrono::steady_clock::now();
	auto t1 = t0;
	double time_total = 0;

// load the nvm files
	std::vector<mo3d::NVM_Model> nvm1;
	mo3d::NVMReader::readFile(FLAGS_nvm1.c_str(), nvm1, true);
	CHECK(nvm1.size() > 0) << " no models found in <" << FLAGS_nvm1 << ">";

	std::vector<mo3d::NVM_Model> nvm2;
	mo3d::NVMReader::readFile(FLAGS_nvm2.c_str(), nvm2, true);
	CHECK(nvm2.size() > 0) << " no models found in <" << FLAGS_nvm2 << ">";

// --------------------------------------------------

	std::vector<theia::Camera> oldCameras;
	std::vector<theia::Camera> newCameras;
	std::map<int, int> oldCamIdx2New;

	for (int oldIdx = 0; oldIdx < nvm1[0].cameras.size(); oldIdx++) {
		const mo3d::NVM_Camera& c = nvm1[0].cameras[oldIdx];
		theia::Camera oldCam = nvmToTheia(c);
		theia::Camera newCam;

		// look for the camera in the new file
		bool found = false;
		for (int newIdx = 0; newIdx < nvm2[0].cameras.size(); newIdx++) {
			const mo3d::NVM_Camera& cn = nvm2[0].cameras[newIdx];
			if (stlplus::basename_part(cn.filename).compare(stlplus::basename_part(c.filename))
					== 0) {
				found = true;
				newCam = nvmToTheia(cn);
				oldCamIdx2New[oldIdx] = newIdx;
			}
		}

		CHECK(found) << c.filename << " not found in nvm2";

		oldCameras.emplace_back(oldCam);
		newCameras.emplace_back(newCam);
	}

	std::vector<int> newCamerasNewIdx;
	for (int newIdx = 0; newIdx < nvm2[0].cameras.size(); newIdx++) {
		const mo3d::NVM_Camera& cn = nvm2[0].cameras[newIdx];
		bool found = false;
		for (int oldIdx = 0; oldIdx < nvm1[0].cameras.size(); oldIdx++) {
			const mo3d::NVM_Camera& c = nvm1[0].cameras[oldIdx];

			if (stlplus::basename_part(cn.filename).compare(stlplus::basename_part(c.filename))
					== 0) {
				found = true;
			}
		}
		if (found == false)
			newCamerasNewIdx.push_back(newIdx);
	}

	LOG(INFO)<< "Cameras mapped, have " << newCamerasNewIdx.size() << " new cameras";

	//////////////////////////////////////////////////////////////////

	// load the patches
	std::vector<mo3d::Ppatch3d> patches;

	if (FLAGS_ismesh) {
		CHECK(loadPatchesMeshPly(oldCameras, FLAGS_patches, patches))
																				<< "unable to load ply file <"
																				<< FLAGS_patches
																				<< ">";
		;
	} else {
		CHECK(loadPatchesPly(FLAGS_patches.c_str(), patches)) << "unable to load ply file <"
																		<< FLAGS_patches << ">";

	}

	// --------------------------------------------------------

	// create a tree for the patches
	Eigen::Vector3f min, max, dist;
	getBoundingBox(patches, min, max);
	dist = max - min;
	float width = std::max(dist[0], std::max(dist[1], dist[2]));

	//////////////////////////////////////////////////////////////////

	if (false) {
		// project points into a debug image
		int ii = 0;
		for (const auto& c : oldCameras) {
			cimg_library::CImg<float> img(c.ImageWidth() / 4, c.ImageHeight() / 4, 1, 3);
			img.fill(0);

			for (auto& patch : patches) {
				Eigen::Vector2d pixel;
				double d = c.ProjectPoint(patch->center_.cast<double>(), &pixel);
				pixel /= 4.0;

				if (d < 0) {
					LOG(WARNING)<< " point behind camera";
					continue;
				}

				// set the pixel
				img.atXY((int) pixel[0], (int) pixel[1], 0) = patch->color_[0];
				img.atXY((int) pixel[0], (int) pixel[1], 1) = patch->color_[1];
				img.atXY((int) pixel[0], (int) pixel[1], 2) = patch->color_[2];

			}

			// save the image
			std::string imgName = stlplus::create_filespec("/tmp", "camera_" + std::to_string(ii++),
					".jpg");
			img.save(imgName.c_str());
		}
	}

	t0 = std::chrono::steady_clock::now();
	std::vector<Eigen::Vector3f> oldPoints;
	std::vector<Eigen::Vector3f> newPoints;
	std::vector<Eigen::Matrix<unsigned char, 3, 1>> colors;
	std::vector<float> newPointErrors;

	// now re-triangulate the points
	for (auto& patch : patches) {
		std::vector<Eigen::Vector3d> rays;
		std::vector<Eigen::Vector2d> points;
		std::vector<Eigen::Matrix<double, 3, 4> > poses;
		std::vector<Eigen::Vector3d> origins;
		for (auto imgIdx : patch->images_) {
			Eigen::Vector2d pixel;
			if (oldCameras[imgIdx].ProjectPoint(patch->center_.cast<double>(), &pixel) < 0) {
				LOG(ERROR)<< "negative depth in triangulation";
				continue;
			}

			rays.emplace_back(newCameras[imgIdx].PixelToUnitDepthRay(pixel).normalized());
			origins.emplace_back(newCameras[imgIdx].GetPosition());
			points.emplace_back(pixel);
			poses.emplace_back();
			newCameras[imgIdx].GetProjectionMatrix(&poses.back());
		}

		// now triangulate
		Eigen::Vector4d newPoint;
//		if (!theia::TriangulateMidpoint(origins, rays, &newPoint))
//			newPoint << Eigen::Vector4d::Zero();
		if (!theia::TriangulateNView(poses, points, &newPoint))
			newPoint << Eigen::Vector4d::Ones();
		newPoint.head(3) = newPoint.hnormalized();
		newPoint[3] = 1.0;

		// get the reprojection error
		float err = 0;
		Eigen::Vector2d pixel;
		for (int ii = 0; ii < patch->images_.size(); ii++) {
			newCameras[patch->images_[ii]].ProjectPoint(newPoint, &pixel);
			err += (pixel - points[ii]).norm();
		}
		err = err / (patch->images_.size());

		oldPoints.emplace_back(patch->center_.hnormalized());
		newPoints.emplace_back(newPoint.cast<float>().hnormalized());
		newPointErrors.emplace_back(err);
		colors.emplace_back(patch->color_.cast<unsigned char>());

//		patch->center_ = newPoint.cast<float>();
	}
	t1 = std::chrono::steady_clock::now();
	double t_triangulation = std::chrono::duration<double>(t1 - t0).count();
	f_timing << "triangulation " << t_triangulation << std::endl;
	time_total += t_triangulation;
	LOG(INFO)<< "points triangulated in " << t_triangulation << " seconds";

	// clean the points by the triangulation error
	std::vector<Eigen::Vector3f> oldPoints_cleared;
	std::vector<Eigen::Vector3f> newPoints_cleared;
	std::vector<Eigen::Matrix<unsigned char, 3, 1>> colors_cleared;
	std::vector<float> newPointErrors_cleared;
	std::vector<mo3d::Ppatch3d> patches_cleared;

	std::vector<int> newIdx2old;
	std::vector<bool> dirtyVertices_orig(oldPoints.size(), false);
	for (size_t ii = 0; ii < oldPoints.size(); ii++) {
		if (newPointErrors[ii] < FLAGS_maxtriang) {
			oldPoints_cleared.emplace_back(oldPoints[ii]);
			newPoints_cleared.emplace_back(newPoints[ii]);
			colors_cleared.emplace_back(colors[ii]);
			newPointErrors_cleared.emplace_back(newPointErrors[ii]);
			patches_cleared.emplace_back(patches[ii]);
			newIdx2old.emplace_back(ii);
		} else {
			dirtyVertices_orig[ii] = true;
		}
	}

	LOG(INFO)<< "Removed " << (oldPoints.size() - oldPoints_cleared.size()) << " due to high reprojection error";

	oldPoints.swap(oldPoints_cleared);
	newPoints.swap(newPoints_cleared);
	colors.swap(colors_cleared);
	newPointErrors.swap(newPointErrors_cleared);
	patches.swap(patches_cleared);

	// create a tree with the old points
	// #####################################################################################
	DynOctTree<std::shared_ptr<PointWithIdx> > oldPtTree((min + max) / 2.0, width);
	const int initLevel = 1;
	for (size_t ii = 0; ii < oldPoints.size(); ii++) {
		oldPtTree.add(std::shared_ptr<PointWithIdx>(new PointWithIdx(oldPoints[ii], ii)),
				width / (1 << initLevel));
	}
	oldPtTree.cellHistogram();

	// #####################################################################################
	t0 = std::chrono::steady_clock::now();
//	std::unique_ptr<fastann::nn_obj<float> > searchTree(
//			fastann::nn_obj_build_exact(oldPoints[0].data(), oldPoints.size(), 3));

	std::vector<Eigen::Vector3f> transformedPoints(oldPoints.size());
	std::vector<Sophus::Sim3f> applied_transformations(oldPoints.size());

	if (!FLAGS_skip_estimation) {

		// loop through the nodes of the tree
		const float maxError = width * FLAGS_maxdist;
		// collect the leafs
		Leaf_iterator<std::shared_ptr<PointWithIdx>> end = oldPtTree.end();
		std::vector<Leaf<std::shared_ptr<PointWithIdx>>*> leafs;
		for (Leaf_iterator<std::shared_ptr<PointWithIdx>> leaf = oldPtTree.begin(); leaf != end;
				leaf++) {
			leafs.emplace_back(&(*leaf));
		}

		// process the leafs
		float weighted_mean = 0;
		std::vector<Sophus::Sim3f> transformations(oldPoints.size());
		for (int ii = 0; ii < leafs.size(); ii++) {
			weighted_mean += estimateCell(leafs[ii], oldPoints, newPoints, transformations,
					maxError);
		}

		LOG(INFO)<< "Total mean error: " << (weighted_mean / oldPoints.size()) << " [max = " << maxError << "]";
		oldPtTree.cellHistogram();

		LOG(INFO)<< "Boost => " << (oldPtTree.getRoot()->nrLeafs() * 1.0 / oldPoints.size());

		// re-transform new points using the new estimated transformations
		for (size_t ii = 0; ii < oldPoints.size(); ii++) {
			transformedPoints[ii] = transformations[ii] * oldPoints[ii];
			applied_transformations[ii] = transformations[ii];
		}

		/*		// 2nd option => use interpolation to blur transformations on cell corners
		 for (size_t ii = 0; ii < oldPoints.size(); ii++) {
		 const int nrNeighs = 4;
		 Eigen::VectorXf dists(nrNeighs);
		 Eigen::Matrix<unsigned int, Eigen::Dynamic, 1> idxs(nrNeighs, 1);
		 searchTree->search_knn(oldPoints[ii].data(), 1, nrNeighs, idxs.data(), dists.data());

		 // Collect the transformations
		 Eigen::Matrix<float, Eigen::Dynamic, 7> localTransforms(nrNeighs, 7);
		 for (int jj = 0; jj < nrNeighs; jj++) {
		 localTransforms.row(jj) = transformations[idxs[jj]].log().transpose();
		 }

		 Eigen::Matrix<float, 7, 1> meanSimV = localTransforms.colwise().mean();

		 transformedPoints[ii] = Sophus::Sim3f::exp(meanSimV) * oldPoints[ii];
		 applied_transformations[ii] = Sophus::Sim3f::exp(meanSimV);
		 }
		 */

	} else {
		// we skipped recursive transformation estimation
		for (size_t ii = 0; ii < oldPoints.size(); ii++) {
			transformedPoints[ii] = newPoints[ii];
		}

	}

	// update timing
	t1 = std::chrono::steady_clock::now();
	double t_transformation = std::chrono::duration<double>(t1 - t0).count();
	f_timing << "transformation " << t_transformation << std::endl;
	time_total += t_transformation;

	// ######################################################################
	// nearest neighbors consistency test
	t0 = std::chrono::steady_clock::now();
//	std::unique_ptr<fastann::nn_obj<float> > searchTree_after(
//			fastann::nn_obj_build_exact(transformedPoints[0].data(), transformedPoints.size(), 3));

	// ######################################################################

	// create an octree with the points before
	std::vector<std::shared_ptr<PointWithIdx>> ptsIdxBefore(oldPoints.size());
	for (size_t ii = 0; ii < ptsIdxBefore.size(); ii++)
		ptsIdxBefore[ii] = std::shared_ptr<PointWithIdx>(new PointWithIdx(oldPoints[ii], ii));
	getBoundingBox(ptsIdxBefore, min, max);
	dist = max - min;
	width = std::max(dist[0], std::max(dist[1], dist[2]));
	DynOctTree<std::shared_ptr<PointWithIdx> > beforeOctree((min + max) / 2.0, width);
	for (size_t ii = 0; ii < ptsIdxBefore.size(); ii++)
		beforeOctree.addSingleCell(ptsIdxBefore[ii]);
	LOG(INFO)<< "Octree before: ";
	beforeOctree.cellHistogram();

	// create an octree with the points after
	std::vector<std::shared_ptr<PointWithIdx>> ptsIdxAfter(newPoints.size());
	for (size_t ii = 0; ii < ptsIdxAfter.size(); ii++)
		ptsIdxAfter[ii] = std::shared_ptr<PointWithIdx>(
				new PointWithIdx(transformedPoints[ii], ii));
	getBoundingBox(ptsIdxAfter, min, max);
	dist = max - min;
	width = std::max(dist[0], std::max(dist[1], dist[2]));
	DynOctTree<std::shared_ptr<PointWithIdx> > afterOctree((min + max) / 2.0, width);
	for (size_t ii = 0; ii < ptsIdxAfter.size(); ii++)
		afterOctree.addSingleCell(ptsIdxAfter[ii]);
	LOG(INFO)<< "Octree after: ";
	afterOctree.cellHistogram();

	// create pcl octrees
	pcl::PointCloud<pcl::PointXYZ>::Ptr pcl_cloud_before(new pcl::PointCloud<pcl::PointXYZ>);
	pcl_cloud_before->width = 1;
	pcl_cloud_before->height = oldPoints.size();
	pcl_cloud_before->points.resize(oldPoints.size());
	for (int ii = 0; ii < oldPoints.size(); ii++) {
		pcl_cloud_before->points[ii].x = oldPoints[ii][0];
		pcl_cloud_before->points[ii].y = oldPoints[ii][1];
		pcl_cloud_before->points[ii].z = oldPoints[ii][2];
	}
	pcl::octree::OctreePointCloudSearch<pcl::PointXYZ> pcl_beforeOctree(width / (1 << 8));
	pcl_beforeOctree.setInputCloud(pcl_cloud_before);
	pcl_beforeOctree.addPointsFromInputCloud();

	pcl::PointCloud<pcl::PointXYZ>::Ptr pcl_cloud_after(new pcl::PointCloud<pcl::PointXYZ>);
	pcl_cloud_after->width = 1;
	pcl_cloud_after->height = transformedPoints.size();
	pcl_cloud_after->points.resize(transformedPoints.size());
	for (int ii = 0; ii < transformedPoints.size(); ii++) {
		pcl_cloud_after->points[ii].x = transformedPoints[ii][0];
		pcl_cloud_after->points[ii].y = transformedPoints[ii][1];
		pcl_cloud_after->points[ii].z = transformedPoints[ii][2];
	}
	pcl::octree::OctreePointCloudSearch<pcl::PointXYZ> pcl_afterOctree(width / (1 << 8));
	pcl_afterOctree.setInputCloud(pcl_cloud_after);
	pcl_afterOctree.addPointsFromInputCloud();

	t1 = std::chrono::steady_clock::now();
	double t_octrees = std::chrono::duration<double>(t1 - t0).count();
	LOG(INFO)<< "built trees in " << t_octrees << " seconds";

	// check for inconsistency
	std::vector<float> consistency(transformedPoints.size());
	unsigned int consistentPoints = 0;
	int wrongNN = 0;
	std::atomic<int> ptCnt(0);
#pragma omp parallel for
	for (size_t ii = 0; ii < transformedPoints.size(); ii++) {
//		LOG(INFO) << "nn for " << ii << " / " << transformedPoints.size();
		const int nrNeighs = FLAGS_nn;

		ptCnt++;
		int myCnt = ptCnt.load();
		if (myCnt % (transformedPoints.size() / 20) == 0)
			LOG(INFO)<< std::round(myCnt * 100.0 / transformedPoints.size()) << "%";

//		std::vector<float> dists(nrNeighs);
//		Eigen::Matrix<unsigned int, Eigen::Dynamic, 1> idxs(nrNeighs, 1);
//		// get the 4 nearest neighbours to the point
//		searchTree->search_knn(oldPoints[ii].data(), 1, nrNeighs, idxs.data(), dists.data());
//
////		for (int jj = 0; jj < neighs.size(); jj++)
////			LOG(INFO) << "FLANN neigh " << jj  << " dist " << dists[jj] << " idx " << idxs[jj] << " pt [" << oldPoints[idxs[jj]].transpose() << "]";
//
//		std::vector<float> dists_after(nrNeighs);
//		Eigen::Matrix<unsigned int, Eigen::Dynamic, 1> idxs_after(nrNeighs, 1);
//		searchTree_after->search_knn(transformedPoints[ii].data(), 1, nrNeighs, idxs_after.data(),
//				dists_after.data());

			// octree nn
//		std::vector<std::shared_ptr<PointWithIdx> > neighsO;
//		kNearestNeighbour(beforeOctree, oldPoints[ii], nrNeighs, neighsO);
//		std::vector<std::shared_ptr<PointWithIdx> > neighsA;
//		kNearestNeighbour(afterOctree, transformedPoints[ii], nrNeighs, neighsA);

			// pcl nn
		pcl::PointXYZ searchPoint1;
		searchPoint1.x = oldPoints[ii][0];
		searchPoint1.y = oldPoints[ii][1];
		searchPoint1.z = oldPoints[ii][2];
		std::vector<int> pointIdxNKNSearchBefore;
		std::vector<float> pointNKNSquaredDistanceBefore;
		CHECK(pcl_beforeOctree.nearestKSearch (searchPoint1, nrNeighs, pointIdxNKNSearchBefore, pointNKNSquaredDistanceBefore) > 0)
																																				<< "nn search failed";
		pcl::PointXYZ searchPoint2;
		searchPoint2.x = transformedPoints[ii][0];
		searchPoint2.y = transformedPoints[ii][1];
		searchPoint2.z = transformedPoints[ii][2];
		std::vector<int> pointIdxNKNSearchAfter;
		std::vector<float> pointNKNSquaredDistanceAfter;
		CHECK(pcl_afterOctree.nearestKSearch (searchPoint2, nrNeighs, pointIdxNKNSearchAfter, pointNKNSquaredDistanceAfter) > 0)
																																			<< " nn search failed";

//		for (int jj = 0; jj < nrNeighs; jj++)
////			CHECK_EQ(neighs[jj]->idx_, idxs[jj]) << "nearest neighbours are not same (ii: " << ii << ", jj: " << jj << ")";
//			wrongNN += (pointIdxNKNSearchBefore[jj] != idxs[jj]);
////			wrongNN += (neighsO[jj]->idx_ != idxs[jj]);
//		for (int jj = 0; jj < nrNeighs; jj++)
////			CHECK_EQ(neighs[jj]->idx_, idxs[jj]) << "nearest neighbours are not same (ii: " << ii << ", jj: " << jj << ")";
////			wrongNN += (neighsA[jj]->idx_ != idxs_after[jj]);
//			wrongNN += (pointIdxNKNSearchAfter[jj] != idxs_after[jj]);

		int changedNeighs = 0;
		for (int jj = 0; jj < nrNeighs; jj++)
			changedNeighs += (pointIdxNKNSearchBefore[jj] != pointIdxNKNSearchAfter[jj]);
//			changedNeighs += (neighsA[jj]->idx_ != neighsO[jj]->idx_);
//			changedNeighs += (idxs[jj] != idxs_after[jj]);

		consistency[ii] = changedNeighs * 1.0 / nrNeighs;

//		float changedNeighs = idxs.cwiseNotEqual(idxs_after).count() * 1.0 / nrNeighs;
//		consistency.push_back(changedNeighs);

//		// get the transformations of the neighbours
//		Eigen::Matrix<double, Eigen::Dynamic, 6> t(nrNeighs,6);
//		for (size_t jj = 0; jj < nrNeighs; jj++)
//			t.row(jj) << transformations[idxs[jj]].log().transpose().cast<double>();
//
//
//		Eigen::Matrix<double,6,1> avg = t.colwise().sum().transpose()/nrNeighs;
//		Eigen::Matrix<double,6,1> avg2 = (t.cwiseProduct(t)).colwise().sum().transpose() / nrNeighs;
//		Eigen::Matrix<double,6,1> stdev = (avg2 - avg.cwiseProduct(avg)).cwiseSqrt();
//
//		for (int j = 0; j < 6 ; j++)
//			stdev[j] = (tmean[j] > 1e-5) ? stdev[j] / tmean[j] : 0;
//
//		consistency.push_back(stdev.norm());
	}
	LOG(INFO)<< wrongNN << " out of " << (oldPoints.size() * FLAGS_nn * 2) << " were wrong";

	// update timing
	t1 = std::chrono::steady_clock::now();
	double t_consistency = std::chrono::duration<double>(t1 - t0).count();
	f_timing << "consistency " << t_consistency << std::endl;
	time_total += t_consistency;

	LOG(INFO)<< "consistency estimated";

//	LOG(INFO)<< (consistentPoints * 1.0 / oldPoints.size()) << " consistent points";

	// create a ply file for the new points
	yaply::PlyFile ply;
	ply["vertex"].nrElements = newPoints.size();
	ply["vertex"].setScalars("x,y,z", newPoints[0].data());
	ply["vertex"].setScalars("red,green,blue", colors[0].data());
	ply["vertex"].setScalars("scalar_consistency", consistency.data());
	ply["vertex"].setScalars("scalar_reprojection", newPointErrors.data());
	ply.save(stlplus::create_filespec(FLAGS_outdir, "newPoints", "ply").c_str(), false);

	ply["vertex"].setScalars("x,y,z", oldPoints[0].data());
	ply.save(stlplus::create_filespec(FLAGS_outdir, "oldPoints", "ply").c_str(), false);

	ply["vertex"].setScalars("x,y,z", transformedPoints[0].data());
	ply.save(stlplus::create_filespec(FLAGS_outdir, "transformedPoints", "ply").c_str(), false);

	// clean some stuff
	newPoints.clear();
	colors.clear();
	newPointErrors.clear();
	oldPoints.clear();

	// clear the trees
	pcl_afterOctree.deleteTree();
	pcl_beforeOctree.deleteTree();
	pcl_cloud_after->clear();
	pcl_cloud_before->clear();


	// ####################################################################
	// now start hpmvs
	LOG(INFO)<< " #################################################### ";
	LOG(INFO)<< " ###########    Start Progressive MVS    ############ ";
	LOG(INFO)<< " #################################################### ";

	mo3d::HpmvsOptions options;
	options.OUTFOLDER = FLAGS_outdir;
	mo3d::Scene scene;
	scene.addCameras(nvm2[0], options);
	scene.extractCoVisiblilty(nvm2[0], options);

	// create optimizers
	const int nrThreads = omp_get_max_threads();

	std::vector<mo3d::PatchOptimizer> optimizers;
	for (int ii = 0; ii < nrThreads; ii++)
		optimizers.emplace_back(options, &scene);

	// get the dirty patches
	t0 = std::chrono::steady_clock::now();
	for (int ii = 0; ii < patches.size(); ii++) {
		bool dirty = false;
		// first we have to convert the patch
		mo3d::Patch3d& patch = *patches[ii];
		patch.center_ = transformedPoints[ii].homogeneous();
		Eigen::Vector3f n = patch.normal_.hnormalized();
		patch.normal_ = (applied_transformations[ii].rxso3() * n).normalized().homogeneous();

		// run through the visible images
		std::vector<int> newVImages;
		for (int jj = 0; jj < patch.images_.size(); jj++) {
			int oldIdx = patch.images_[jj];
			if (oldCamIdx2New.find(oldIdx) == oldCamIdx2New.end()) { // camera was removed
				if (jj == 0) // it was the reference idx
					dirty = true;
			} else {
				newVImages.push_back(oldCamIdx2New[oldIdx]);
			}
		}
		patch.images_.swap(newVImages);

		// if we have to less images, we mark as dirty as well
		if (patch.images_.size() < options.MIN_IMAGES_PER_PATCH)
			dirty = true;

		// check consistency
		if (!dirty && consistency[ii] > 0.5)
			dirty = true;

		// check the visibility of the patch in the new cameras
		if (!FLAGS_skip_dirtyfromvisible && !dirty && patch.images_.size() <= options.MAX_IMAGES_PER_PATCH - 1) {
			for (int camIdx : newCamerasNewIdx) {
				Eigen::Vector3f p = scene.cameras_[camIdx].project(patch.center_, 0);
				if (p[0] > 0 && p[0] < scene.images_[camIdx].getWidth() && p[1] > 0
						&& p[1] < scene.images_[camIdx].getHeight()) {
					// TODO also check normal!
					dirty = true;
					break;
				}
			}
		}

		patch.dirty_ = dirty;
		patch.expanded_ = !dirty; // assume expanded if not dirty
		patch.flatness_ = (dirty) ? -1.0 : 0;
	}

	// clear some stuff
	transformedPoints.clear();
	consistency.clear();

	// update timing
	t1 = std::chrono::steady_clock::now();
	double t_dirtymask = std::chrono::duration<double>(t1 - t0).count();
	f_timing << "dirtymask " << t_dirtymask << std::endl;
	time_total += t_dirtymask;

	int nrDirty = 0;
	for (const auto& p : patches)
		nrDirty += (p->dirty_ == true);
	LOG(INFO)<< nrDirty << " / " << patches.size() << " are dirty";

	// check the dirty patches
	t0 = std::chrono::steady_clock::now();
	patches_cleared.clear();
//#pragma omp parallel for
	for (int ii = 0; ii < patches.size(); ii++) {
		int originalIdx = newIdx2old[ii];
		if (!patches[ii]->dirty_) {
			patches_cleared.emplace_back(patches[ii]);
			continue;
		}

		Eigen::Vector3f posBefore = patches[ii]->center_.hnormalized();

		// optimization
		const int threadId = omp_get_thread_num();
		if (!optimizers[threadId].optimize(*patches[ii])) {
			dirtyVertices_orig[originalIdx] = true;
			continue;
		}

		// make sure we did not moved to far...
		if ((patches[ii]->center_.hnormalized() - posBefore).norm() > patches[ii]->scale_3dx_ * 2) {
			dirtyVertices_orig[originalIdx] = true;
			continue;
		}

		patches_cleared.emplace_back(patches[ii]);
	}
	LOG(INFO)<< patches_cleared.size() << " / " << patches.size() << " are valid after optimization";
	patches.swap(patches_cleared);
	patches_cleared.clear();

	// output a list of dirty vertices

	std::ofstream f_dirtyvertices(
			stlplus::create_filespec(FLAGS_outdir,
					stlplus::basename_part(FLAGS_patches) + "_dirtyVertices", "txt"));
	for (int ii = 0; ii < dirtyVertices_orig.size(); ii++)
		if (dirtyVertices_orig[ii])
			f_dirtyvertices << ii << std::endl;
	f_dirtyvertices.close();

	if (FLAGS_ismesh) {
		saveColoredMesh(FLAGS_patches, dirtyVertices_orig);
	}

	// update timing
	t1 = std::chrono::steady_clock::now();
	double t_dirtyoptim = std::chrono::duration<double>(t1 - t0).count();
	f_timing << "dirtyoptim " << t_dirtyoptim << std::endl;
	time_total += t_dirtyoptim;

	// create a tree for the new patches
//	Eigen::Vector3f min, max, dist;
	getBoundingBox(patches, min, max);
	for (const auto& p : nvm2[0].points) {
		min[0] = std::min(min[0], (float) p.xyz[0]);
		min[1] = std::min(min[1], (float) p.xyz[1]);
		min[2] = std::min(min[2], (float) p.xyz[2]);

		max[0] = std::max(max[0], (float) p.xyz[0]);
		max[1] = std::max(max[1], (float) p.xyz[1]);
		max[2] = std::max(max[2], (float) p.xyz[2]);
	}

	dist = max - min;
	width = std::max(dist[0], std::max(dist[1], dist[2]));

	Branch<mo3d::Ppatch3d>* oldRoot_p = scene.patchTree_.swapRoot(
			new Branch<mo3d::Ppatch3d>((min + max) / 2.0, width));
	delete oldRoot_p;

	// and add them
	for (auto& p : patches) {
		scene.patchTree_.add(p, p->scale_3dx_);
		scene.setDepths(*p, false);
	}
	scene.patchTree_.cellHistogram();

	// check the new points from the nvm2
	int newFromSfM = 0;
	for (int ii = 0; ii < nvm2[0].points.size(); ii++) {
		Eigen::Vector3f pos = nvm2[0].points[ii].xyz.cast<float>();
		Leaf<mo3d::Ppatch3d>* leaf = scene.patchTree_.at(pos[0], pos[1], pos[2]);
		if (scene.patchTree_.nodeLevel(leaf) > 9 || !leaf->data.empty())
			continue;

		const mo3d::NVM_Point& pt = nvm2[0].points[ii];
		mo3d::Ppatch3d ppatch(new mo3d::Patch3d);
		ppatch->center_ = pt.xyz.cast<float>().homogeneous();

		// only accept point with enough images
		if (pt.measurements.size() < options.MIN_IMAGES_PER_PATCH)
			continue;

		Eigen::Vector4f commonCenter = Eigen::Vector4f::Zero();

		// add all measurements to the patch
		for (const mo3d::NVM_Measurement& m : pt.measurements) {
//			int idx = dict_[m.imgIndex];
			int idx = m.imgIndex;
			if (idx < 0)
				continue;

			// check visibility of that patch:
			Eigen::Vector3f p = scene.cameras_[idx].project(ppatch->center_, options.START_LEVEL);

			const int margin = 2;
			if (p[0] < margin || p[1] < margin
					|| p[0] >= scene.images_[idx].getWidth(options.START_LEVEL) - margin
					|| p[1] >= scene.images_[idx].getHeight(options.START_LEVEL) - margin) {
				continue;
			}

			ppatch->images_.push_back(idx);
			commonCenter += scene.cameras_[idx].center_;
		}

		if (ppatch->images_.size() < 2)
			continue; // not enough images attached

		commonCenter = commonCenter / (float) ppatch->images_.size();
		ppatch->normal_ = commonCenter - ppatch->center_;
		ppatch->normal_ = scene.cameras_[ppatch->images_[0]].center_ - ppatch->center_; // just first
		ppatch->normal_.normalize();
		ppatch->normal_[3] = 0.0;
		ppatch->scale_3dx_ = scene.cameras_[ppatch->images_[0]].getScale(ppatch->center_,
				options.START_LEVEL);
		ppatch->scale_3dx_ = std::max(ppatch->scale_3dx_, width / (1 << 10));

		// optimization
		const int threadId = omp_get_thread_num();
		if (!optimizers[threadId].optimize(*ppatch))
			continue;

		// make sure we did not moved to far...
		if ((ppatch->center_.head(3) - pt.xyz.cast<float>()).norm() > ppatch->scale_3dx_ * 2)
			continue;

		// also do some additional tests
		bool patchGood = (scene.patchTree_.at(ppatch->x(), ppatch->y(), ppatch->z()) == leaf);
		patchGood = patchGood
				&& scene.depthTests(*ppatch, options.DEPTH_TEST_FACTOR)
						>= options.MIN_IMAGES_PER_PATCH
				&& scene.viewBlockTest(*ppatch, options.DEPTH_TEST_FACTOR)
						< options.MIN_IMAGES_PER_PATCH;

		// add the patch
		scene.patchTree_.add(ppatch, ppatch->scale_3dx_);
		scene.setDepths(*ppatch, false);
		newFromSfM++;
	}
	LOG(INFO)<< "added " << newFromSfM << " points from sfm to the dense model";
	scene.patchTree_.cellHistogram();

	// clear nvm models
	nvm1.clear();
	nvm2.clear();

	// divide the initial tree into multiple subtrees
	std::vector<std::shared_ptr<DynOctTree<mo3d::Ppatch3d> > > subTrees;
	getSubTrees(scene.patchTree_, subTrees, FLAGS_subtrees);

	// initialize the CellProcessors
	std::vector<std::unique_ptr<mo3d::CellProcessor> > cellProcessors;
	std::function<void(mo3d::Ppatch3d, const float)> borderPatchFn = std::bind(
			&mo3d::CellProcessor::distributeBorderCell, &cellProcessors, std::placeholders::_1,
			std::placeholders::_2);
	for (int ii = 0; ii < subTrees.size(); ii++) {
		cellProcessors.emplace_back(new mo3d::CellProcessor(&scene, options));
		cellProcessors.back()->initFromTree(subTrees[ii].get(), &borderPatchFn, FLAGS_skip_clean);
	}

	// some  timing
	const auto start = std::chrono::steady_clock::now();
	double procTime;

	// process the queues level by level
	const int maxPrio = 200;
	for (int prio = 0; prio < maxPrio; prio += 1) {

		std::atomic<uint32_t> trees_changed(0);

#pragma omp parallel for schedule(dynamic)
		for (int ii = 0; ii < subTrees.size(); ii++) {
			const int threadId = omp_get_thread_num();
			if (cellProcessors[ii]->processQueue(&optimizers[threadId], prio))
				trees_changed++;
		}

		if (trees_changed.load() > 0 && ((int) prio) % 10 < 3) {
			scene.patchTree_.toExtPly(
					stlplus::create_filespec(options.OUTFOLDER, "patches-" + std::to_string(prio),
							"ply").c_str(), true);
			procTime =
					std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
			f_timing << "L" << prio << " " << scene.patchTree_.getRoot()->nrLeafs() << " "
					<< procTime << std::endl;

			LOG(INFO)<< "prio " << prio << " finished => saved to <" << stlplus::create_filespec(options.OUTFOLDER, "patches-" + std::to_string(prio),
					"ply") << ">";

			if (prio >= FLAGS_maxprio)
				break;

		}

		// are we finished ?
		bool moreWork = false;
		for (int ii = 0; ii < cellProcessors.size(); ii++)
			moreWork |= cellProcessors[ii]->haveWork();

		if (moreWork == false)
			break;
	}
	scene.patchTree_.cellHistogram();
	const auto end = std::chrono::steady_clock::now();
	procTime = std::chrono::duration<double>(end - start).count();
	LOG(INFO)<< "Done within " << procTime << " seconds";

	// update timing
	f_timing << "hpmvs " << procTime << std::endl;
	time_total += procTime;
	f_timing << "total " << time_total << std::endl;

	scene.patchTree_.toExtPly(
			stlplus::create_filespec(options.OUTFOLDER,
					"patches-final", "ply").c_str());

}

