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

#ifndef INCLUDE_HPMVS_ESTIMATESE3_HPP_
#define INCLUDE_HPMVS_ESTIMATESE3_HPP_

#include <Eigen/Dense>
#include <sophus/se3.hpp>
#include <sophus/sim3.hpp>
#include <glog/logging.h>

#include <theia/solvers/estimator.h>
#include <theia/sfm/create_and_initialize_ransac_variant.h>

namespace mo3d {

template<typename SCALAR = float>
bool Sim3fromCorrspondence(std::vector<Eigen::Matrix<SCALAR, 3, 1>>& pl,
		std::vector<Eigen::Matrix<SCALAR, 3, 1>>& pr, Sophus::Sim3Group<SCALAR>& sim3) {

	Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "[", "]",
			"np.array([", "])");

	// calculate the rigid transformation out of at least 3 correspondences
	// following the approach presented by Horn 1987

	if (pl.size() < 3 || pl.size() != pr.size())
		return false;

	Eigen::MatrixX3f ml(pl.size(), 3);
	Eigen::MatrixX3f mr(pr.size(), 3);
	for (int ii = 0; ii < pl.size(); ii++) {
		ml.row(ii) << pl[ii].transpose();
		mr.row(ii) << pr[ii].transpose();
	}

	VLOG(3) << "pl = " << ml.format(CommaInitFmt);
	VLOG(3) << "pr = " << mr.format(CommaInitFmt);

	// get the centroids of the points
	Eigen::Matrix<SCALAR, 3, 1> cl(Eigen::Matrix<SCALAR, 3, 1>::Zero());
	for (const auto& p : pl)
		cl += p;
	cl /= pl.size();

	Eigen::Matrix<SCALAR, 3, 1> cr(Eigen::Matrix<SCALAR, 3, 1>::Zero());
	for (const auto& p : pr)
		cr += p;
	cr /= pr.size();

	VLOG(3) << "cl = " << cl.format(CommaInitFmt);
	VLOG(3) << "cr = " << cr.format(CommaInitFmt);

	// get the scale
	double sl(0), sr(0);
	for (size_t ii = 0; ii < pl.size(); ii++) {
		sl += (pl[ii] - cl).squaredNorm();
		sr += (pr[ii] - cr).squaredNorm();
	}
	SCALAR scale = std::sqrt(sr / sl);

	// get help matrix
	Eigen::Matrix<SCALAR, 3, 3> m(Eigen::Matrix<SCALAR, 3, 3>::Zero());
	for (size_t ii = 0; ii < pl.size(); ii++)
		m += (pl[ii] - cl) * (pr[ii] - cr).transpose();

	VLOG(3) << "M = " << m.format(CommaInitFmt);

	// get 2nd help matrix
	Eigen::Matrix<SCALAR, 4, 4> n(Eigen::Matrix<SCALAR, 4, 4>::Zero());

	n(0, 0) = m(0, 0) + m(1, 1) + m(2, 2); // Sxx + Syy + Szz
	n(0, 1) = m(1, 2) - m(2, 1); // Syz - Szy
	n(0, 2) = m(2, 0) - m(0, 2); // Szx - Sxz
	n(0, 3) = m(0, 1) - m(1, 0); // Sxy - Syx

	n(1, 1) = m(0, 0) - m(1, 1) - m(2, 2); // Sxx - Syy - Szz
	n(1, 2) = m(0, 1) + m(1, 0); // Sxy + Syx
	n(1, 3) = m(2, 0) + m(0, 2); // Szx + Sxz

	n(2, 2) = -m(0, 0) + m(1, 1) - m(2, 2); // -Sxx + Syy - Szz
	n(2, 3) = m(1, 2) + m(2, 1); // Syz + Szy
	n(3, 3) = -m(0, 0) - m(1, 1) + m(2, 2); // -Sxx - Syy + Szz

	// update the lower part of the triangular matrix
	n(1, 0) = n(0, 1);
	n(2, 0) = n(0, 2);
	n(2, 1) = n(1, 2);
	n(3, 0) = n(0, 3);
	n(3, 1) = n(1, 3);
	n(3, 2) = n(2, 3);

	VLOG(3) << "N = " << n.format(CommaInitFmt);

	// now get the eigenvalues & -vectors
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix<SCALAR, 4, 4>> es(n);

	// get the largest positive eigenvalue
	if (es.eigenvalues()[3] <= 0 || es.info() != Eigen::Success)
		return false;

	VLOG(3) << "EigenVectors = " << es.eigenvectors().format(CommaInitFmt);
	VLOG(3) << "EigenValues  = " << es.eigenvalues().format(CommaInitFmt);

	// now we have a quaternion for the rotation
	Eigen::Matrix<SCALAR, 4, 1> vq = es.eigenvectors().col(3);
	vq *= scale;
	Eigen::Quaternion<SCALAR> qr(vq[0], vq[1], vq[2], vq[3]);
	sim3.rxso3().quaternion() = qr;

//	// get the error
//	std::vector<SCALAR> e(pl.size(),0);
//	SCALAR sum = 0;
//	for (int ii = 0; ii < pl.size(); ii++){
//		e[ii] = ((pr[ii] - cr ) - (sim3.so3() * (pl[ii] - cl))).norm();
//		e[ii] = e[ii] * e[ii];
//		sum += e[ii];
//	}
//	sum /= pl.size();
//	VLOG(3) << "mean residual = " << sum;

	// get the translation
	sim3.translation() = cr - sim3.rxso3() * cl;

	return true;
}

template<typename SCALAR = float>
bool SE3fromCorrspondence(std::vector<Eigen::Matrix<SCALAR, 3, 1>>& pl,
		std::vector<Eigen::Matrix<SCALAR, 3, 1>>& pr, Sophus::SE3Group<SCALAR>& se3) {

	Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "[", "]",
			"np.array([", "])");

	// calculate the rigid transformation out of at least 3 correspondences
	// following the approach presented by Horn 1987

	if (pl.size() < 3 || pl.size() != pr.size())
		return false;

	Eigen::MatrixX3f ml(pl.size(), 3);
	Eigen::MatrixX3f mr(pr.size(), 3);
	for (int ii = 0; ii < pl.size(); ii++) {
		ml.row(ii) << pl[ii].transpose();
		mr.row(ii) << pr[ii].transpose();
	}

	VLOG(3) << "pl = " << ml.format(CommaInitFmt);
	VLOG(3) << "pr = " << mr.format(CommaInitFmt);

	// get the centroids of the points
	Eigen::Matrix<SCALAR, 3, 1> cl(Eigen::Matrix<SCALAR, 3, 1>::Zero());
	for (const auto& p : pl)
		cl += p;
	cl /= pl.size();

	Eigen::Matrix<SCALAR, 3, 1> cr(Eigen::Matrix<SCALAR, 3, 1>::Zero());
	for (const auto& p : pr)
		cr += p;
	cr /= pr.size();

	VLOG(3) << "cl = " << cl.format(CommaInitFmt);
	VLOG(3) << "cr = " << cr.format(CommaInitFmt);

	// get help matrix
	Eigen::Matrix<SCALAR, 3, 3> m(Eigen::Matrix<SCALAR, 3, 3>::Zero());
	for (size_t ii = 0; ii < pl.size(); ii++)
		m += (pl[ii] - cl) * (pr[ii] - cr).transpose();

	VLOG(3) << "M = " << m.format(CommaInitFmt);

	// get 2nd help matrix
	Eigen::Matrix<SCALAR, 4, 4> n(Eigen::Matrix<SCALAR, 4, 4>::Zero());

	n(0, 0) = m(0, 0) + m(1, 1) + m(2, 2); // Sxx + Syy + Szz
	n(0, 1) = m(1, 2) - m(2, 1); // Syz - Szy
	n(0, 2) = m(2, 0) - m(0, 2); // Szx - Sxz
	n(0, 3) = m(0, 1) - m(1, 0); // Sxy - Syx

	n(1, 1) = m(0, 0) - m(1, 1) - m(2, 2); // Sxx - Syy - Szz
	n(1, 2) = m(0, 1) + m(1, 0); // Sxy + Syx
	n(1, 3) = m(2, 0) + m(0, 2); // Szx + Sxz

	n(2, 2) = -m(0, 0) + m(1, 1) - m(2, 2); // -Sxx + Syy - Szz
	n(2, 3) = m(1, 2) + m(2, 1); // Syz + Szy
	n(3, 3) = -m(0, 0) - m(1, 1) + m(2, 2); // -Sxx - Syy + Szz

	// update the lower part of the triangular matrix
	n(1, 0) = n(0, 1);
	n(2, 0) = n(0, 2);
	n(2, 1) = n(1, 2);
	n(3, 0) = n(0, 3);
	n(3, 1) = n(1, 3);
	n(3, 2) = n(2, 3);

	VLOG(3) << "N = " << n.format(CommaInitFmt);

	// now get the eigenvalues & -vectors
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix<SCALAR, 4, 4>> es(n);

	// get the largest positive eigenvalue
	if (es.eigenvalues()[3] <= 0 || es.info() != Eigen::Success)
		return false;

	VLOG(3) << "EigenVectors = " << es.eigenvectors().format(CommaInitFmt);
	VLOG(3) << "EigenValues  = " << es.eigenvalues().format(CommaInitFmt);

	// now we have a quaternion for the rotation
	auto& vq = es.eigenvectors().col(3);
	Eigen::Quaternion<SCALAR> qr(vq[0], vq[1], vq[2], vq[3]);
	se3.so3().setQuaternion(qr);

	// get the error
	std::vector<SCALAR> e(pl.size(), 0);
	SCALAR sum = 0;
	for (int ii = 0; ii < pl.size(); ii++) {
		e[ii] = ((pr[ii] - cr) - (se3.so3() * (pl[ii] - cl))).norm();
		e[ii] = e[ii] * e[ii];
		sum += e[ii];
	}
	sum /= pl.size();
	VLOG(3) << "mean residual = " << sum;

	// get the translation
	se3.translation() = cr - se3.so3() * cl;

	return true;
}

class Sim3Estimator: public theia::Estimator<unsigned int, Sophus::Sim3f> {

	const Eigen::Vector3f* pl_;
	const Eigen::Vector3f* pr_;

public:
	Sim3Estimator(const std::vector<Eigen::Vector3f>& pl, const std::vector<Eigen::Vector3f>& pr) {
		pl_ = pl.data();
		pr_ = pr.data();
	}

	// we need at least 3 correspondences for a proper sim3
	double SampleSize() const {
		return 3;
	}

	// estimate a candidate sim3 from correspondences
	bool EstimateModel(const std::vector<unsigned int>& idxs,
			std::vector<Sophus::Sim3f>* sim3s) const {

		// FIXME this is stupid!
		std::vector<Eigen::Vector3f> pls(idxs.size()), prs(idxs.size());
		for (int ii = 0; ii < idxs.size(); ii++) {
			pls[ii] = pl_[idxs[ii]];
			prs[ii] = pr_[idxs[ii]];
		}

		Sophus::Sim3f sim3;
		if (Sim3fromCorrspondence(pls, prs, sim3)) {
			sim3s->emplace_back(sim3);
			return true;
		}
		return false;
	}

	// The error for a point given a plane model is the point-to-plane distance.
	double Error(const unsigned int& idx, const Sophus::Sim3f& sim3) const {
		return (pr_[idx] - sim3 * pl_[idx]).norm();
	}

private:
	Sim3Estimator(const Sim3Estimator&);
	void operator=(const Sim3Estimator&);

};

float Sim3_ransac(std::vector<Eigen::Vector3f>& pl, std::vector<Eigen::Vector3f>& pr,
		std::vector<unsigned int>& pairs, Sophus::Sim3f& sim3, const float error_thresh, theia::RansacSummary& ransac_summary) {

	theia::RansacParameters ransac_params;
	ransac_params.error_thresh = error_thresh;
	Sim3Estimator sim3estimator(pl, pr);

	std::unique_ptr<theia::SampleConsensusEstimator<Sim3Estimator> > ransac =
			theia::CreateAndInitializeRansacVariant(theia::RansacType::RANSAC, ransac_params, sim3estimator);
	// Estimate the dominant plane.
	if (!ransac->Estimate(pairs, &sim3, &ransac_summary)){
		LOG(WARNING) << "Ransac failed";
		return 0;
	}

	// get the error
	return ransac_summary.confidence;

}

}

#endif /* INCLUDE_HPMVS_ESTIMATESE3_HPP_ */
