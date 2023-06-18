//#######################################################################################
//#                                                                                     #
//#                              CLOUDCOMPARE PLUGIN: qTreeIso                          #
//#                                                                                     #
//#        This program is free software; you can redistribute it and/or modify         #
//#        it under the terms of the GNU General Public License as published by         #
//#        the Free Software Foundation; version 2 or later of the License.             #
//#                                                                                     #
//#        This program is distributed in the hope that it will be useful,              #
//#        but WITHOUT ANY WARRANTY; without even the implied warranty of               #
//#        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                 #
//#        GNU General Public License for more details.                                 #
//#                                                                                     #
//#        Please cite the following paper if you find this tool helpful                #
//#                                                                                     #
//#        Xi, Z.; Hopkinson, C. 3D Graph-Based Individual-Tree Isolation (Treeiso)     #
//#        from Terrestrial Laser Scanning Point Clouds. Remote Sens. 2022, 14, 6116.   #
//#        https://doi.org/10.3390/rs14236116                                           #
//#                                                                                     #
//#		   Our work relies on the cut-pursuit algorithm, please also consider citing:   #
//#        Landrieu, L.; Obozinski, G. Cut Pursuit: Fast Algorithms to Learn Piecewise  #
//#        Constant Functions on General Weighted Graphs. SIAM J. Imaging Sci.          #
//#        2017, 10, 1724–1766.                                                         #
//#                                                                                     #
//#                                     Copyright ©                                     #
//#                  Artemis Lab, Department of Geography & Environment                 #
//#                            University of Lethbridge, Canada                         #
//#                                                                                     #
//#                                                                                     #
//#                           Zhouxin Xi and Chris Hopkinson;                           #
//#                    truebelief2010@gmail.com; c.hopkinson@uleth.ca                   #
//#                                                                                     #
//#######################################################################################

// A Matlab version shared via:
// https://github.com/truebelief/artemis_treeiso

//TreeIso
#include "TreeIso.h"
#include "TreeIsoHelper.h"

//CC
#include <ccMainAppInterface.h>
#include <ccQtHelpers.h>

//qCC_db
#include <ccPointCloud.h>

//Qt
#include <QProgressDialog>
#include <QCoreApplication>
#include <QElapsedTimer>

//system
#include <cmath>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <iostream>

//Boost
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>

typedef boost::geometry::model::d2::point_xy<double> point_xy;
typedef boost::geometry::model::polygon<point_xy> polygon;
typedef boost::geometry::model::multi_point<point_xy> multi_point;
typedef boost::geometry::model::multi_polygon <polygon> multi_polygon;

typedef std::vector<float> Vec3d;

//custom
using namespace CP;

//scalar field names
static const char InitSegsSFName[] = "init_segs";
static const char IntermedSegsSFName[] = "intermediate_segs";
static const char FinalSegsSFName[] = "final_segs";

template <	class result_t = std::chrono::milliseconds,
			class clock_t = std::chrono::steady_clock,
			class duration_t = std::chrono::milliseconds >
static auto Since(std::chrono::time_point<clock_t, duration_t> const& start)
{
	return std::chrono::duration_cast<result_t>(clock_t::now() - start);
}

bool TreeIso::Init_seg_pcd(ccPointCloud* pc, const unsigned min_nn1, const float regStrength1, const float PR_DECIMATE_RES1, QProgressDialog* progressDlg/*=nullptr*/)
{
	if (!pc)
	{
		assert(false);
		return false;
	}

	if (progressDlg)
	{
		progressDlg->setRange(0, 100);
	}

	auto start = std::chrono::steady_clock::now();     // start timer

	const unsigned pointCount = pc->size();

	std::vector<Vec3d> pc_vec;
	toTranslatedVector(pc, pc_vec);

	if (progressDlg)
	{
		progressDlg->setValue(10);
		QCoreApplication::processEvents();
	}

	std::vector<size_t> ia, ic;
	std::vector<Vec3d> pc_dec, pc_sub;
	decimate_vec(pc_vec, PR_DECIMATE_RES1, pc_dec);
	unique_index_by_rows(pc_dec, ia, ic);
	get_subset(pc_vec, ia, pc_sub);

	if (progressDlg)
	{
		progressDlg->setValue(30);
		QCoreApplication::processEvents();
	}

	const unsigned K = (min_nn1 - 1);
	Vec3d edgeWeight;
	std::vector<uint32_t> Eu;
	std::vector<uint32_t> Ev;
	std::vector<uint32_t> in_component;
	std::vector<std::vector<uint32_t>> components;

	perform_cut_pursuit(K, regStrength1, pc_sub, edgeWeight, Eu, Ev, in_component, components);

	if (progressDlg)
	{
		progressDlg->setValue(90);
		QCoreApplication::processEvents();
	}

	//export segments as a new scalar field
	int outSFIndex = pc->getScalarFieldIndexByName(InitSegsSFName);
	if (outSFIndex < 0)
	{
		outSFIndex = pc->addScalarField(InitSegsSFName);
		if (outSFIndex < 0)
		{
			ccLog::Error("[TreeIso] Not enough memory!");
			return false;
		}
	}
	CCCoreLib::ScalarField* outSF = pc->getScalarField(outSFIndex);
	outSF->fill(CCCoreLib::NAN_VALUE);

	std::vector<uint32_t> clusterIdx(pointCount);
	for (unsigned i = 0; i < pointCount; ++i)
	{
		clusterIdx[i] = in_component[ic[i]];
		outSF->setValue(i, in_component[ic[i]]);
	}
	outSF->computeMinAndMax();
	pc->colorsHaveChanged();
	pc->setCurrentDisplayedScalarField(outSFIndex);
	pc->showSF(true);

	auto elapsed = Since(start).count() / 1000;
	ccLog::Print(QString("[TreeIso] Init segs took %1 seconds").arg(elapsed));

	if (progressDlg)
	{
		progressDlg->setValue(100);
		QCoreApplication::processEvents();
	}

	return true;
}

bool TreeIso::Intermediate_seg_pcd(ccPointCloud* pc, const unsigned PR_MIN_NN2, const float PR_REG_STRENGTH2, const float PR_DECIMATE_RES2, const float PR_MAX_GAP, QProgressDialog* progressDlg/*=nullptr*/)
{
	if (!pc)
	{
		assert(false);
		return false;
	}

	if (progressDlg)
	{
		progressDlg->setRange(0, 100);
	}

	auto start = std::chrono::steady_clock::now();

	unsigned pointCount = pc->size();

	int initSFIndex = pc->getScalarFieldIndexByName(InitSegsSFName);
	if (initSFIndex < 0)
	{
		ccLog::Error("[TreeIso] Please run initial segmentation first!");
		return false;
	}
	CCCoreLib::ScalarField* initSF = pc->getScalarField(initSFIndex);

	initSF->size();
	std::vector<uint32_t> in_component;
	try
	{
		in_component.resize(pointCount);
	}
	catch (const std::bad_alloc&)
	{
		ccLog::Error("[TreeIso] Not enough memory");
		return false;
	}
	for (unsigned i = 0; i < pointCount; ++i)
	{
		in_component[i] = initSF->getValue(i);
	}
	std::vector<std::vector<uint32_t>> clusterVGroup;
	unique_group(in_component, clusterVGroup);

	std::vector<Vec3d> pc_vec;
	toTranslatedVector(pc, pc_vec);

	size_t n_clusters = clusterVGroup.size();

	std::vector<Vec3d> clusterCentroids(n_clusters);
	std::vector<std::vector<Vec3d>> currentClusterDecs(n_clusters);
	std::vector<std::vector<size_t>> currentClusterDecsICs(n_clusters);

	for (size_t i = 0; i < n_clusters; ++i)
	{
		std::vector<Vec3d> currentClusterPos;
		get_subset(pc_vec, clusterVGroup[i], currentClusterPos);

		if (currentClusterPos.size() > 1)
		{
			std::vector<size_t> ia, ic;
			std::vector<Vec3d> currentClusterPosDec;
			decimate_vec(currentClusterPos, PR_DECIMATE_RES2, currentClusterPosDec);
			unique_index_by_rows(currentClusterPosDec, ia, ic);

			std::vector<Vec3d> currentClusterPosUnq;
			get_subset(currentClusterPos, ia, currentClusterPosUnq);
			currentClusterDecs[i] = currentClusterPosUnq;

			mean_col(currentClusterPos, clusterCentroids[i]);
			currentClusterDecsICs[i] = ic;
		}
		else
		{
			currentClusterDecs[i] = currentClusterPos;
			clusterCentroids[i] = currentClusterPos[0];
			currentClusterDecsICs[i] = std::vector<size_t>(1, 0);
		}
	}

	if (progressDlg)
	{
		progressDlg->setValue(10);
		QCoreApplication::processEvents();
	}

	std::vector<std::vector<uint32_t>> minIdxsC;
	std::vector<Vec3d> minIdxsD;
	knn_cpp_nearest_neighbors(clusterCentroids, PR_MIN_NN2, minIdxsC, minIdxsD, 8);

	if (progressDlg)
	{
		progressDlg->setValue(15);
		QCoreApplication::processEvents();
	}

	size_t n_centroids = minIdxsC.size();
	if (n_centroids == 0)
	{
		assert(false);
		return false;
	}
	size_t n_K = minIdxsC[0].size();
	std::vector<Eigen::MatrixXf> currentClusterDecMats;
	std::vector<knncpp::KDTreeMinkowskiX<float, knncpp::EuclideanDistance<float>>> knn_kdtrees;
	for (size_t i = 0; i < n_centroids; ++i)
	{
		std::vector<Vec3d> currentClusterDec = currentClusterDecs[i];
		Eigen::MatrixXf currentClusterDecMat(currentClusterDec[0].size(), currentClusterDec.size());
		for (size_t k = 0; k < currentClusterDec.size(); k++)
		{
			currentClusterDecMat.col(k) = Eigen::VectorXf::Map(&currentClusterDec[k][0], currentClusterDec[k].size());
		}
		currentClusterDecMats.push_back(currentClusterDecMat);
	}

	if (progressDlg)
	{
		progressDlg->setValue(20);
		QCoreApplication::processEvents();
	}

	std::vector<Vec3d> nnDists;
	nnDists.resize(n_centroids, Vec3d(n_K));
	for (size_t i = 0; i < n_centroids; ++i)
	{
		knncpp::KDTreeMinkowskiX<float, knncpp::EuclideanDistance<float>> knn_kdtree(currentClusterDecMats[minIdxsC[i][0]]);
		knn_cpp_build(knn_kdtree);
		for (size_t j = 1; j < n_K; ++j)
		{
			if (minIdxsD[i][j] > 0)
			{
				std::vector<std::vector<size_t>> min_c;
				std::vector<Vec3d> min_D2;
				float min_D = knn_cpp_query_min_d(knn_kdtree, currentClusterDecMats[minIdxsC[i][j]], 1);
				nnDists[i][j] = min_D;
			}
		}
	}

	if (progressDlg)
	{
		progressDlg->setValue(40);
		QCoreApplication::processEvents();
	}

	std::vector<Vec3d> currentClusterDecsFlat;
	std::vector<size_t> currentClusterDecsFlatIndex;
	std::vector<size_t> currentClusterDecsICsFlatIndex;
	std::vector<uint32_t> currentClusterDecsIDs;

	for (auto& vec : clusterVGroup)
	{
		currentClusterDecsIDs.insert(currentClusterDecsIDs.end(), std::make_move_iterator(vec.begin()), std::make_move_iterator(vec.end()));
	}
	for (auto& vec : currentClusterDecs)
	{
		currentClusterDecsFlat.insert(currentClusterDecsFlat.end(), std::make_move_iterator(vec.begin()), std::make_move_iterator(vec.end()));
		currentClusterDecsICsFlatIndex.insert(currentClusterDecsICsFlatIndex.end(), vec.size());
	}

	for (unsigned i = 0; i < currentClusterDecs.size(); i++)
	{
		for (unsigned j = 0; j < currentClusterDecs[i].size(); ++j)
		{
			currentClusterDecsFlatIndex.push_back(i);
		}
	}

	if (progressDlg)
	{
		progressDlg->setValue(50);
		QCoreApplication::processEvents();
	}

	size_t nNodes = currentClusterDecsFlat.size();
	unsigned nKs = PR_MIN_NN2;

	std::vector<std::vector<uint32_t>> minIdxs;
	std::vector<Vec3d> Ds;
	knn_cpp_nearest_neighbors(currentClusterDecsFlat, PR_MIN_NN2, minIdxs, Ds, 8);

	Vec3d edgeWeight;
	std::vector<uint32_t> Eu;
	std::vector<uint32_t> Ev;

	for (size_t i = 0; i < minIdxs.size(); i++)
	{
		size_t currentNode = currentClusterDecsFlatIndex[i];
		Vec3d currentDists = nnDists[currentNode];
		for (size_t j = 1; j < minIdxs[0].size(); ++j)
		{
			size_t nnNode = currentClusterDecsFlatIndex[minIdxs[i][j]];
			const std::vector<uint32_t>& nnCand = minIdxsC[currentNode];
			auto it = std::find(nnCand.begin(), nnCand.end(), nnNode);
			if (it != nnCand.end())
			{
				float nnDist = currentDists[it - nnCand.begin()];
				if (nnDist < PR_MAX_GAP)
				{
					Eu.push_back(static_cast<uint32_t>(i));
					Ev.push_back(minIdxs[i][j]);
					edgeWeight.push_back(10 / ((nnDist + 0.1) / 0.01));
				}
			}
		}
	}

	if (progressDlg)
	{
		progressDlg->setValue(55);
		QCoreApplication::processEvents();
	}

	std::vector<uint32_t> in_component2d;
	perform_cut_pursuit2d(PR_MIN_NN2, PR_REG_STRENGTH2, currentClusterDecsFlat, edgeWeight, Eu, Ev, in_component2d);

	if (progressDlg)
	{
		progressDlg->setValue(80);
		QCoreApplication::processEvents();
	}

	std::vector<uint32_t> currentClusterDecsICsReverse;
	std::vector<Vec3d> currentClusterDecsReverse;
	size_t ia_counter = 0;
	for (size_t i = 0; i < currentClusterDecs.size(); i++)
	{
		size_t N = currentClusterDecsICs[i].size();
		std::vector<uint32_t> v(in_component2d.begin() + ia_counter, in_component2d.begin() + ia_counter + currentClusterDecsICsFlatIndex[i]);
		for (size_t j = 0; j < N; j++)
		{
			currentClusterDecsICsReverse.push_back(v[currentClusterDecsICs[i][j]]);
		}
		ia_counter += currentClusterDecsICsFlatIndex[i];
	}

	std::vector<uint32_t> currentClusterDecsIDsSortedIdx;
	std::vector<uint32_t> currentClusterDecsIDsSorted;

	sort_indexes(currentClusterDecsIDs, currentClusterDecsIDsSortedIdx, currentClusterDecsIDsSorted);
	std::vector<uint32_t> currentClusterDecsICsReverseSorted;
	currentClusterDecsICsReverseSorted.reserve(currentClusterDecsICsReverse.size());
	std::transform(currentClusterDecsIDsSortedIdx.begin(), currentClusterDecsIDsSortedIdx.end(), std::back_inserter(currentClusterDecsICsReverseSorted),
		[&](uint32_t i) { return currentClusterDecsICsReverse[i]; });

	//export segments as a new scalar field
	int outSFIndex = pc->getScalarFieldIndexByName(IntermedSegsSFName);
	if (outSFIndex < 0)
	{
		outSFIndex = pc->addScalarField(IntermedSegsSFName);
		if (outSFIndex < 0)
		{
			ccLog::Error("[TreeIso] Not enough memory!");
			return false;
		}
	}
	CCCoreLib::ScalarField* outSF = pc->getScalarField(outSFIndex);
	outSF->fill(CCCoreLib::NAN_VALUE);

	std::vector<uint32_t> groupIdx(pointCount);
	for (unsigned i = 0; i < pointCount; ++i)
	{
		groupIdx[i] = currentClusterDecsICsReverseSorted[i];
		outSF->setValue(i, currentClusterDecsICsReverseSorted[i]);
	}
	outSF->computeMinAndMax();
	pc->colorsHaveChanged();
	pc->setCurrentDisplayedScalarField(outSFIndex);
	pc->showSF(true);

	if (progressDlg)
	{
		progressDlg->setValue(100);
		QCoreApplication::processEvents();
	}

	auto elapsed = Since(start).count() / 1000;
	ccLog::Print(QString("[TreeIso] Intermediate segs took %1 seconds").arg(elapsed));

	return true;
}

bool TreeIso::Final_seg_pcd(ccPointCloud* pc, const unsigned PR_MIN_NN3, const float PR_REL_HEIGHT_LENGTH_RATIO, const float PR_VERTICAL_WEIGHT, QProgressDialog* progressDlg/*=nullptr*/)
{
	if (!pc)
	{
		assert(false);
		return false;
	}

	if (progressDlg)
	{
		progressDlg->setRange(0, 100);
	}

	auto start = std::chrono::steady_clock::now();

	unsigned pointCount = pc->size();
	std::vector<Vec3d> pc_vec;
	toTranslatedVector(pc, pc_vec);

	int initIdx = pc->getScalarFieldIndexByName(InitSegsSFName);
	if (initIdx < 0)
	{
		ccLog::Error("[TreeIso] Please run initial segmentation first!");
		return false;
	}
	int groupIdx = pc->getScalarFieldIndexByName(IntermedSegsSFName);
	if (groupIdx < 0)
	{
		ccLog::Error("[TreeIso] Please run intermeditate segmentation first!");
		return false;
	}

	CCCoreLib::ScalarField* initSF = pc->getScalarField(initIdx);
	initSF->size();
	std::vector<uint32_t> segs_init_ids;
	segs_init_ids.resize(pointCount);
	for (unsigned i = 0; i < pointCount; ++i)
	{
		segs_init_ids[i] = initSF->getValue(i);
	}

	CCCoreLib::ScalarField* groupSF = pc->getScalarField(groupIdx);
	if (!groupSF)
	{
		assert(false);
		return false;
	}

	groupSF->size();
	std::vector<uint32_t> segs_group_ids;
	segs_group_ids.resize(pointCount);

	for (unsigned i = 0; i < pointCount; ++i)
	{
		segs_group_ids[i] = groupSF->getValue(i);
	}
	
	std::vector<std::vector<uint32_t>> initVGroup;
	std::vector<uint32_t> initU;
	std::vector<uint32_t> initUI;
	unique_group(segs_init_ids, initVGroup, initU, initUI);

	size_t n_init_clusters = initVGroup.size();
	std::vector<Vec3d> clusterCentroids(n_init_clusters);

	for (size_t i = 0; i < n_init_clusters; ++i)
	{
		std::vector<Vec3d> clusterPts;
		get_subset(pc_vec, initVGroup[i], clusterPts);
		Vec3d clusterCentroid;
		mean_col(clusterPts, clusterCentroid);
		clusterCentroids[i] = clusterCentroid;
		
		std::vector<uint32_t> segs_group_id;
		get_subset(segs_group_ids, initVGroup[i], segs_group_id);
		float seg_group_mode = mode_col(segs_group_id);
		for (int j = 0; j < segs_group_id.size(); ++j) {
			segs_group_ids[initVGroup[i][j]] = seg_group_mode;
		}
	}

	std::vector<uint32_t> cluster_ids;
	get_subset(segs_group_ids, initUI, cluster_ids);

	std::vector<std::vector<uint32_t>> clusterVGroup;
	std::vector<uint32_t> clusterU0;
	unique_group(cluster_ids, clusterVGroup, clusterU0);

	if (clusterVGroup.size() > 1)
	{
		auto start = std::chrono::steady_clock::now();
		size_t n_clusters = clusterVGroup.size();

		std::vector<unsigned> mergedRemainIds;

		size_t n_to_merge_ids = 1;
		int n_prev_merge_ids = -1;
		size_t iter = 1;

		std::vector<std::vector<uint32_t>> groupVGroup;
		std::vector<uint32_t> groupU;

		std::vector<std::vector<float>> centroid2DFeatures;
		Eigen::MatrixXf groupCentroidsToMerge;
		Eigen::MatrixXf groupCentroidsRemain;

		while ((n_to_merge_ids != 0) && (static_cast<int>(n_to_merge_ids) != n_prev_merge_ids))
		{
			if (iter > 1)
			{
				n_prev_merge_ids = static_cast<int>(n_to_merge_ids);
			}
			unique_group(segs_group_ids, groupVGroup, groupU);
			size_t nGroups = groupVGroup.size();
			centroid2DFeatures.resize(nGroups, std::vector<float>(2));

			std::vector<float> zFeatures;
			zFeatures.resize(nGroups);
			std::vector<float> lenFeatures;
			lenFeatures.resize(nGroups);

			std::vector<polygon> groupHulls;
			for (size_t i = 0; i < nGroups; ++i)
			{
				std::vector<Vec3d> groupPts;
				get_subset(pc_vec, groupVGroup[i], groupPts);
				Vec3d groupCentroids;
				mean_col(groupPts, groupCentroids);
				centroid2DFeatures[i][0] = groupCentroids[0];
				centroid2DFeatures[i][1] = groupCentroids[1];

				std::vector<float> minPts;
				min_col(groupPts, minPts);
				zFeatures[i] = minPts[2];
				std::vector<float> maxPts;
				max_col(groupPts, maxPts);
				lenFeatures[i] = maxPts[2] - minPts[2];

				polygon hull;
				multi_point conv_points;
				for (const auto& p : groupPts)
				{
					conv_points.push_back(point_xy(p[0], p[1]));
				}

				boost::geometry::convex_hull(conv_points, hull);
				groupHulls.push_back(hull);
			}

			size_t knncpp_nn = (PR_MIN_NN3 < n_clusters ? PR_MIN_NN3 : n_clusters);
			std::vector<std::vector<uint32_t>> groupNNIdxC;
			std::vector<Vec3d> groupNNCDs;
			knn_cpp_nearest_neighbors(centroid2DFeatures, knncpp_nn, groupNNIdxC, groupNNCDs, 8);
			Vec3d mds;
			mean_col(groupNNCDs, mds);
			float sigmaD = mds[1];

			std::vector<size_t> toMergeIds;
			std::vector<size_t> toMergeCandidateIds;
			Vec3d toMergeCandidateMetrics;
			for (size_t i = 0; i < nGroups; ++i)
			{
				const std::vector<uint32_t>& nnGroupId = groupNNIdxC[i];
				Vec3d nnGroupZ(nnGroupId.size());
				Vec3d nnGroupLen(nnGroupId.size());
				for (size_t j = 0; j < nnGroupId.size(); ++j)
				{
					nnGroupZ[j] = zFeatures[nnGroupId[j]];
					nnGroupLen[j] = lenFeatures[nnGroupId[j]];
				}
				float minZ = min_col(nnGroupZ);
				float minLen = min_col(nnGroupLen);

				float currentGroupRelHt = (zFeatures[i] - minZ) / lenFeatures[i];

				if (abs(currentGroupRelHt) > PR_REL_HEIGHT_LENGTH_RATIO)
				{
					if (iter == 1)
					{
						float initialLenRatio = lenFeatures[i] / median_col(lenFeatures);
						if (initialLenRatio > 1.5)
						{
							toMergeIds.push_back(i);
						}
						toMergeCandidateIds.push_back(i);
						toMergeCandidateMetrics.push_back(initialLenRatio);
					}
					else
					{
						toMergeIds.push_back(i);
					}

				}
			}
			if ((iter == 1) && (toMergeCandidateMetrics.size() == 0))
			{
				break;
			}
			if ((iter == 1) && (toMergeIds.size() == 0))
			{
				size_t cand_ind = arg_max_col(toMergeCandidateMetrics);
				toMergeIds.push_back(toMergeCandidateIds[cand_ind]);
			}

			std::vector<size_t> remainIds;
			std::vector<size_t> allIds(nGroups);
			std::iota(std::begin(allIds), std::end(allIds), 0);
			std::set_difference(allIds.begin(), allIds.end(), toMergeIds.begin(), toMergeIds.end(), std::inserter(remainIds, remainIds.begin()));

			get_subset(centroid2DFeatures, toMergeIds, groupCentroidsToMerge);
			get_subset(centroid2DFeatures, remainIds, groupCentroidsRemain);
			size_t knncpp_nn2 = (PR_MIN_NN3 < remainIds.size() ? PR_MIN_NN3 : remainIds.size());


			knncpp::KDTreeMinkowskiX<float, knncpp::EuclideanDistance<float>> knn_kdtree(groupCentroidsRemain);
			knn_cpp_build(knn_kdtree);
			std::vector<std::vector<size_t>> groupNNIdx;
			std::vector<Vec3d> groupNNIdxDists;

			knn_cpp_query(knn_kdtree, groupCentroidsToMerge, knncpp_nn2, groupNNIdx, groupNNIdxDists);

			size_t nToMergeIds = toMergeIds.size();
			size_t nRemainIds = remainIds.size();
			for (size_t i = 0; i < nToMergeIds; ++i)
			{
				size_t toMergeId = toMergeIds[i];

				Eigen::MatrixXf currentClusterCentroids;
				get_subset(clusterCentroids, clusterVGroup[toMergeId], currentClusterCentroids);

				size_t nNNs = groupNNIdx.size();

				std::vector<float> scores;
				std::vector<size_t> filteredRemainIds;
				std::vector<float> min3DSpacings;

				for (size_t j = 0; j < nNNs; ++j)
				{
					size_t remainId = remainIds[groupNNIdx[j][i]];

					float lineSegs2 = zFeatures[toMergeId] + lenFeatures[toMergeId] - zFeatures[remainId];
					float lineSegs1 = zFeatures[remainId] + lenFeatures[remainId] - zFeatures[toMergeId];

					float verticalOverlapRatio = (lineSegs2 > lineSegs1 ? lineSegs1 : lineSegs2) / (lineSegs1 > lineSegs2 ? lineSegs1 : lineSegs2);
					float horizontalOverlapRatio;
					if ((boost::geometry::num_points(groupHulls[toMergeId]) > 3) & (boost::geometry::num_points(groupHulls[remainId]) > 3))
					{
						multi_polygon intersection;
						boost::geometry::intersection(groupHulls[toMergeId], groupHulls[remainId], intersection);
						float intersect_area = boost::geometry::area(intersection);
						float area1 = boost::geometry::area(groupHulls[toMergeId]);
						float area2 = boost::geometry::area(groupHulls[remainId]);
						horizontalOverlapRatio = intersect_area / (area1 < area2 ? area1 : area2);
					}
					else
					{
						horizontalOverlapRatio = 0.0;
					}

					Eigen::MatrixXf nnClusterCentroids;
					get_subset(clusterCentroids, clusterVGroup[remainId], nnClusterCentroids);

					knncpp::KDTreeMinkowskiX<float, knncpp::EuclideanDistance<float>> knn_kdtree2(nnClusterCentroids);
					knn_cpp_build(knn_kdtree2);
					std::vector<std::vector<size_t>> min3D_idx;
					std::vector<Vec3d> min3D_dists;

					knn_cpp_query(knn_kdtree2, currentClusterCentroids, 1, min3D_idx, min3D_dists);
					float min3DSpacing = min_col(min3D_dists[0]);
					min3DSpacings.push_back(min3DSpacing);

					Eigen::MatrixXf nnClusterCentroids2D = nnClusterCentroids.block(0, 0, 2, nnClusterCentroids.cols());
					Eigen::MatrixXf currentClusterCentroids2D = currentClusterCentroids.block(0, 0, 2, currentClusterCentroids.cols());
					Eigen::VectorXf nnClusterCentroids2Dmean = nnClusterCentroids2D.rowwise().mean();
					Eigen::VectorXf currentClusterCentroids2Dmean = currentClusterCentroids2D.rowwise().mean();
					float min2DSpacing = sqrt((nnClusterCentroids2Dmean[0] - currentClusterCentroids2Dmean[0]) * (nnClusterCentroids2Dmean[0] - currentClusterCentroids2Dmean[0]) + (nnClusterCentroids2Dmean[1] - currentClusterCentroids2Dmean[1]) * (nnClusterCentroids2Dmean[1] - currentClusterCentroids2Dmean[1]));

					verticalOverlapRatio = verticalOverlapRatio > 0 ? verticalOverlapRatio : 0;

					float score = exp(-(1 - horizontalOverlapRatio) * (1 - horizontalOverlapRatio) - PR_VERTICAL_WEIGHT * (1 - verticalOverlapRatio) * (1 - verticalOverlapRatio) - ((min3DSpacing < min2DSpacing ? min3DSpacing : min2DSpacing) / sigmaD) * ((min3DSpacing < min2DSpacing ? min3DSpacing : min2DSpacing) / sigmaD));
					scores.push_back(score);

					filteredRemainIds.push_back(remainId);
				}


				std::vector<size_t> scoreSortI;
				std::vector<float> scoreSort;

				sort_indexes<float, size_t>(scores, scoreSortI, scoreSort);
				float score_highest = scoreSort[0];
				if (score_highest == 0)
				{
					continue;
				}
				std::vector<float> scoreSortRatio;

				scoreSortRatio.reserve(scoreSort.size());
				std::transform(scoreSort.begin(), scoreSort.end(), std::back_inserter(scoreSortRatio),
					[score_highest](float value) { return value / score_highest; });

				std::vector<std::size_t> scoreSortCandidateIdx;
				for (std::size_t i = 0; i < scoreSortRatio.size(); ++i)
				{
					if (scoreSortRatio[i] > 0.7)
					{
						scoreSortCandidateIdx.push_back(i);
					}
				}
				size_t nScoreSortCandidateIdx = scoreSortCandidateIdx.size();
				unsigned mergeNNId = 0;
				if (nScoreSortCandidateIdx == 1)
				{
					mergeNNId = groupU[filteredRemainIds[scoreSortI[scoreSortCandidateIdx[0]]]];

				}
				else if (nScoreSortCandidateIdx > 1)
				{
					std::vector<float> min3DSpacingsFiltered;
					std::vector<size_t> scoreSortIFiltered;
					get_subset(scoreSortI, scoreSortCandidateIdx, scoreSortIFiltered);
					get_subset(min3DSpacings, scoreSortIFiltered, min3DSpacingsFiltered);

					size_t filterMinSpacingIdx = arg_min_col(min3DSpacingsFiltered);
					mergeNNId = groupU[filteredRemainIds[scoreSortI[filterMinSpacingIdx]]];
				}
				else
				{
					continue;
				}

				std::vector<uint32_t> currentVGroup = groupVGroup[toMergeIds[i]];
				for (unsigned k = 0; k < currentVGroup.size(); ++k)
				{
					segs_group_ids[currentVGroup[k]] = mergeNNId;
				}
				mergedRemainIds.push_back(mergeNNId);
				mergedRemainIds.push_back(groupU[toMergeIds[i]]);
			}

			n_to_merge_ids = toMergeIds.size();
			get_subset(segs_group_ids, initUI, cluster_ids);
			unique_group(cluster_ids, clusterVGroup, clusterU0);
			n_clusters = clusterVGroup.size();
			++iter;
		}

		unique_group(segs_group_ids, groupVGroup, groupU);
		for (unsigned j = 0; j < groupVGroup.size(); ++j)
		{
			std::vector<uint32_t> currentVGroup = groupVGroup[j];
			for (unsigned k = 0; k < currentVGroup.size(); ++k)
			{
				segs_group_ids[currentVGroup[k]] = j + 1;
			}
		}

		//export segments as a new scalar field
		int outSFIndex = pc->getScalarFieldIndexByName(FinalSegsSFName);
		if (outSFIndex < 0)
		{
			outSFIndex = pc->addScalarField(FinalSegsSFName);
			if (outSFIndex < 0)
			{
				ccLog::Error("[TreeIso] Not enough memory!");
				return false;
			}
		}
		CCCoreLib::ScalarField* outSF = pc->getScalarField(outSFIndex);
		outSF->fill(CCCoreLib::NAN_VALUE);

		std::vector<uint32_t> groupIdx(pointCount);
		for (unsigned i = 0; i < pointCount; ++i)
		{
			outSF->setValue(i, segs_group_ids[i]);
		}
		outSF->computeMinAndMax();
		pc->colorsHaveChanged();
		pc->setCurrentDisplayedScalarField(outSFIndex);
		pc->showSF(true);

		if (progressDlg)
		{
			progressDlg->setValue(100);
			QCoreApplication::processEvents();
		}

		auto elapsed = Since(start).count() / 1000;
		ccLog::Print(QString("[TreeIso] Final segs took: %1 seconds !!!").arg(elapsed));
	}

	return true;
}

//1. initial 3D segmentation
bool TreeIso::Init_seg(	const unsigned min_nn1,
						const float regStrength1,
						const float PR_DECIMATE_RES1,
						ccMainAppInterface* app,
						QProgressDialog* progressDlg)
{
	if (!app)
	{
		assert(false);
		return false;
	}

	const ccHObject::Container& selectedEntities = app->getSelectedEntities();
	if (selectedEntities.empty())
	{
		assert(false);
		app->dispToConsole("[TreeIso] Select at least one cloud", ccMainAppInterface::ERR_CONSOLE_MESSAGE);
		return false;
	}
	ccHObject* ent = selectedEntities[0];

	if (!ent || !ent->isA(CC_TYPES::POINT_CLOUD))
	{
		app->dispToConsole("[TreeIso] Not a point cloud", ccMainAppInterface::ERR_CONSOLE_MESSAGE);
		return false;
	}
	ccPointCloud* pointCloud = static_cast<ccPointCloud*>(ent);

	if (Init_seg_pcd(pointCloud, min_nn1, regStrength1, PR_DECIMATE_RES1, progressDlg))
	{
		ent->redrawDisplay();
		return true;
	}
	else
	{
		return false;
	}
}

//2. Bottom-up 2D segmentation
bool TreeIso::Intermediate_seg(	const unsigned PR_MIN_NN2,
								const float PR_REG_STRENGTH2,
								const float PR_DECIMATE_RES2,
								const float PR_MAX_GAP,
								ccMainAppInterface* app,
								QProgressDialog* progressDlg)
{
	if (!app)
	{
		assert(false);
		return false;
	}

	const ccHObject::Container& selectedEntities = app->getSelectedEntities();

	if (selectedEntities.empty())
	{
		assert(false);
		app->dispToConsole("[TreeIso] Select at least one cloud", ccMainAppInterface::ERR_CONSOLE_MESSAGE);
		return false;
	}

	ccHObject* ent = selectedEntities[0];

	if (!ent || !ent->isA(CC_TYPES::POINT_CLOUD))
	{
		app->dispToConsole("[TreeIso] Not a point cloud", ccMainAppInterface::ERR_CONSOLE_MESSAGE);
		return false;
	}

	ccPointCloud* pointCloud = static_cast<ccPointCloud*>(ent);

	if (Intermediate_seg_pcd(pointCloud, PR_MIN_NN2, PR_REG_STRENGTH2, PR_DECIMATE_RES2, PR_MAX_GAP, progressDlg))
	{
		ent->redrawDisplay();
		return true;
	}
	else
	{
		return false;
	}
}

//3. Global edge refinement
bool TreeIso::Final_seg(const unsigned PR_MIN_NN3,
						const float PR_REL_HEIGHT_LENGTH_RATIO,
						const float PR_VERTICAL_WEIGHT,
						ccMainAppInterface* app,
						QProgressDialog* progressDlg)
{
	if (!app)
	{
		assert(false);
		return false;
	}

	const ccHObject::Container& selectedEntities = app->getSelectedEntities();

	if (selectedEntities.empty())
	{
		assert(false);
		app->dispToConsole("[TreeIso] Select at least one cloud", ccMainAppInterface::ERR_CONSOLE_MESSAGE);
		return false;
	}

	ccHObject* ent = selectedEntities[0];
	if (!ent || !ent->isA(CC_TYPES::POINT_CLOUD))
	{
		app->dispToConsole("[TreeIso] Not a point cloud", ccMainAppInterface::ERR_CONSOLE_MESSAGE);
		return false;
	}
	ccPointCloud* pointCloud = static_cast<ccPointCloud*>(ent);

	if (Final_seg_pcd(pointCloud, PR_MIN_NN3, PR_REL_HEIGHT_LENGTH_RATIO, PR_VERTICAL_WEIGHT, progressDlg))
	{
		ent->redrawDisplay();
		return true;
	}
	else
	{
		return false;
	}
}

void knn_cpp_build(knncpp::KDTreeMinkowskiX<float, knncpp::EuclideanDistance<float>>& kdtree, unsigned n_thread/*=0*/)
{
	kdtree.setBucketSize(16);
	kdtree.setSorted(true);
	kdtree.setThreads(n_thread);
	kdtree.build();
}

void knn_cpp_query(knncpp::KDTreeMinkowskiX<float, knncpp::EuclideanDistance<float>>& kdtree, Eigen::MatrixXf& query_points, size_t k, std::vector <std::vector<size_t>>& res_idx, std::vector<Vec3d>& res_dists)
{
	res_idx.clear();
	res_dists.clear();

	knncpp::Matrixi indices;
	Eigen::MatrixXf distances;

	kdtree.query(query_points, k, indices, distances);
	res_idx.resize(indices.rows(), std::vector<size_t>(indices.cols()));
	for (Eigen::Index i = 0; i < indices.rows(); i++)
	{
		res_idx[i] = std::vector<size_t>(indices.row(i).data(), indices.row(i).data() + indices.cols());
	}
	res_dists.resize(distances.rows(), Vec3d(distances.cols()));
	for (Eigen::Index i = 0; i < distances.rows(); i++)
	{
		res_dists[i] = Vec3d(distances.row(i).data(), distances.row(i).data() + distances.cols());
	}
}

float knn_cpp_query_min_d(knncpp::KDTreeMinkowskiX<float, knncpp::EuclideanDistance<float>>& kdtree, Eigen::MatrixXf& query_points, size_t k)
{
	knncpp::Matrixi indices;
	Eigen::MatrixXf distances;

	kdtree.query(query_points, k, indices, distances);
	return distances.minCoeff();
}

void knn_cpp_nearest_neighbors(	const std::vector<Vec3d>& dataset,
								size_t k,
								std::vector<std::vector<uint32_t>>& res_idx,
								std::vector<Vec3d>& res_dists,
								unsigned n_thread)
{
	if (dataset.empty())
	{
		assert(false);
		return;
	}

	Eigen::MatrixXf mat(dataset[0].size(), dataset.size());
	for (size_t i = 0; i < dataset.size(); i++)
	{
		mat.col(i) = Eigen::VectorXf::Map(&dataset[i][0], dataset[i].size());
	}

	knncpp::KDTreeMinkowskiX<float, knncpp::EuclideanDistance<float>> kdtree(mat);
	kdtree.setBucketSize(16);
	kdtree.setSorted(true);
	kdtree.setThreads(n_thread);
	kdtree.build();
	knncpp::Matrixi indices;
	Eigen::MatrixXf distances;

	kdtree.query(mat, k, indices, distances);

	res_idx.resize(indices.cols(), std::vector<uint32_t>(indices.rows()));
	for (Eigen::Index i = 0; i < indices.cols(); i++)
	{
		res_idx[i] = std::vector<uint32_t>(indices.col(i).data(), indices.col(i).data() + indices.rows());
	}
	res_dists.resize(distances.cols(), Vec3d(distances.rows()));
	for (Eigen::Index i = 0; i < distances.cols(); i++)
	{
		res_dists[i] = Vec3d(distances.col(i).data(), distances.col(i).data() + distances.rows());
	}
}

template <typename T>
void unique_group(std::vector<T>& arr, std::vector<std::vector<T>>& u_group, std::vector<T>& arr_unq, std::vector<T>& ui)
{
	arr_unq.clear();
	ui.clear();
	u_group.clear();

	if (arr.empty())
	{
		assert(false);
		return;
	}

	std::vector<T> arr_sorted_idx;
	std::vector<T> arr_sorted;
	sort_indexes(arr, arr_sorted_idx, arr_sorted);

	const size_t n = arr.size();

	ui.push_back(arr_sorted_idx[0]);
	std::size_t counter = 0;

	std::vector<T> ut;
	ut.push_back(arr_sorted_idx[0]);

	//detect the location (before sorted) where a row in the sorted array is different from the previous row (as ia), and add one for the reverse index as ic
	for (std::size_t i = 1; i < n; ++i)
	{
		if (arr_sorted[i] != arr_sorted[i - 1])
		{
			ui.push_back(arr_sorted_idx[i]);
			arr_unq.push_back(arr_sorted[i]);

			u_group.push_back(ut);
			ut.clear();
			ut.push_back(arr_sorted_idx[i]);
			counter++;
		}
		else
		{
			ut.push_back(arr_sorted_idx[i]);
		}
	}
	u_group.push_back(ut);
}

template <typename T>
void unique_group(std::vector<T>& arr, std::vector<std::vector<T>>& u_group, std::vector<T>& arr_unq)
{
	arr_unq.clear();
	u_group.clear();

	if (arr.empty())
	{
		assert(false);
		return;
	}

	std::vector<T> arr_sorted_idx;
	std::vector<T> arr_sorted;
	sort_indexes(arr, arr_sorted_idx, arr_sorted);

	const size_t n = arr.size();
	arr_unq.push_back(arr_sorted[0]);
	std::size_t counter = 0;

	std::vector<T>* ut = new std::vector<T>();
	ut->push_back(arr_sorted_idx[0]);

	//detect the location (before sorted) where a row in the sorted array is different from the previous row (as ia), and add one for the reverse index as ic
	for (std::size_t i = 1; i < n; ++i)
	{

		if (arr_sorted[i] != arr_sorted[i - 1])
		{
			arr_unq.push_back(arr_sorted[i]);

			u_group.push_back(*ut);
			ut = new std::vector<T>();
			ut->push_back(arr_sorted_idx[i]);
			counter++;
		}
		else
		{
			ut->push_back(arr_sorted_idx[i]);
		}
	}
	u_group.push_back(*ut);
}

template <typename T>
void unique_group(std::vector<T>& arr, std::vector<std::vector<T>>& u_group)
{
	u_group.clear();

	if (arr.empty())
	{
		assert(false);
		return;
	}

	std::vector<T> arr_sorted_idx;
	std::vector<T> arr_sorted;
	sort_indexes(arr, arr_sorted_idx, arr_sorted);

	const size_t n = arr.size();
	std::size_t counter = 0;

	std::vector<T> ut;
	ut.push_back(arr_sorted_idx[0]);

	//detect the location (before sorted) where a row in the sorted array is different from the previous row (as ia), and add one for the reverse index as ic
	for (std::size_t i = 1; i < n; ++i)
	{
		if (arr_sorted[i] != arr_sorted[i - 1])
		{
			u_group.push_back(ut);
			ut.clear();
			ut.push_back(arr_sorted_idx[i]);
			counter++;
		}
		else
		{
			ut.push_back(arr_sorted_idx[i]);
		}
	}
	u_group.push_back(ut);
}

template <typename T1, typename T2>
void get_subset(std::vector<T1>& arr, std::vector<T2>& indices, std::vector<T1>& arr_sub)
{
	arr_sub.clear();
	for (const auto& idx : indices)
	{
		arr_sub.push_back(arr[idx]);
	}
}

template <typename T1, typename T2>
void get_subset(const std::vector<std::vector<T1>>& arr, const std::vector<T2>& indices, Eigen::MatrixXf& arr_sub)
{
	arr_sub.setZero();

	if (arr.empty())
	{
		assert(false);
		return;
	}

	arr_sub.resize(arr[0].size(), indices.size());

	for (size_t i = 0; i < indices.size(); ++i)
	{
		for (size_t j = 0; j < arr[0].size(); ++j)
		{
			arr_sub(j,i) = arr[indices[i]][j];
		}
	}
}

template <typename T1, typename T2>
void get_subset(const std::vector<std::vector<T1>>& arr, const std::vector<T2>& indices, std::vector<std::vector<T1>>& arr_sub)
{
	arr_sub.clear();

	if (arr.empty() || indices.empty())
	{
		return;
	}

	arr_sub.resize(indices.size());
	for (size_t i = 0; i < indices.size(); ++i)
	{
		arr_sub[i] = arr[indices[i]];
	}
}

template <typename T1, typename T2>
bool get_subset(ccPointCloud* pcd, std::vector<T2>& indices, std::vector<std::vector<T1>>& arr_sub)
{
	arr_sub.clear();
	arr_sub.resize(indices.size(), std::vector<T1>(3));
	for (size_t i = 0; i < indices.size(); ++i)
	{
		const CCVector3* vec = pcd->getPoint(indices[i]);
		arr_sub[i][0] = vec->x;
		arr_sub[i][1] = vec->y;
		arr_sub[i][2] = vec->z;
	}
	return true;
}

//convert ccPointCloud to a vector of points (centering the points about their gravity center)
template <typename T>
void toTranslatedVector(const ccPointCloud* pc, std::vector<std::vector<T>>& y)
{
	y.clear();

	const unsigned pointCount = pc->size();
	if (pointCount == 0)
	{
		assert(false);
		return;
	}
	y.resize(pointCount, std::vector<T>(3));

	for (unsigned i = 0; i < pointCount; ++i)
	{
		const CCVector3* pv = pc->getPoint(i);
		y[i] = {	static_cast<float>(pv->x),
					static_cast<float>(pv->y),
					static_cast<float>(pv->z) };
	}

	std::vector<T> y_mean;
	mean_col(y, y_mean);

	for (unsigned i = 0; i < pointCount; ++i)
	{
		y[i][0] -= y_mean[0];
		y[i][1] -= y_mean[1];
		y[i][2] -= y_mean[2];
	}
}

bool perform_cut_pursuit(	const uint32_t K,
							const float regStrength,
							const std::vector<Vec3d>& pc,
							std::vector<float>& edgeWeight,
							std::vector<uint32_t>& Eu,
							std::vector<uint32_t>& Ev,
							std::vector<uint32_t>& in_component,
							std::vector<std::vector<uint32_t>>& components )
{
	const uint32_t pointCount = static_cast<uint32_t>(pc.size());
	if (pointCount == 0)
	{
		return false;
	}

	std::vector<std::vector<uint32_t>> nn_idx;
	std::vector<Vec3d> nn_D;
	knn_cpp_nearest_neighbors(pc, K + 1, nn_idx, nn_D, 8);

	const uint32_t nNod = pointCount;
	const uint32_t nObs = 3;
	const uint32_t nEdg = pointCount * K;
	const uint32_t cutoff = 0;
	const float mode = 1.0f;
	const float speed = 0;
	const float weight_decay = 0;
	const float verbose = 0;

	if (edgeWeight.size() == 0)
	{
		edgeWeight.resize(nEdg);
		std::fill(edgeWeight.begin(), edgeWeight.end(), 1.0);
	}

	//minus average
	std::vector<double> y_avg(nObs, 0.0);

	if (Eu.size() == 0)
	{
		Eu.resize(nEdg);
		Ev.resize(nEdg);
		for (unsigned i = 0; i < pointCount; ++i)
		{
			y_avg[0] += pc[i][0];
			y_avg[1] += pc[i][1];
			y_avg[2] += pc[i][2];

			for (unsigned j = 0; j < K; ++j)
			{
				Eu[i * K + j] = i;
				Ev[i * K + j] = nn_idx[i][j + 1];
			}
		}
	}
	else
	{
		for (unsigned i = 0; i < pointCount; ++i)
		{
			y_avg[0] += pc[i][0];
			y_avg[1] += pc[i][1];
			y_avg[2] += pc[i][2];
		}
	}

	y_avg[0] /= pointCount;
	y_avg[1] /= pointCount;
	y_avg[2] /= pointCount;

	std::vector<std::vector<float>> y = pc;
	for (unsigned i = 0; i < pointCount; ++i)
	{
		y[i][0] -= static_cast<float>(y_avg[0]);
		y[i][1] -= static_cast<float>(y_avg[1]);
		y[i][2] -= static_cast<float>(y_avg[2]);
	}

	std::vector<float> nodeWeight(pointCount, 1.0);
	std::vector<Vec3d> solution(pointCount, std::vector<float>(K));
	CP::cut_pursuit<float>(nNod, nEdg, nObs, y, Eu, Ev, edgeWeight, nodeWeight, solution, in_component, components, regStrength, cutoff, mode, speed, weight_decay, verbose);
	return true;
}

void perform_cut_pursuit2d(	const uint32_t K,
							const float regStrength,
							const std::vector<Vec3d>& pc_vec,
							Vec3d& edgeWeight,
							std::vector<uint32_t>& Eu,
							std::vector<uint32_t>& Ev,
							std::vector<uint32_t>& in_component )
{
	//build graph for cut-pursuit
	const uint32_t pointCount = static_cast<uint32_t>(pc_vec.size());

	const uint32_t nNod = pointCount;
	const uint32_t nObs = 2;
	const uint32_t nEdgMax = pointCount * K;
	const uint32_t cutoff = 0;
	const float mode = 1;
	const float speed = 0;
	const float weight_decay = 0;
	const float verbose = 0;

	if (edgeWeight.size() == 0)
	{
		edgeWeight.resize(nEdgMax);
		std::fill(edgeWeight.begin(), edgeWeight.end(), 1.0);
	}
	const uint32_t nEdg = static_cast<uint32_t>(edgeWeight.size());

	//minus average
	Vec3d y_avg(nObs, 0.0);

	if (Eu.size() == 0)
	{
		std::vector<std::vector<uint32_t>> nn_idx;
		std::vector<Vec3d> nn_D;
		knn_cpp_nearest_neighbors(pc_vec, K + 1, nn_idx, nn_D, 8);
		Eu.resize(nEdg);
		Ev.resize(nEdg);
		for (uint32_t i = 0; i < pointCount; ++i)
		{
			y_avg[0] += pc_vec[i][0];
			y_avg[1] += pc_vec[i][1];
			for (unsigned j = 0; j < K; ++j)
			{
				Eu[i * K + j] = i;
				Ev[i * K + j] = nn_idx[i][j + 1];
			}
		}
	}
	else
	{
		for (unsigned i = 0; i < pointCount; ++i)
		{
			y_avg[0] += pc_vec[i][0];
			y_avg[1] += pc_vec[i][1];
		}
	}

	y_avg[0] /= pointCount;
	y_avg[1] /= pointCount;

	std::vector<Vec3d> y = pc_vec;
	for (unsigned i = 0; i < pointCount; ++i)
	{
		y[i][0] -= static_cast<float>(y_avg[0]);
		y[i][1] -= static_cast<float>(y_avg[1]);
	}

	std::vector<float> nodeWeight(pointCount, 1.0);
	std::vector<std::vector<uint32_t> > components;
	std::vector<std::vector<float>> solution(pointCount, std::vector<float>(K));
	CP::cut_pursuit<float>(nNod, nEdg, nObs, y, Eu, Ev, edgeWeight, nodeWeight, solution, in_component, components, regStrength, cutoff, mode, speed, weight_decay, verbose);
}

template <typename T>
size_t arg_min_col(std::vector<T>& arr)
{
	auto min_element_it = std::min_element(arr.begin(), arr.end());
	std::size_t min_index = std::distance(arr.begin(), min_element_it);
	return min_index;
}

template <typename T>
size_t arg_max_col(std::vector<T>& arr)
{
	auto max_element_it = std::max_element(arr.begin(), arr.end());
	std::size_t max_index = std::distance(arr.begin(), max_element_it);
	return max_index;
}

template <typename T>
void min_col(std::vector<std::vector<T>>& arr, std::vector<T>& min_vals)
{
	min_vals.clear();
	if (arr.empty())
	{
		return;
	}

	min_vals = arr[0];
	for (size_t j = 1; j < arr.size(); j++)
	{
		const auto& row = arr[j];
		for (size_t i = 0; i < row.size(); i++)
		{
			min_vals[i] = std::min(min_vals[i], row[i]);
		}
	}
}

template <typename T>
T min_col(std::vector<T>& arr)
{
	auto result_it = std::min_element(arr.begin(), arr.end());
	T result=*result_it;
	return result;
}

template <typename T>
T mean_col(std::vector<T>& arr)
{
	if (arr.empty())
	{
		return std::numeric_limits<T>::quite_NaN();
	}

	double sum = 0.0;
	for (const T& value : arr)
	{
		sum += value;
	}

	return static_cast<T>(sum /= arr.size());
}

template <typename T>
T median_col(std::vector<T>& arr)
{
	size_t n = arr.size();
	if (n == 0)
	{
		return std::numeric_limits<T>::quiet_NaN();
	}
	// Sort the vector
	std::sort(arr.begin(), arr.end());
	// Calculate the median
	if (n % 2 == 0)
	{
		return (arr[n / 2 - 1] + arr[n / 2]) / 2;
	}
	else
	{
		return arr[n / 2];
	}	
}

template <typename T>
T mode_col(std::vector<T>& arr) {
	std::unordered_map<T, int> freq;
	for (const auto& val : arr) freq[val]++;
	return std::max_element(freq.begin(), freq.end(), [](const auto& a, const auto& b) { return a.second < b.second; })->first;
}

template <typename T>
void max_col(std::vector<std::vector<T>>& arr, std::vector<T>& max_vals)
{
	max_vals.clear();
	if (arr.empty())
	{
		return;
	}

	max_vals = arr[0];
	for (size_t j = 1; j < arr.size(); j++)
	{
		const auto& row = arr[j];
		for (size_t i = 0; i < row.size(); i++)
		{
			max_vals[i] = std::max(max_vals[i], row[i]);
		}
	}
}

template <typename T>
void mean_col(std::vector<std::vector<T>>& arr, std::vector<T>& mean_vals)
{
	mean_vals.clear();
	if (arr.empty())
	{
		return;
	}

	std::vector<double> sums(arr[0].size(), 0.0);

	for (const std::vector<T>& element : arr)
	{
		for (size_t i = 0; i < arr[0].size(); i++)
		{
			sums[i] += element[i];
		}
	}

	mean_vals = std::vector<T>(arr[0].size());
	for (size_t i = 0; i < arr[0].size(); i++)
	{
		mean_vals[i] = static_cast<T>(sums[i] / arr.size());
	}
}

template <typename T>
void decimate_vec(std::vector<std::vector<T>>& arr, T res, std::vector<std::vector<T>>& vec_dec)
{
	vec_dec.clear();

	if (res == 0)
	{
		assert(false);
		return;
	}

	size_t num_rows = arr.size();
	if (num_rows == 0)
	{
		assert(false);
		return;
	}
	size_t num_cols = arr[0].size();

	std::vector<T> arr_min;
	min_col(arr, arr_min);
	vec_dec.resize(num_rows, std::vector<T>(num_cols));

	for (unsigned i = 0; i < num_rows; ++i)
	{
		for (unsigned j = 0; j < num_cols; ++j)
		{
			vec_dec[i][j] = std::floor((arr[i][j] - arr_min[j]) / res) + 1;
		}
	}
}

//return index
template <typename T>
void sort_indexes_by_row(std::vector<std::vector<T>>& v, std::vector<size_t>& idx, std::vector<std::vector<T>>& v_sorted)
{
	idx.clear();
	v_sorted.clear();

	// initialize original index locations
	size_t m = v.size();
	if (m == 0)
	{
		return;
	}
	size_t n = v[0].size();

	idx.resize(m);
	std::iota(idx.begin(), idx.end(), 0);

	//std::vector<std::vector<T>> v_sorted;
	v_sorted.resize(m, std::vector<T>(n));

	// sort indexes based on comparing values in v
	// using std::stable_sort instead of std::sort
	// to avoid unnecessary index re-orderings
	// when v contains elements of equal values 
	std::stable_sort(idx.begin(), idx.end(),
		[&v](size_t i1, size_t i2)
		{
			//bool vt = true;
			for (size_t k = 0; k < v[0].size(); ++k)
			{
				if (v[i1][k] == v[i2][k])
				{
					continue;
				}
				else
				{
					return v[i1][k] < v[i2][k];
				}
			}
			return false;
		});

	for (size_t i = 0; i < idx.size(); ++i)
	{
		v_sorted[i].resize(n);
		for (size_t k = 0; k < n; ++k)
		{
			v_sorted[i][k] = v[idx[i]][k];
		}
	}
}

//return index
template <typename ValueType, typename IndexType>
void sort_indexes(std::vector<ValueType>& v, std::vector<IndexType>& idx, std::vector<ValueType>& v_sorted)
{
	idx.clear();
	v_sorted.clear();

	size_t m = v.size();
	if (m == 0)
	{
		return;
	}
	idx.resize(m);
	std::iota(idx.begin(), idx.end(), 0);

	v_sorted.resize(m);

	std::stable_sort(idx.begin(), idx.end(),
		[&v](IndexType i1, IndexType i2)
		{
			return v[i1] < v[i2];
		});

	for (size_t i = 0; i < idx.size(); ++i)
	{
		v_sorted[i] = v[idx[i]];
	}
}

template <typename T>
void unique_index_by_rows(std::vector<std::vector<T>>& arr, std::vector<size_t>& ia, std::vector<size_t>& ic)
{
	ia.clear();
	ic.clear();

	std::vector<std::vector<T>> arr_sorted;
	arr_sorted.resize(arr.size());
	std::vector<size_t> sort_idx;

	//sort array first and get indices
	sort_indexes_by_row(arr, sort_idx, arr_sorted);

	const size_t num_rows = arr_sorted.size();
	const size_t num_cols = arr_sorted[0].size();

	ic.resize(num_rows);
	ia.push_back(sort_idx[0]);
	std::size_t counter = 0;
	ic[sort_idx[0]] = counter;

	//detect the location (before sorted) where a row in the sorted array is different from the previous row (as ia), and add one for the reverse index as ic
	for (std::size_t i = 1; i < num_rows; ++i)
	{
		bool diff = false;
		for (std::size_t k = 0; k < num_cols; ++k)
		{
			if (arr_sorted[i][k] != arr_sorted[i - 1][k])
			{
				diff = true;
				break;
			}
		}

		if (diff)
		{
			ia.push_back(sort_idx[i]);
			counter++;
		}
		ic[sort_idx[i]] = counter;
	}
}
