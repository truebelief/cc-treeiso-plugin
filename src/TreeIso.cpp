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

// Matlab and python versions shared via:
// https://github.com/truebelief/artemis_treeiso

//TreeIso
#include "TreeIso.h"
#include "TreeIsoHelper.hpp"

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
#include <unordered_map>
#include <algorithm>

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


struct BBox {
	float minX, minY, maxX, maxY;

	float area() const {
		return (maxX - minX) * (maxY - minY);
	}

	static BBox from_points(const std::vector<Vec3d>& points) {
		BBox bbox = {
			std::numeric_limits<float>::max(),
			std::numeric_limits<float>::max(),
			-std::numeric_limits<float>::max(),
			-std::numeric_limits<float>::max()
		};

		for (const auto& p : points) {
			bbox.minX = std::min(bbox.minX, p[0]);
			bbox.minY = std::min(bbox.minY, p[1]);
			bbox.maxX = std::max(bbox.maxX, p[0]);
			bbox.maxY = std::max(bbox.maxY, p[1]);
		}
		return bbox;
	}

	static float overlap_ratio(const BBox& a, const BBox& b) {
		float intersectX = std::max(0.0f,
			std::min(a.maxX, b.maxX) - std::max(a.minX, b.minX));
		float intersectY = std::max(0.0f,
			std::min(a.maxY, b.maxY) - std::max(a.minY, b.minY));
		float intersectArea = intersectX * intersectY;
		float minArea = std::min(a.area(), b.area());
		return intersectArea / minArea;
	}
};


bool TreeIso::Init_seg_pcd(ccPointCloud* pc, const unsigned PR_MIN_NN1, const float PR_REG_STRENGTH1, const float PR_DECIMATE_RES1, const unsigned PR_THREAD1, QProgressDialog* progressDlg/*=nullptr*/)
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
	to_translated_vector(pc, pc_vec);

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

	const unsigned K = (PR_MIN_NN1 - 1);
	Vec3d edgeWeight;
	std::vector<uint32_t> Eu;
	std::vector<uint32_t> Ev;
	std::vector<uint32_t> in_component;
	std::vector<std::vector<uint32_t>> components;

	perform_cut_pursuit(K, 3, PR_REG_STRENGTH1, pc_sub, edgeWeight, Eu, Ev, in_component, PR_THREAD1);

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
		outSF->setValue(i, (int)in_component[ic[i]]);
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


bool TreeIso::Intermediate_seg_pcd(ccPointCloud* pc, const unsigned PR_MIN_NN2, const float PR_REG_STRENGTH2, const float PR_DECIMATE_RES2, const float PR_MAX_GAP, const unsigned PR_THREADS2, QProgressDialog* progressDlg/*=nullptr*/)
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

	std::vector<index_t> in_component;

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
	std::vector<Vec3d> pc_vec;
	to_translated_vector(pc, pc_vec);

	//std::vector<Vec3d> pc_vec;
	//load_initseg_points("F:\\treeiso\\data\\JP10_plot_2cm_test2_treeiso_dec_res.txt", pc_vec, in_component);
	//std::cout << "Loaded " << pc_vec.size() << " points" << std::endl;
	//unsigned pointCount = pc_vec.size();

	std::vector<std::vector<index_t>> clusterVGroup;
	unique_group(in_component, clusterVGroup);

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
	knn_cpp_nearest_neighbors(clusterCentroids, PR_MIN_NN2, minIdxsC, minIdxsD, PR_THREADS2);

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
		knn_cpp_build(knn_kdtree, PR_THREADS2);
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


	std::vector<std::vector<uint32_t>> minIdxs;
	std::vector<Vec3d> Ds;
	knn_cpp_nearest_neighbors(currentClusterDecsFlat, PR_MIN_NN2, minIdxs, Ds, PR_THREADS2);

	std::vector<float> edgeWeight;
	std::vector<index_t> Eu;
	std::vector<index_t> Ev;

	Eu.resize(minIdxs.size() + 1);
	Eu[0] = 0;

	Ev.clear();
	edgeWeight.clear();
	Ev.reserve(minIdxs.size() * minIdxs[0].size());
	edgeWeight.reserve(minIdxs.size() * minIdxs[0].size());

	for (size_t i = 0; i < minIdxs.size(); ++i)
	{
		size_t currentNode = currentClusterDecsFlatIndex[i];
		Vec3d currentDists = nnDists[currentNode];

		size_t edges_for_point = 0;
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
					//Eu.push_back(static_cast<index_t>(i));
					Ev.push_back(static_cast<index_t>(minIdxs[i][j]));
					edgeWeight.push_back(10 / ((nnDist + 0.001) / 0.01) * PR_REG_STRENGTH2);//when there is edge weight other than 1, this PR_REG_STRENGTH2 will be ignored; so multiply PR_REG_STRENGTH2 here
					edges_for_point++;
				}
			}
		}
		Eu[i + 1] = Eu[i] + static_cast<index_t>(edges_for_point);
	}

	if (progressDlg)
	{
		progressDlg->setValue(55);
		QCoreApplication::processEvents();
	}

	std::vector<index_t> in_component2d;
	perform_cut_pursuit(PR_MIN_NN2, 2, PR_REG_STRENGTH2, currentClusterDecsFlat, edgeWeight, Eu, Ev, in_component2d, 0);



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
	to_translated_vector(pc, pc_vec);

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

			std::vector<BBox> groupBBoxes;
			for (size_t i = 0; i < nGroups; ++i) {
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

				groupBBoxes.push_back(BBox::from_points(groupPts));
			}

			size_t knncpp_nn = (PR_MIN_NN3 < n_clusters ? PR_MIN_NN3 : n_clusters);
			std::vector<std::vector<uint32_t>> groupNNIdxC;
			std::vector<Vec3d> groupNNCDs;
			knn_cpp_nearest_neighbors(centroid2DFeatures, knncpp_nn, groupNNIdxC, groupNNCDs, 1);//threads=0 means optimal

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
			knn_cpp_build(knn_kdtree, 1);//threads=0 means optimal
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


				std::vector<float> scores;
				std::vector<size_t> filteredRemainIds;
				std::vector<float> min3DSpacings;

				for (size_t j = 0; j < knncpp_nn2; ++j)
				{
					size_t remainId = remainIds[groupNNIdx[i][j]];
					float lineSegs1 = zFeatures[remainId] + lenFeatures[remainId] - zFeatures[toMergeId];
					float lineSegs2 = zFeatures[toMergeId] + lenFeatures[toMergeId] - zFeatures[remainId];
					float verticalOverlapRatio = (lineSegs2 > lineSegs1 ? lineSegs1 : lineSegs2) / (lineSegs1 > lineSegs2 ? lineSegs1 : lineSegs2);
					float horizontalOverlapRatio = BBox::overlap_ratio(groupBBoxes[toMergeId], groupBBoxes[remainId]);

					Eigen::MatrixXf nnClusterCentroids;
					get_subset(clusterCentroids, clusterVGroup[remainId], nnClusterCentroids);

					knncpp::KDTreeMinkowskiX<float, knncpp::EuclideanDistance<float>> knn_kdtree2(nnClusterCentroids);
					knn_cpp_build(knn_kdtree2, 1);//threads=0 means optimal
					std::vector<std::vector<size_t>> min3D_idx;
					std::vector<Vec3d> min3D_dists;

					knn_cpp_query(knn_kdtree2, currentClusterCentroids, 1, min3D_idx, min3D_dists);
					float min3DSpacing = std::min_element(min3D_dists.begin(), min3D_dists.end(),
						[](const Vec3d& a, const Vec3d& b) { return a[0] < b[0]; })->operator[](0);

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
	return true;
}

//1. initial 3D segmentation
bool TreeIso::Init_seg(const unsigned PR_MIN_NN1,
	const float PR_REG_STRENGTH1,
	const float PR_DECIMATE_RES1,
	const unsigned PR_THREADS1,
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

	if (Init_seg_pcd(pointCloud, PR_MIN_NN1, PR_REG_STRENGTH1, PR_DECIMATE_RES1, PR_THREADS1, progressDlg))
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
								const unsigned PR_THREADS2,
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

	if (Intermediate_seg_pcd(pointCloud, PR_MIN_NN2, PR_REG_STRENGTH2, PR_DECIMATE_RES2, PR_MAX_GAP, PR_THREADS2, progressDlg))
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
	if (n_thread > 0) kdtree.setThreads(n_thread);
	kdtree.build();
}

void knn_cpp_query(knncpp::KDTreeMinkowskiX<float, knncpp::EuclideanDistance<float>>& kdtree, Eigen::MatrixXf& query_points, size_t k, std::vector<std::vector<size_t>>& res_idx, std::vector<Vec3d>& res_dists)
{
	res_idx.clear();
	res_dists.clear();
	knncpp::Matrixi indices;
	Eigen::MatrixXf distances;
	kdtree.query(query_points, k, indices, distances);

	// Change to column-wise access
	res_idx.resize(indices.cols(), std::vector<size_t>(indices.rows()));
	for (Eigen::Index i = 0; i < indices.cols(); i++)
	{
		res_idx[i] = std::vector<size_t>(indices.col(i).data(), indices.col(i).data() + indices.rows());
	}

	res_dists.resize(distances.cols(), Vec3d(distances.rows()));
	for (Eigen::Index i = 0; i < distances.cols(); i++)
	{
		res_dists[i] = Vec3d(distances.col(i).data(), distances.col(i).data() + distances.rows());
	}
}

float knn_cpp_query_min_d(knncpp::KDTreeMinkowskiX<float, knncpp::EuclideanDistance<float>>& kdtree, Eigen::MatrixXf& query_points, size_t k)
{
	knncpp::Matrixi indices;
	Eigen::MatrixXf distances;

	kdtree.query(query_points, k, indices, distances);
	return distances.minCoeff();
}

void knn_cpp_nearest_neighbors(const std::vector<Vec3d>& dataset,
	size_t k,
	std::vector<std::vector<uint32_t>>& res_idx,
	std::vector<Vec3d>& res_dists,
	unsigned n_thread) {

	if (dataset.empty()) {
		return;
	}

	Eigen::MatrixXf mat(dataset[0].size(), dataset.size());
	for (size_t i = 0; i < dataset.size(); i++) {
		mat.col(i) = Eigen::VectorXf::Map(&dataset[i][0], dataset[i].size());
	}

	knncpp::KDTreeMinkowskiX<float, knncpp::EuclideanDistance<float>> kdtree(mat);
	kdtree.setBucketSize(16);
	kdtree.setSorted(true);
	if (n_thread > 0) kdtree.setThreads(n_thread);
	kdtree.build();
	knncpp::Matrixi indices;
	Eigen::MatrixXf distances;

	kdtree.query(mat, k, indices, distances);

	res_idx.resize(indices.cols(), std::vector<uint32_t>(indices.rows()));
	for (Eigen::Index i = 0; i < indices.cols(); i++) {
		res_idx[i] = std::vector<uint32_t>(indices.col(i).data(), indices.col(i).data() + indices.rows());
	}
	res_dists.resize(distances.cols(), Vec3d(distances.rows()));
	for (Eigen::Index i = 0; i < distances.cols(); i++) {
		res_dists[i] = Vec3d(distances.col(i).data(), distances.col(i).data() + distances.rows());
	}
}

void build_knn_graph(const std::vector<Vec3d>& points, size_t k,
	std::vector<index_t>& first_edge,
	std::vector<index_t>& adj_vertices,
	std::vector<float>& edge_weights, float regStrength1, unsigned n_thread) {

	const size_t n_points = points.size();
	std::vector<std::vector<uint32_t>> nn_idx;
	std::vector<Vec3d> nn_D;

	// Get k+1 nearest neighbors (first one will be the point itself)
	knn_cpp_nearest_neighbors(points, k + 1, nn_idx, nn_D, n_thread);

	// Initialize edge arrays
	first_edge.resize(n_points + 1);
	first_edge[0] = 0;

	adj_vertices.clear();
	edge_weights.clear();
	adj_vertices.reserve(n_points * k);
	edge_weights.reserve(n_points * k);

	// Build graph structure
	for (size_t i = 0; i < n_points; i++) {
		size_t edges_for_point = 0;
		// Skip first neighbor (self)
		for (size_t j = 1; j < k + 1; j++) {
			adj_vertices.push_back(nn_idx[i][j]);
			// Convert distance to weight - for d0, we want high weights for close points
			float dist = std::sqrt(nn_D[i][j]) + 1e-6f;
			edge_weights.push_back(std::exp(-dist * dist) * regStrength1); // Gaussian weight
			edges_for_point++;
		}
		first_edge[i + 1] = first_edge[i] + static_cast<index_t>(edges_for_point);
	}
}

bool perform_cut_pursuit(const unsigned K,
	size_t D,
	const float regStrength,
	const std::vector<Vec3d>& pc_vec,
	std::vector<float>& edge_weights,
	std::vector<index_t>& Eu,
	std::vector<index_t>& Ev,
	std::vector<index_t>& in_component,
	const unsigned threads
) {

	using CP = Cp_d0_dist<float, index_t, comp_t>;

	const index_t pointCount = static_cast<index_t>(pc_vec.size());
	if (pointCount == 0) {
		return false;
	}

	if (Eu.size() == 0) {
		// Build graph structure using efficient kNN
		build_knn_graph(pc_vec, K, Eu, Ev, edge_weights, regStrength, threads);
	}

	const index_t E = static_cast<index_t>(Ev.size());

	// Convert point cloud to flat array for observations
	std::vector<float> Y(pointCount * D);
	for (size_t i = 0; i < pointCount; i++) {
		for (size_t d = 0; d < D; d++) {
			Y[i * D + d] = pc_vec[i][d];
		}
	}

	// Create cut pursuit instance
	CP* cp = new CP(pointCount, E, Eu.data(), Ev.data(), Y.data(), D);

	// Rest of the implementation stays the same
	cp->set_edge_weights(edge_weights.data(), regStrength);
	cp->set_loss(cp->quadratic_loss());
	cp->set_cp_param(1e-4, 20, 1000);
	//cp->set_min_comp_weight(10.0);//optional

	comp_t rV = 1;
	cp->set_components(rV, nullptr);
	cp->cut_pursuit();

	// Get components assignment
	const comp_t* comp_assign;
	const index_t* first_vertex;
	const index_t* comp_list;
	rV = cp->get_components(&comp_assign, &first_vertex, &comp_list);

	// Copy results, converting back to unsigned types for output
	in_component.resize(pointCount);
	for (uint32_t i = 0; i < pointCount; i++) {
		in_component[i] = static_cast<index_t>(comp_assign[i]);
	}

	delete cp;
	return true;
}

// Load points from file: an example
void load_initseg_points(const std::string& filename, std::vector<Vec3d>& points, std::vector<index_t>& in_component) {
	std::ifstream file(filename);
	if (!file.is_open()) {
		std::cerr << "Error opening file: " << filename << std::endl;
		return;
	}

	float x, y, z, s;
	while (file >> x >> y >> z >> s) {
		Vec3d point = { x, y, z };
		points.push_back(point);
		in_component.push_back(s);
	}
}

