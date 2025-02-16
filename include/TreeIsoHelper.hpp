#pragma once

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

//Local
#include "knncpp.h"
#include "cp_d0_dist.hpp"

//Eigen
#include <Eigen/Dense>

//STL
#include <vector>
#include <functional>
#include <cstddef>
#include <unordered_map>

using namespace std;

class ccPointCloud;
typedef std::vector<float> Vec3d;

typedef uint32_t index_t;   // For vertex and edge indices
typedef uint16_t comp_t;    // For component indices


void knn_cpp_build(knncpp::KDTreeMinkowskiX<float, knncpp::EuclideanDistance<float>>& kdtree, unsigned n_thread = 0);
void knn_cpp_query(knncpp::KDTreeMinkowskiX<float, knncpp::EuclideanDistance<float>>& kdtree, Eigen::MatrixXf& query_points, size_t k, std::vector <std::vector<size_t>>& res_idx, std::vector <std::vector<float>>& res_dists);
float knn_cpp_query_min_d(knncpp::KDTreeMinkowskiX<float, knncpp::EuclideanDistance<float>>& kdtree, Eigen::MatrixXf& query_points, size_t k);
void build_knn_graph(const std::vector<Vec3d>& points, size_t k, std::vector<index_t>& first_edge, std::vector<index_t>& adj_vertices, std::vector<float>& edge_weights, float regStrength1 = 1.0, unsigned n_thread = 8);
void knn_cpp_nearest_neighbors(const std::vector<Vec3d>& dataset, size_t k, std::vector<std::vector<uint32_t>>& res_idx, std::vector<Vec3d>& res_dists, unsigned n_thread);
bool detectGroundByQuantile(const std::vector<std::vector<float>>& points, std::function<void(int)> progressCallBack, float groundThresholdFraction = 0.1f, float ratioThreshold = 0.1f);

void load_initseg_points(const std::string& filename, std::vector<Vec3d>& points, std::vector<index_t>& in_component);
bool perform_cut_pursuit(const unsigned K, size_t D, const float regStrength, const std::vector<Vec3d>& pc_vec, std::vector<float>& edge_weights, std::vector<index_t>& Eu, std::vector<index_t>& Ev, std::vector<index_t>& in_component, const unsigned threads, std::function<void(int)> progressCallback);

template <typename T>
size_t arg_min_col(const std::vector<T>& arr) {
	return std::distance(arr.begin(), std::min_element(arr.begin(), arr.end()));
}

template <typename T>
size_t arg_max_col(const std::vector<T>& arr) {
	return std::distance(arr.begin(), std::max_element(arr.begin(), arr.end()));
}

template <typename T>
void min_col(const std::vector<std::vector<T>>& arr, std::vector<T>& min_vals) {
	if (arr.empty()) {
		min_vals.clear();
		return;
	}

	min_vals = arr[0];
	for (const auto& row : arr) {
		std::transform(min_vals.begin(), min_vals.end(), row.begin(),
			min_vals.begin(), [](const T& a, const T& b) { return std::min(a, b); });
	}
}

template <typename T>
T min_col(const std::vector<T>& arr) {
	return arr.empty() ? std::numeric_limits<T>::quiet_NaN()
		: *std::min_element(arr.begin(), arr.end());
}

template <typename T>
T mean_col(const std::vector<T>& arr) {
	if (arr.empty()) return std::numeric_limits<T>::quiet_NaN();
	return static_cast<T>(std::accumulate(arr.begin(), arr.end(), 0.0) / arr.size());
}

template <typename T>
T median_col(std::vector<T>& arr) {
	if (arr.empty()) return std::numeric_limits<T>::quiet_NaN();

	const size_t n = arr.size();
	const size_t mid = n / 2;
	std::nth_element(arr.begin(), arr.begin() + mid, arr.end());

	if (n % 2 == 0) {
		const T right = arr[mid];
		std::nth_element(arr.begin(), arr.begin() + mid - 1, arr.end());
		return (arr[mid - 1] + right) / 2;
	}
	return arr[mid];
}

template <typename T>
T mode_col(const std::vector<T>& arr) {
	if (arr.empty()) return std::numeric_limits<T>::quiet_NaN();

	std::unordered_map<T, size_t> freq;
	for (const auto& val : arr) ++freq[val];
	return std::max_element(freq.begin(), freq.end(),
		[](const auto& a, const auto& b) { return a.second < b.second; })->first;
}

template <typename T>
void max_col(const std::vector<std::vector<T>>& arr, std::vector<T>& max_vals) {
	if (arr.empty()) {
		max_vals.clear();
		return;
	}

	max_vals = arr[0];
	for (const auto& row : arr) {
		std::transform(max_vals.begin(), max_vals.end(), row.begin(),
			max_vals.begin(), [](const T& a, const T& b) { return std::max(a, b); });
	}
}

template <typename T>
void mean_col(const std::vector<std::vector<T>>& arr, std::vector<T>& mean_vals) {
	if (arr.empty()) {
		mean_vals.clear();
		return;
	}

	const size_t cols = arr[0].size();
	mean_vals.resize(cols);
	std::fill(mean_vals.begin(), mean_vals.end(), T{});

	for (const auto& row : arr) {
		std::transform(mean_vals.begin(), mean_vals.end(), row.begin(),
			mean_vals.begin(), std::plus<T>());
	}

	const T size = static_cast<T>(arr.size());
	std::transform(mean_vals.begin(), mean_vals.end(), mean_vals.begin(),
		[size](T val) { return val / size; });
}

template <typename T>
void decimate_vec(const std::vector<std::vector<T>>& arr, T res, std::vector<std::vector<T>>& vec_dec) {
	if (arr.empty() || res <= T{}) {
		vec_dec.clear();
		return;
	}

	std::vector<T> arr_min;
	min_col(arr, arr_min);

	vec_dec.resize(arr.size(), std::vector<T>(arr[0].size()));
	for (size_t i = 0; i < arr.size(); ++i) {
		std::transform(arr[i].begin(), arr[i].end(), arr_min.begin(),
			vec_dec[i].begin(),
			[res](T val, T min) { return std::floor((val - min) / res) + T{ 1 }; });
	}
}


template <typename T>
void sort_indexes_by_row(const std::vector<std::vector<T>>& v, std::vector<size_t>& idx,
	std::vector<std::vector<T>>& v_sorted) {
	if (v.empty()) {
		idx.clear();
		v_sorted.clear();
		return;
	}

	const size_t rows = v.size();
	const size_t cols = v[0].size();

	idx.resize(rows);
	std::iota(idx.begin(), idx.end(), 0);

	std::stable_sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {
		return std::lexicographical_compare(v[i1].begin(), v[i1].end(),
			v[i2].begin(), v[i2].end());
		});

	v_sorted.resize(rows);
	for (size_t i = 0; i < rows; ++i) {
		v_sorted[i] = v[idx[i]];
	}
}

template <typename ValueType, typename IndexType>
void sort_indexes(const std::vector<ValueType>& v, std::vector<IndexType>& idx,
	std::vector<ValueType>& v_sorted) {
	if (v.empty()) {
		idx.clear();
		v_sorted.clear();
		return;
	}

	idx.resize(v.size());
	std::iota(idx.begin(), idx.end(), 0);

	std::stable_sort(idx.begin(), idx.end(),
		[&v](IndexType i1, IndexType i2) { return v[i1] < v[i2]; });

	v_sorted.resize(v.size());
	std::transform(idx.begin(), idx.end(), v_sorted.begin(),
		[&v](IndexType i) { return v[i]; });
}


template <typename T>
struct VectorHash {
	std::size_t operator()(const std::vector<T>& vec) const {
		std::size_t seed = vec.size();
		for (const auto& elem : vec) {
			seed ^= std::hash<T>{}(elem)+0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
		}
		return seed;
	}
};
template <typename T>
void unique_index_by_rows(const std::vector<std::vector<T>>& arr,
	std::vector<size_t>& ia,
	std::vector<size_t>& ic) {
	if (arr.empty()) {
		ia.clear();
		ic.clear();
		return;
	}

	ia.clear();
	ic.resize(arr.size());
	// Map each unique row to its unique group index.
	std::unordered_map<std::vector<T>, size_t, VectorHash<T>> row_to_group;

	size_t unique_counter = 0;
	for (size_t i = 0; i < arr.size(); ++i) {
		auto it = row_to_group.find(arr[i]);
		if (it == row_to_group.end()) {
			// First occurrence of this unique row.
			row_to_group[arr[i]] = unique_counter;
			ia.push_back(i);
			ic[i] = unique_counter;
			++unique_counter;
		}
		else {
			ic[i] = it->second;
		}
	}
}

template <typename T>
void to_translated_vector(const ccPointCloud* pc, std::vector<std::vector<T>>& y) {
	if (!pc || pc->size() == 0) {
		y.clear();
		return;
	}

	const size_t pointCount = pc->size();
	y.resize(pointCount, std::vector<T>(3));

	std::vector<T> y_mean(3, 0);
	for (size_t i = 0; i < pointCount; ++i) {
		const CCVector3* pv = pc->getPoint(i);
		y[i] = { static_cast<T>(pv->x), static_cast<T>(pv->y), static_cast<T>(pv->z) };
		std::transform(y_mean.begin(), y_mean.end(), y[i].begin(), y_mean.begin(), std::plus<T>());
	}

	std::transform(y_mean.begin(), y_mean.end(), y_mean.begin(),
		[pointCount](T val) { return val / pointCount; });

	for (auto& point : y) {
		std::transform(point.begin(), point.end(), y_mean.begin(), point.begin(), std::minus<T>());
	}
}


template <typename T>
void unique_group(const std::vector<T>& arr, std::vector<std::vector<T>>& u_group,
	std::vector<T>& arr_unq, std::vector<T>& ui) {
	if (arr.empty()) {
		arr_unq.clear();
		ui.clear();
		u_group.clear();
		return;
	}

	std::vector<T> arr_sorted_idx;
	std::vector<T> arr_sorted;
	sort_indexes(arr, arr_sorted_idx, arr_sorted);

	arr_unq.clear();
	ui.clear();
	u_group.clear();

	ui.push_back(arr_sorted_idx[0]);
	std::vector<T> current_group = { arr_sorted_idx[0] };

	for (size_t i = 1; i < arr.size(); ++i) {
		if (arr_sorted[i] != arr_sorted[i - 1]) {
			ui.push_back(arr_sorted_idx[i]);
			arr_unq.push_back(arr_sorted[i - 1]);
			u_group.push_back(std::move(current_group));
			current_group = { arr_sorted_idx[i] };
		}
		else {
			current_group.push_back(arr_sorted_idx[i]);
		}
	}

	arr_unq.push_back(arr_sorted.back());
	u_group.push_back(std::move(current_group));
}

// Overloaded versions with fewer return parameters
template <typename T>
void unique_group(const std::vector<T>& arr, std::vector<std::vector<T>>& u_group,
	std::vector<T>& arr_unq) {
	std::vector<T> ui;
	unique_group(arr, u_group, arr_unq, ui);
}

template <typename T>
void unique_group(const std::vector<T>& arr, std::vector<std::vector<T>>& u_group) {
	std::vector<T> arr_unq, ui;
	unique_group(arr, u_group, arr_unq, ui);
}

template <typename T1, typename T2>
void get_subset(const std::vector<T1>& arr,const std::vector<T2>& indices,std::vector<T1>& arr_sub) {
	arr_sub.resize(indices.size());
	for (size_t i = 0; i < indices.size(); ++i) {
		arr_sub[i] = arr[indices[i]];
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
			arr_sub(j, i) = arr[indices[i]][j];
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
