# Changelog

All notable changes to this project will be documented in this file.

The format is based on [this branch](https://github.com/truebelief/cc-treeiso-plugin/tree/revert-7-optimize_treeiso_speed).


## - 2023-05-13
### Cleanup
- Adjust the Boost and Eigen library dependency on CMake

## - 2023-05-20
### Cleanup
- Major cleanup to integrate into CloudCompare

## - 2023-01-15
### Added
- Icon added and better progress dialog management

## - 2023-06-18
### Cleanup
- Remove two useless pragma once
- Use standard uint64 type
- Fix case of a file name in CMakeLists.txt

## - 2023-06-30
### Cleanup
- Fixes after a scan of the CloudCompare project by Coverity.

## - 2025-02-01
### Cleanup
- Removed the Boost dependency by adopting a newer cut-pursuit version (boosting speed by 5–15×).
- Fixed an issue with the nearest neighbor search caused by an incorrect 2D array order in knn_cpp.
- Enabled thread count customization for knn-cpp (minor speed boost).
- Optimized array calculations with improved use of STL algorithms.

## - 2025-02-10
### Tested
- Tested on several TLS plots

## - 2025-02-11
### Added
- Improved progress bar increment
- Added rapid detection of ground points based on quantile distribution
- Added warning of ground points before continuing initial segmentation
