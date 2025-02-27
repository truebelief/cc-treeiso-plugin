﻿#pragma once

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

class ccMainAppInterface;
class ccPointCloud;
class QProgressDialog;

class TreeIso
{
public:

	static bool Init_seg(const unsigned PR_MIN_NN1,const float PR_REG_STRENGTH1,const float PR_DECIMATE_RES1,const unsigned PR_THREADS1,ccMainAppInterface* app, std::function<void(int)> progressCallBack);
	static bool Intermediate_seg(const unsigned PR_MIN_NN2, const float PR_REG_STRENGTH2, const float PR_DECIMATE_RES2, const float PR_MAX_GAP, const unsigned PR_THREADS2, ccMainAppInterface* app, std::function<void(int)> progressCallBack);
	static bool Final_seg(const unsigned PR_MIN_NN3, const float PR_REL_HEIGHT_LENGTH_RATIO, const float PR_VERTICAL_WEIGHT, ccMainAppInterface* app, std::function<void(int)> progressCallBack);

	static bool Init_seg_pcd(ccPointCloud* pc, const unsigned PR_MIN_NN1, const float PR_REG_STRENGTH1, const float PR_DECIMATE_RES1, const unsigned PR_THREADS1, std::function<void(int)> progressCallback);
	static bool Intermediate_seg_pcd(ccPointCloud* pc, const unsigned PR_MIN_NN2, const float PR_REG_STRENGTH2, const float PR_DECIMATE_RES2, const float PR_MAX_GAP, const unsigned PR_THREADS2, std::function<void(int)> progressCallback);
	static bool Final_seg_pcd(ccPointCloud* pc, const unsigned PR_MIN_NN3, const float PR_REL_HEIGHT_LENGTH_RATIO, const float PR_VERTICAL_WEIGHT, std::function<void(int)> progressCallback);
};
