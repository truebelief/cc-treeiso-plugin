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

// A Matlab version shared via:
// https://github.com/truebelief/artemis_treeiso

#include "ccStdPluginInterface.h"

class ccTreeIsoDlg;

//! TreeIso plugin
class qTreeIso : public QObject, public ccStdPluginInterface
{
	Q_OBJECT
	Q_INTERFACES( ccPluginInterface ccStdPluginInterface )
	
	Q_PLUGIN_METADATA( IID "cccorp.cloudcompare.plugin.qTreeISO" FILE "../info.json" )

public:

	//! Default constructor
	explicit qTreeIso(QObject* parent = nullptr);
	//! Destructor
	virtual ~qTreeIso() = default;

	//inherited from ccStdPluginInterface
	void onNewSelection(const ccHObject::Container& selectedEntities) override;
	QList<QAction *> getActions() override;
	void registerCommands(ccCommandLineInterface* cmd) override;

	//! Parameters
	struct Parameters
	{
		float reg_strength1 = 1.0f; //lambda1
		int min_nn1 = 5; //K1:key parameter
		float decimate_res1 = 0.05f;

		int reg_strength2 = 20; //lambda2:key parameter
		int min_nn2 = 20; //K2:key parameter
		float decimate_res2 = 0.1f;
		float max_gap = 2.0f;

		float rel_height_length_ratio = 0.5f; //rho
		float vertical_weight = 0.5; //w:key parameter
	};

protected:

	//! Slot called when associated ation is triggered
	void doAction();

	void init_segs(const Parameters& parameters, QWidget* parent = nullptr);
	void intermediate_segs(const Parameters& parameters, QWidget* parent = nullptr);
	void final_segs(const Parameters& parameters, QWidget* parent = nullptr);

protected: // members

	//! Associated action
	QAction* m_action;
};
