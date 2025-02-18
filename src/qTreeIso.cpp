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

#include "qTreeIso.h"

//Qt
#include <QApplication>
#include <QProgressDialog>
#include <QMainWindow>
#include <QComboBox>
#include <QElapsedTimer>
#include <QMessageBox>

//Local
#include "ccTreeIsoDlg.h"
#include "qTreeIsoCommands.h"

//System
#include <iostream>
#include <vector>
#include <string>
#include <assert.h>

//qCC_db
#include <ccGenericPointCloud.h>
#include <ccPointCloud.h>
#include <ccHObjectCaster.h>
#include <ccOctree.h>
#include <ccMesh.h>

//TreeIso
#include <TreeIso.h>

qTreeIso::qTreeIso(QObject* parent)
	: QObject(parent)
	, ccStdPluginInterface(":/CC/plugin/qTreeIso/info.json")
	, m_action(nullptr)
{
}	

void qTreeIso::onNewSelection(const ccHObject::Container& selectedEntities)
{
	if (m_action)
	{
		bool hasCloud = false;
		for (ccHObject* entity : selectedEntities)
		{
			if (entity && entity->isA(CC_TYPES::POINT_CLOUD))
			{
				hasCloud = true;
				break;
			}
		}
		m_action->setEnabled(hasCloud);
	}
}

QList<QAction *> qTreeIso::getActions()
{
	if (!m_action)
	{
		m_action = new QAction(getName(), this);
		m_action->setToolTip(getDescription());
		m_action->setIcon(getIcon());

		//connect appropriate signal
		connect(m_action, &QAction::triggered, this, &qTreeIso::doAction);
	}

	return { m_action };
}

void qTreeIso::doAction()
{
	Parameters parameters;
	ccTreeIsoDlg treeisoDlg(m_app ? m_app->getMainWindow() : nullptr);

	connect(treeisoDlg.pushButtonInitSeg, &QPushButton::clicked, [&]
		{
			parameters.min_nn1 = treeisoDlg.spinBoxK1->value();
			parameters.reg_strength1 = treeisoDlg.doubleSpinBoxLambda1->value();;
			parameters.decimate_res1 = treeisoDlg.doubleSpinBoxDecRes1->value();
			parameters.threads1 = treeisoDlg.doubleSpinBoxThreads1->value();
			
			init_segs(parameters, &treeisoDlg);
		});

	connect(treeisoDlg.pushButtonInterSeg, &QPushButton::clicked, [&]
		{
			parameters.min_nn2 = treeisoDlg.spinBoxK2->value();
			parameters.reg_strength2 = treeisoDlg.doubleSpinBoxLambda2->value();
			parameters.decimate_res2 = treeisoDlg.doubleSpinBoxDecRes2->value();
			parameters.max_gap = treeisoDlg.doubleSpinBoxMaxGap->value();
			parameters.threads2 = treeisoDlg.doubleSpinBoxThreads2->value();

			intermediate_segs(parameters, &treeisoDlg);

		});

	connect(treeisoDlg.pushButtonReseg, &QPushButton::clicked, [&]
		{
			parameters.rel_height_length_ratio = treeisoDlg.doubleSpinBoxRelHLRatio->value();
			parameters.vertical_weight = treeisoDlg.doubleSpinBoxVWeight->value();

			final_segs(parameters, &treeisoDlg);
		});
	treeisoDlg.pushButtonInitSeg->setEnabled(true);

	if (treeisoDlg.exec())
	{
		if (m_app)
		{
			m_app->refreshAll();
		}
		else
		{
			// m_app should have already been initialized by CC when plugin is loaded!
			assert(false);
		}
	}
}

void qTreeIso::init_segs(const Parameters& parameters, QWidget* parent/*=nullptr*/)
{
	// display the progress dialog
	QProgressDialog* progressDlg = new QProgressDialog(parent);
	progressDlg->setWindowTitle("TreeIso Step 1. Initial segmentation");
	progressDlg->setLabelText(tr("Computing...."));
	progressDlg->setCancelButton(nullptr);
	progressDlg->setRange(0, 100); // infinite progress bar
	progressDlg->show();
	
	std::function<void(int)> progressCallback = [progressDlg](int percent) {
		progressDlg->setValue(percent);
		QApplication::processEvents(); // Allow UI updates
	};

	if (!TreeIso::Init_seg(parameters.min_nn1, parameters.reg_strength1, parameters.decimate_res1, parameters.threads1, m_app, progressCallback))
	{
		progressCallback(0);
		m_app->dispToConsole("Discontinued in initial segmentation.", ccMainAppInterface::ERR_CONSOLE_MESSAGE);
		progressDlg->close();
		QApplication::processEvents();
		m_app->updateUI();
		m_app->refreshAll();
		return;
	}

	progressDlg->close();
	QApplication::processEvents();

	m_app->updateUI();
	m_app->refreshAll();
}

void qTreeIso::intermediate_segs(const Parameters& parameters, QWidget* parent/*=nullptr*/)
{
	// display the progress dialog
	QProgressDialog* progressDlg = new QProgressDialog(parent);
	progressDlg->setWindowTitle("TreeIso Step 2. Interim segmentation");
	progressDlg->setLabelText(tr("Computing...."));
	progressDlg->setCancelButton(nullptr);
	progressDlg->setRange(0, 100); // infinite progress bar
	progressDlg->show();
	
	std::function<void(int)> progressCallback = [progressDlg](int percent) {
		progressDlg->setValue(percent);
		QApplication::processEvents(); // Allow UI updates
	};
	
	if (!TreeIso::Intermediate_seg(parameters.min_nn2, parameters.reg_strength2, parameters.decimate_res2, parameters.max_gap, parameters.threads2, m_app, progressCallback))
	{
		progressCallback(0);
		m_app->dispToConsole("Error in interim segmentation!", ccMainAppInterface::ERR_CONSOLE_MESSAGE);
		progressDlg->close();
		QApplication::processEvents();
		m_app->updateUI();
		m_app->refreshAll();
		return;
	}

	progressDlg->hide();
	QApplication::processEvents();

	m_app->updateUI();
	m_app->refreshAll();
}

void qTreeIso::final_segs(const Parameters& parameters, QWidget* parent/*=nullptr*/)
{
	// display the progress dialog
	QProgressDialog* progressDlg = new QProgressDialog(parent);
	progressDlg->setWindowTitle("TreeIso Step 3. Final segmentation");
	progressDlg->setLabelText(tr("Computing...."));
	progressDlg->setCancelButton(nullptr);
	progressDlg->setRange(0, 100); // infinite progress bar
	progressDlg->show();
	
	std::function<void(int)> progressCallback = [progressDlg](int percent) {
		progressDlg->setValue(percent);
		QApplication::processEvents(); // Allow UI updates
	};

	if (!TreeIso::Final_seg(parameters.min_nn2, parameters.rel_height_length_ratio, parameters.vertical_weight, m_app, progressCallback))
	{
		progressCallback(0);
		m_app->dispToConsole("Error in final segmentation!", ccMainAppInterface::ERR_CONSOLE_MESSAGE);
		progressDlg->close();
		QApplication::processEvents();
		m_app->updateUI();
		m_app->refreshAll();
		return;
	}

	progressDlg->close();
	QApplication::processEvents();

	m_app->updateUI();
	m_app->refreshAll();
}

void qTreeIso::registerCommands(ccCommandLineInterface* cmd)
{
	if (!cmd)
	{
		assert(false);
		return;
	}
	cmd->registerCommand(ccCommandLineInterface::Command::Shared(new CommandTreeIso));
}
