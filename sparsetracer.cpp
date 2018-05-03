/* SparseTracer
*	function: trace single neuron
*	Author : zhouhang, lishiwei
*	2015-10-28
*/
#include <QFileDialog>
#include <QMessageBox>
#include <iostream>
#include <QDir>
#include <QDateTime>
#include <QDialog>
#ifdef _WIN32
#include <io.h>
#include <direct.h>
#include <ctime>
#include <fstream>
#define GLOG_NO_ABBREVIATED_SEVERITIES
#include <glog/logging.h>
using namespace google;
#pragma warning(disable : 4996)
#else
#include <unistd.h>
#include <sys/stat.h>
#include<sys/types.h>
#include<dirent.h>
#endif
#include "sparsetracer.h"
#include "ui_sparsetracer.h"
#include "DockWidget/binarysettingsdockwidget.h"
#include "Dialog/openimagedialog.h"
#include "Function/IO/imagereader.h"
#include "Function/IO/somareader.h"
#include "Render/neuroglwidget.h"
#include "Function/binaryfilter.h"
#include "ngtypes/soma.h"
#include "Function/IO/treeReader.h"
#include "Function/IO/somawriter.h"
#include "ngtypes/tree.h"
#include "Function/IO/treewriter.h"
#include "Function/Trace/SparseTrace/sparsetracefilter.h"
#include "Function/Trace/SparseTrace/largesparsetracefilter.h"
#include "Function/IO/imagewriter.h"
#include "Function/IO/MOSTDReader.h"
#include "Function/IO/HDF5Reader.h"
#include "Function/IO/HDF5Writer.h"
#include "Function/volumealgo.h"
#include "../ngtypes/ParamPack.h"
#include "Function/IO/tracestackreader.h"
#include "Function/IO/tracestackwriter.h"
#include "Function/Trace/traceutil.h"
#include "Function/Trace/NeuroGPSTreeFilter.h"
#include "Function/Trace/SparseTrace/sparsetraceaxonfilter.h"
#include "Function/DendriteLevelImageMaker.h"
#include "Function/NGUtility.h"
#include "Function/probefilter.h"
#include "Function/ImageScaleFiter.h"
#include "TreeChecker.h"

#define QQUESTION(title, tip) QMessageBox::question(this, title, tip, QMessageBox::Yes | QMessageBox::No, QMessageBox::No)

void SignalHandle(const char* data, int size)
{
	std::ofstream fs("glog_dump.log", std::ios::app);
	std::string str = std::string(data, size);
	fs << str;
	fs.close();
	LOG(ERROR) << str;
}

SparseTracer::SparseTracer(QWidget *parent) :
QMainWindow(parent),
ui(new Ui::SparseTracer)
{
	paramPack = NGParamPack(new ParamPack());

	ui->setupUi(this);
	this->resize(QSize(1024, 720));

	
	sumSteps_ = 0;
	manualLabelCurve_ = SMARTTREE(new std::vector<VectorVec5d>);
	//initialize glog
	google::InitGoogleLogging("test program");
	QString strLogPath = QDir::currentPath() + "/LogInfo/";
	QDir dir; dir.mkpath(strLogPath);
	google::SetLogDestination(google::GLOG_INFO, strLogPath.toStdString().c_str());//
	FLAGS_stderrthreshold = 3;

	oldEditToolbarTag = newEditToolbarTag = EDITTOOLBARTAG::NONE;
	toolBarButtonGroup = new QActionGroup(this);
	toolBarButtonGroup->addAction(ui->actionChoose_Line);
	toolBarButtonGroup->addAction(ui->actionChoose_Vertex);
	toolBarButtonGroup->addAction(ui->actionDraw_Line);
	toolBarButtonGroup->addAction(ui->actionPick_Soma);
	toolBarButtonGroup->addAction(ui->actionSelect_Tree);
	ui->actionRun->setEnabled(false);
	ui->actionStop->setEnabled(false);
	ui->actionBackTrack->setEnabled(false);
	ui->action2DView->setEnabled(false);
	ui->actionChoose_Line->setEnabled(false);
	ui->actionChoose_Vertex->setEnabled(false);
	ui->actionCut_Vertex->setEnabled(false);
	ui->actionDelete_Line->setEnabled(false);
	ui->actionDelete_Soma->setEnabled(false);
	ui->actionDelete_Tree->setEnabled(false);
	ui->actionDraw_Line->setEnabled(false);
	ui->actionNGTree_Trace->setEnabled(false);
	ui->actionPick_Soma->setEnabled(false);
	ui->actionSelect_Tree->setEnabled(false);
	ui->actionTrain->setEnabled(false);
	ui->actionVisible->setChecked(true);
	ui->actionZoom_in->setEnabled(false);
	ui->actionZoom_out->setEnabled(false);
	ui->actionTreeChecker->setEnabled(false);
	//ui->actionTreeChecker->setEnabled(false);
	settingDock = new BinarySettingsDockWidget(this,paramPack);
	settingDock->SetParam(paramPack);
	settingDock->UpdateParam();
	this->addDockWidget(Qt::LeftDockWidgetArea, settingDock, Qt::Vertical);
	settingDock->setMaximumWidth(250);
	paramPack->separateTree = std::make_shared<NeuronPopulation>();
	glbox = new NeuroGLWidget(this);
	glbox->SetInputNeuronPopulation(paramPack->separateTree);
	glbox->SetParam(paramPack);
	glbox->setMinimumSize(100, 100);
	glbox->show();
	glbox->SetClippingDistance_Slot(paramPack->thicknesss_);
	ui->centralWidget->layout()->setMargin(0);
	ui->centralWidget->layout()->addWidget(glbox);
	ui->toolBar_2->setVisible(false);

	QFont statusFont("Arial", 12, QFont::Bold);
	ui->statusBar->setFont(statusFont);
	paramPack->defaultDir = ".";

	timer = new QTimer(this);
	connect(timer, SIGNAL(timeout()), this, SLOT(TimingSave_Slot()));
	runningTimer_ = new QTimer(this);
	connect(runningTimer_, SIGNAL(timeout()), this, SLOT(AddRunningTime_Slot()));
	connect(glbox, SIGNAL(BoxWidgetChanged_Signal()), this, SLOT(SetROIByBoxWidget_Slot()));
	connect(glbox, SIGNAL(ChangeMoveCursor_Signal(bool)), this, SLOT(ChangeMoveCursor_Slot(bool)));
	connect(glbox, SIGNAL(clippingDistanceChanged_Signal(int)), settingDock, SLOT(SetThickness_Slot(int)));
	connect(glbox, SIGNAL(Set2DProject_Signal(int)), this, SLOT(SetProjectMaxThickness_Slot(int)));
	connect(glbox, SIGNAL(Set2DView_Signal(bool)), this, SLOT(Set2DViewStatus_Slot(bool)));
	connect(glbox, SIGNAL(SetDraw_Signal(bool)), this, SLOT(SetDrawStatus_Slot(bool)));
	connect(settingDock, SIGNAL(ApplyReadImage_Signal()), this, SLOT(ApplyReadNewImage_Slot()));
	connect(settingDock, SIGNAL(DisplayNewOpac_Signal()), glbox, SLOT(UpdateVolumeOpac_Slot()));
	connect(settingDock, SIGNAL(MaxBoundNumChanged_Signal(int)), this, SLOT(MaxBoundNumChanged_Slot(int)));
	connect(settingDock, SIGNAL(ROIChanged_Signal()), this, SLOT(UpdateGLBoxWidget_Slot()));
	connect(settingDock, SIGNAL(ThicknessChanged_Signal(int)), glbox, SLOT(SetClippingDistance_Slot(int)));
	connect(toolBarButtonGroup, SIGNAL(triggered(QAction*)), this, SLOT(EditActionGroup_triggered(QAction*)));
	connect(settingDock, SIGNAL(ClearCache_Signal()), this, SLOT(ClearMOSTDCache_Slot()));
	//tree widget
	connect(settingDock, SIGNAL(ActiveTree_Signal(int)), this, SLOT(ActivateTree_Slot(int)));
	connect(settingDock, SIGNAL(CompareTree_Signal(int)), this, SLOT(CompareTree_Slot(int)));
	connect(settingDock, SIGNAL(ToggleTreeVisible_Signal(int)), this, SLOT(ToggleTreeVisible_Slot(int)));
	connect(settingDock, SIGNAL(GotoDiff_Signal(int)), this, SLOT(GotoDiff_Slot(int)));//TODO

	connect(settingDock, SIGNAL(ToggleTraverse_Signal(bool)), this, SLOT(ToggleTraverse_Slot(bool)));
	//connect(settingDock, SIGNAL(Annotate_Signal()), this, SLOT(Annotate_Slot()));
	connect(settingDock, SIGNAL(BackTraverse_Signal()), this, SLOT(BackTraverse_Slot()));
	connect(settingDock, SIGNAL(NextTraverse_Signal()), this, SLOT(NextTraverse_Slot()));
	connect(settingDock, SIGNAL(ToggleShowDiffArrow_Signal(bool)), this, SLOT(ToggleShowDiffArrow_Slot(bool)));

	connect(glbox, SIGNAL(SetTraversePos_Signal(int, int)), this, SLOT(StartTraverseFromHere_Slot(int, int)));
	connect(settingDock, SIGNAL(ResetTraverse_Signal()), this, SLOT(ResetTraverse_Slot()));

	connect(glbox, SIGNAL(NeedUpdateTreeWidget()), settingDock, SLOT(UpdateTreeWidget_Slot()));
	connect(glbox, SIGNAL(NeedUpdateCheckWidget()), settingDock, SLOT(UpdateCheckWidget_Slot()));

	connect(settingDock, SIGNAL(ToggleShowTraverseFlag_Signal()), glbox, SLOT(ToggleShowTraverseFlag_Slot()));

	connect(glbox, SIGNAL(NextTraverse_Signal()), this, SLOT(NextTraverse_Slot()));

	connect(settingDock, SIGNAL(OpacAdjustApply_Signal()), this, SLOT(OpacAdjustApply_Slot()));
	connect(settingDock, SIGNAL(PreviewOpacAdjust_Signal(bool)), this, SLOT(PreviewOpacAdjust_Slot(bool)));
	connect(settingDock, SIGNAL(PreviewBinary_Signal(bool)), this, SLOT(PreviewBinary_Slot(bool)));
	connect(settingDock, SIGNAL(PreviewAxonBinary_Signal(bool)), this, SLOT(PreviewAxonBinary_Slot(bool)));

	//connect(this, SIGNAL(readyToHis()), settingDock, SLOT(UpdateHistogram()));//2018-2-10
	connect(this, SIGNAL(renderCaliber()), glbox, SLOT(renderCaliber_Slot()));

	ui->actionTest->setVisible(true);
}

SparseTracer::~SparseTracer()
{
	delete ui;
	google::ShutdownGoogleLogging();
}

void SparseTracer::on_actionOpenImage_triggered()
{
	if (traverseMode_) {
		printf("traverse mode forbid open file.\n");
		return;
	}
	OpenImageDialog fileDialog(this);
	fileDialog.SetParam(paramPack);
	QString stackFileName;
	NGTraceStackReader stackReader;
	if (fileDialog.exec() == QDialog::Accepted){
		ClearAll();
		imageFileName_ = fileDialog.GetOpenImageName();
		if (imageFileName_.isEmpty()){
			QMessageBox::warning(this, "wrong image file!", "please input right image file.");
			return;
		}
		stackFileName = fileDialog.GetOpenStackName();
	}
	else return;
	if (paramPack->dataMode_ == READMODE::MOSTD) {
		mostdReader = MOSTDReader::New();
		mostdTraverseReader = MOSTDReader::New();
		stackReader = TraceStackReader::New();
		mostdReader->SetParam(paramPack);
		mostdTraverseReader->SetParam(paramPack);
		if (!mostdReader->SetInputFileName(imageFileName_.toStdString())) {
			printf("MOSTD file is invalid.\n");
			return;
		}
		if (!mostdTraverseReader->SetInputFileName(imageFileName_.toStdString())) {
			printf("MOSTD Traverse file is invalid.\n");
			return;
		}
	}
	else if (paramPack->dataMode_ == READMODE::HDF5){
		mostdReader = HDF5Reader::New();
		mostdTraverseReader = HDF5Reader::New();
		stackReader = TraceStackReader::New();
		mostdReader->SetParam(paramPack);
		mostdTraverseReader->SetParam(paramPack);
		if (!mostdReader->SetInputFileName(imageFileName_.toStdString())) {
			printf("MOSTD file is invalid.\n");
			return;
		}
		if (!mostdTraverseReader->SetInputFileName(imageFileName_.toStdString())) {
			printf("MOSTD Traverse file is invalid.\n");
			return;
		}
	}
	else{
		paramPack->xScale_ = SetScale(paramPack->xRes_);
		paramPack->yScale_ = SetScale(paramPack->yRes_);
		settingDock->UpdateImageScale();
	}


	glbox->SetParam(paramPack);
	//GUI
	if (paramPack->dataMode_ == READMODE::MOSTD) {
		paramPack->highOpac_ = mostdReader->GetImageType() == DATATYPE::IMAGE8 ? 255 : 1000;
	}
	else
		paramPack->highOpac_ = fileDialog.GetDataType() == DATATYPE::IMAGE8 ? 255 : 1000;

	settingDock->SetImageRange(paramPack->xRangeMax_, paramPack->yRangeMax_, paramPack->zRangeMax_);
	settingDock->SetApplyReadImageButton(true);
	settingDock->SetDisplayOpacValue(0, paramPack->highOpac_);
	ui->action2DView->setEnabled(true);
	ui->actionBackTrack->setEnabled(true);
	ui->actionChoose_Line->setEnabled(true);
	ui->actionChoose_Vertex->setEnabled(true);
	ui->actionNGTree_Trace->setEnabled(true);
	ui->actionPick_Soma->setEnabled(true);
	ui->actionSelect_Tree->setEnabled(true);
	ui->actionTrain->setEnabled(true);
	ui->actionVisible->setEnabled(true);
	ui->actionTrain->setEnabled(true);
	if (paramPack->dataMode_ == READMODE::MOSTD){
		ui->actionRun->setEnabled(true);
		ui->actionStop->setEnabled(true);
	}
	settingDock->UpdateGUI();

	//
	if (!stackFileName.isEmpty()){
		QStringList swcList = stackFileName.split(';');
		if (paramPack->dataMode_ == READMODE::MOSTD && swcList[0].section('.', -1).compare("sta") == 0) {//read stack file
			stackReader->SetInput(paramPack->separateTree);
			stackReader->SetParam(paramPack);
			stackReader->SetInputFileName(swcList[0].toStdString());
			auto res = stackReader->Update();
			if (!res->success()) {
				NG_ERROR_MESSAGE("cannot read stack file.");
				QMessageBox::warning(this, QString::fromStdString(res->ErrorClassName()), QString::fromStdString(res->ErrorInfo()));
				NG_CREATE_DYNAMIC_CAST(NeuronPopulation, tmp2, paramPack->separateTree);
				tmp2->Clear();
				return;
			}
			ui->actionTreeChecker->setEnabled(true);//need tree
		}
		else if (swcList[0].section('.', -1).compare("swc") == 0){//read SWC
			NGTreeReader treeReader = TreeReader::New();
			std::vector<std::string> swcFile;
			for (int k = 0; k < swcList.size(); ++k) {
				if (swcList[k].isEmpty()) continue;
				swcFile.push_back(swcList[k].toStdString());
			}
			treeReader->SetInputFileNameList(swcFile);
			treeReader->SetParam(paramPack);
			auto res = treeReader->Update();
			if (!res->success()){
				NG_ERROR_MESSAGE("tree reader error");
				QMessageBox::warning(this, QString::fromStdString(res->ErrorClassName()), QString::fromStdString(res->ErrorInfo()));
				return;
			}
			auto tmp1 = treeReader->ReleaseData();//seperateTree = is wrong, it is value transfer.
			NG_CREATE_DYNAMIC_CAST(NeuronPopulation, swcReadRes, tmp1);
			NG_CREATE_DYNAMIC_CAST(NeuronPopulation, tmp2, paramPack->separateTree);
			tmp2->m_pop.swap(swcReadRes->m_pop);
			ui->actionTreeChecker->setEnabled(true);//need tree
		}
		else if (swcList[0].section('.', -1).compare("cbd") == 0)//read cbd
		{
			if (!paramPack->caliberBulidData)
				paramPack->caliberBulidData = std::make_shared<CellVec3d>();
			else {
				paramPack->caliberBulidData->clear();
				paramPack->caliberBulidData = std::make_shared<CellVec3d>();
			}

			for (int i = 0; i < swcList.size(); ++i){
				if (swcList[i].isEmpty())
					continue;
				else{
					FILE * inputFile = fopen(swcList[i].toStdString().c_str(), "r");
					if (!inputFile){
						std::cout << "can not open file:" << std::endl << swcList[i].toStdString() << std::endl;
						break;
					}

					char line[256];
					int id, type, pid, indexNum = 0, falg = 0;
					VectorVec3d temVector;

					while (!feof(inputFile)) {
						Vec3d temPoint;

						if (fgets(line, 256, inputFile) == 0) 
							continue;
						sscanf(line, "%d %d %lf %lf %lf %d %d", &id, &type, &temPoint[0], &temPoint[1], &temPoint[2], &indexNum, &pid);
						temPoint[0] /= paramPack->xRes_;
						temPoint[1] /= paramPack->yRes_;
						temPoint[2] /= paramPack->zRes_;
						if (falg != indexNum){
							paramPack->caliberBulidData->push_back(temVector);
							temVector.clear();
							temVector.push_back(temPoint);
							falg = indexNum;
						}
						else temVector.push_back(temPoint);
					}
					paramPack->caliberBulidData->push_back(temVector);
					fclose(inputFile);
					std::cout << "read completed!" << std::endl << swcList[i].toStdString() << std::endl;

					//read completed,emit signal
					emit renderCaliber();
				}
			}
		}
	}
	else if (paramPack->dataMode_ == READMODE::MOSTD) {//no sta or swc file
		printf("Image have been Read, please set the trace start points.\n");
		return;
	}
	if (paramPack->dataMode_ == READMODE::TIFF3D){
		paramPack->xMin_ = paramPack->yMin_ = paramPack->zMin_ = 0;
		paramPack->xMax_ = paramPack->xRangeMax_;
		paramPack->yMax_ = paramPack->yRangeMax_;
		paramPack->zMax_ = paramPack->zRangeMax_;
	}

	printf("Image and Stack information Reading complete.\n");
	//-------------------display-----------//
	if (paramPack->xMax_ - paramPack->xMin_ > 19 && paramPack->yMax_ - paramPack->yMin_ > 19 && paramPack->zMax_ - paramPack->zMin_ > 9
		&& paramPack->xMin_ > -1 && paramPack->yMin_ > -1 && paramPack->zMin_ > -1){
		QApplication::setOverrideCursor(Qt::WaitCursor);
		if (!ReadImage()) {
			QApplication::restoreOverrideCursor(); return;
		}

		//emit readyToHis();

		printf("Image read complete.\n");
		QApplication::restoreOverrideCursor();
		glbox->SetInputVolume(paramPack->OrigImage, paramPack->OrigImage->DataType());
		glbox->UpdateBoxWidget_Slot();
	}
	else {
		QMessageBox::warning(this, "data size error", "Please check data size !");
	}
	glbox->SetInputNeuronPopulation(paramPack->separateTree);
	glbox->SetClippingDistance_Slot(paramPack->thicknesss_);
	settingDock->UpdateGUI();//again to refresh the result
	settingDock->UpdateTreeWidget();
	if (paramPack->activeTree) {
		if (!paramPack->activeTree->m_suspiciousPositions_.empty()) {
			settingDock->UpdateCheckWidget();
		}
	}
}

void SparseTracer::ApplyReadNewImage_Slot()
{
	if (true == traverseMode_ || (paramPack->dataMode_ != READMODE::MOSTD && paramPack->dataMode_ != READMODE::HDF5)){
		QMessageBox::information(this, "info", "current mode does not support local image.");
		return;
	}
	//get setting box ROI value
	settingDock->GetROIValues(boxWidgetBound);
	paramPack->xMin_ = boxWidgetBound[0]; paramPack->xMax_ = boxWidgetBound[1];
	paramPack->yMin_ = boxWidgetBound[2]; paramPack->yMax_ = boxWidgetBound[3];
	paramPack->zMin_ = boxWidgetBound[4]; paramPack->zMax_ = boxWidgetBound[5];
	if (!CheckReadImageRangeValid()) return;
	QApplication::setOverrideCursor(Qt::WaitCursor);
	if (!ReadImage()) {
		QApplication::restoreOverrideCursor(); return;
	}
	QApplication::restoreOverrideCursor();
	printf("Image Reading complete.\n");
	//-------------------display-----------//
	auto &layerDisplay_ = glbox->GetLayerDisplay();

	settingDock->SetMaxThickness(paramPack->zMax_ - paramPack->zMin_);
	glbox->SetParam(paramPack);
	if (!layerDisplay_.empty() && layerDisplay_[0] > 0) {
		DendriteLevelImageMaker maker;
		maker.SetParam(paramPack);
		maker.SetLayer(layerDisplay_); maker.SetBranch(glbox->GetBranchDisplay());
		maker.SetLayerRadius(layerRadiusDisplay_);
		maker.SetHasSoma(hasSoma_);
		if (!maker.Update()) {
			NG_ERROR_MESSAGE("Layer Image Construct failed.");
			glbox->SetInputVolume(paramPack->OrigImage, paramPack->OrigImage->DataType());
		}
		else{
			glbox->SetInputVolume(paramPack->LayerImage, paramPack->LayerImage->DataType());
			glbox->ToggleActiveTreeLayerDisplay(true);
		}
	}
	else{
		glbox->SetInputVolume(paramPack->OrigImage, paramPack->OrigImage->DataType());
	}

	glbox->BuildCurrentBoundaryCollection();
	glbox->UpdateBoxWidget_Slot();
	glbox->UpdateVolumeOpac_Slot();
	glbox->update();
	ui->action2DView->setEnabled(true);
	ui->actionChoose_Line->setEnabled(true);
	ui->actionChoose_Vertex->setEnabled(true);
	ui->actionVisible->setEnabled(true);
	if (paramPack->dataMode_ == READMODE::MOSTD) {
		if (paramPack->mostdLevel_ != 1) {
			ui->actionRun->setEnabled(false);
			ui->actionBackTrack->setEnabled(false);
			ui->actionStop->setEnabled(false);
		}
		else if (ui->actionRun->isEnabled() == false && paramPack->mostdLevel_ == 1){
			ui->actionRun->setEnabled(true);
			ui->actionBackTrack->setEnabled(true);
			ui->actionStop->setEnabled(true);
		}
	}
}

bool SparseTracer::CheckReadImageRangeValid()
{
	if (paramPack->xMax_ > 1000000 || paramPack->yMax_ > 1000000 || paramPack->zMax_ > 1000000){
		QMessageBox::warning(this, "error reader", "you may not read most d or the image reader is broken\n");
		return false;
	}
	if (paramPack->xMax_ <= paramPack->xMin_ || paramPack->yMax_ <= paramPack->yMin_ || paramPack->zMax_ <= paramPack->zMin_) {
		QMessageBox::warning(this, "warning", "User-defined image range is out of original image range. Please set suitable ordination.");
		return false;
	}
	Int64 threv = 2000;
	double den = std::pow(2.0, paramPack->mostdLevel_ - 1);
	Int64 arrayLen = Int64((paramPack->xMax_ - paramPack->xMin_) * Int64(paramPack->yMax_ - paramPack->yMin_) *Int64(paramPack->zMax_ - paramPack->zMin_) / (den*den*den));
	if (Int64((paramPack->xMax_ - paramPack->xMin_) / den) > threv ||
		Int64((paramPack->yMax_ - paramPack->yMin_) / den) > threv || Int64((paramPack->zMax_ - paramPack->zMin_) / den) > threv
		|| arrayLen > 4000000000){
		QMessageBox::warning(this, "warning", "User-defined image range is out of ordinary GPU memory. Please set suitable ordination.");
		return false;
	}
	return true;
}

bool SparseTracer::SaveSeperateTreeAsOne(QString savePath)
{
	if (savePath.isEmpty() || paramPack->activeTree == NULL || !paramPack->separateTree) {
		return false;
	}
	LOG(INFO) << "Begin SaveSeperateTree";
	std::vector<std::string> fileList;
	NG_CREATE_DYNAMIC_CAST(NeuronPopulation, tmpTree, paramPack->separateTree);
	QString swcPath;
	for (size_t i = 0; i < tmpTree->m_pop.size(); ++i) {
		//save tree iteratively
		swcPath = savePath.section('.', 0, -2) + QString(tr("_%1").arg(i, 3, 10, QChar('0'))) + QString(".swc");
		fileList.push_back(swcPath.toStdString());
	}
	auto writer = TreeWriter::New();
	writer->SetInput(paramPack->separateTree);
	writer->SetParam(paramPack);
	writer->SetOutputFileName(fileList);
	auto res = writer->Update();
	if (!res->success()) {
		NG_ERROR_MESSAGE("tree save failed.");
		QMessageBox::warning(this, QString::fromStdString(res->ErrorClassName()), QString::fromStdString(res->ErrorInfo()));
		return false;
	}
	//if (paramPack->dataMode_ == MOSTD) {
	QString stackPath = savePath.section('.', 0, -2) + QString(".sta");
	NGTraceStackWriter stackWriter = TraceStackWriter::New();
	stackWriter->SetSwcFileNames(fileList);
	stackWriter->SetInput(paramPack->separateTree);
	stackWriter->SetParam(paramPack);
	stackWriter->SetOutputFileName(stackPath.toStdString());
	auto stres = stackWriter->Update();
	if (!stres->success()){
		QMessageBox::warning(this, QString::fromStdString(stres->ErrorClassName()), QString::fromStdString(stres->ErrorInfo()));
		return false;
	}
	printf("save tree as %s. bound pt and bound dir: %llu\n", swcPath.toStdString().c_str(), paramPack->activeTree->m_traceInitInfo.size());
	//LOG(INFO) << "save " << tmpTree->m_pop.size() << "trees ";
	//}
	//else printf("save tree as %s.\n", swcPath.toStdString().c_str());
	printf("save tree as %s.\n", swcPath.toStdString().c_str());
	return true;
}

void SparseTracer::on_actionSaveImage_triggered()
{
	if (paramPack->OrigImage) {
		QString saveName = QFileDialog::getSaveFileName(this, "save current image", paramPack->defaultDir, "TIF (*.tif)");
		paramPack->defaultDir = saveName.section('/', 0, -2);
		NGImageWriter writer = ImageWriter::New();
		writer->SetOutputFileName(saveName.toStdString());
		writer->SetInput(paramPack->OrigImage);//OrigImage
		auto res = writer->Update();
		if (!res->success()){
			QMessageBox::warning(this, QString::fromStdString(res->ErrorClassName()), QString::fromStdString(res->ErrorInfo()));
			return;
		}
		printf("save image:%s\n", (const char*)saveName.toLocal8Bit());
	}
}

void SparseTracer::on_actionSaveBack_triggered()
{
	if (paramPack->BackImage) {
		QString saveName = QFileDialog::getSaveFileName(this, "save current image", paramPack->defaultDir, "TIF (*.tif)");
		paramPack->defaultDir = saveName.section('/', 0, -2);
		NGImageWriter writer = ImageWriter::New();
		writer->SetOutputFileName(saveName.toStdString());
		writer->SetInput(paramPack->BackImage);
		auto res = writer->Update();
		if (!res->success()) {
			QMessageBox::warning(this, QString::fromStdString(res->ErrorClassName()), QString::fromStdString(res->ErrorInfo()));
			return;
		}
		printf("save bakc image:%s\n", (const char*)saveName.toLocal8Bit());
	}
}

void SparseTracer::on_actionSaveSVM_triggered()
{
	if (paramPack->w.rows() == 0) {
		printf("w is not prepared. \n");
		return;
	}
	QString savePath = QFileDialog::getSaveFileName(this, "save SVM",
		paramPack->defaultDir, tr("txt (*.txt)"));
	if (savePath.isEmpty()) {
		printf("exit save svm.\n");
		return;
	}
	paramPack->defaultDir = savePath.section('/', 0, -2);
	FILE* fp = fopen(savePath.toStdString().c_str(), "w");
	if (!fp) {
		NG_ERROR_MESSAGE("save svm file error.");
		return;
	}
	for (int i = 0; i < paramPack->w.rows(); ++i) fprintf(fp, "%lf ", paramPack->w(i));
	fprintf(fp, "\n %lf", paramPack->b);
	fclose(fp);
	printf("save svm as %s\n", savePath.toStdString().c_str());
}

void SparseTracer::on_actionSaveTree_triggered()
{
	if (!paramPack->separateTree){
		printf("no separate tree.\n");
		return;
	}
	NG_CREATE_DYNAMIC_CAST(NeuronPopulation, tmpTree, paramPack->separateTree);
	if (paramPack->activeTree == NULL) {
		QMessageBox::warning(this, "warning", "there is no active tree to save.");
		return;
	}
	if (tmpTree->m_pop.empty()) {
		QMessageBox::warning(this, "warning", "there is no tree to save.");
		return;
	}

	//if (!CheckReadImageRangeValid()) return;

	QString savePath = QFileDialog::getSaveFileName(this, "save Tree", paramPack->defaultDir, tr("swc (*.swc)"));
	if (savePath.isEmpty()) {
		return;
	}
	paramPack->defaultDir = savePath.section('/', 0, -2);
	SaveSeperateTreeAsOne(savePath);
	timer->start(timerForSaveTime_);
}

void SparseTracer::on_actionClear_triggered()
{
	if (QMessageBox::No == QQUESTION("Clear All", "All unsaved data will be cleared, continue?")) {
		//QMessageBox::question(this, "Clear All", "All unsaved data will be cleared, continue?", QMessageBox::Yes | QMessageBox::No, QMessageBox::No)
		QMessageBox::information(this, "Information", "Please restart the trace mission.");
		return;
	}
	ClearAll();
	ui->action2DView->setEnabled(false);
	ui->actionBackTrack->setEnabled(false);
	ui->actionChoose_Line->setEnabled(false);
	ui->actionChoose_Vertex->setEnabled(false);
	ui->actionCut_Vertex->setEnabled(false);
	ui->actionDelete_Line->setEnabled(false);
	ui->actionDelete_Soma->setEnabled(false);
	ui->actionDelete_Tree->setEnabled(false);
	ui->actionDraw_Line->setEnabled(false);
	ui->actionNGTree_Trace->setEnabled(false);
	ui->actionPick_Soma->setEnabled(false);
	ui->actionRun->setEnabled(false);
	ui->actionSelect_Tree->setEnabled(false);
	ui->actionStop->setEnabled(false);
	ui->actionTrain->setEnabled(false);
	ui->actionZoom_in->setEnabled(false);
	ui->actionZoom_out->setEnabled(false);
	LOG(INFO) << "Clear all data.";
}

void SparseTracer::UpdateStatus_Slot()
{
	switch (statusTip_){
	case 0:
		ui->statusBar->showMessage(statusString_ + ".");
		++statusTip_;
		break;
	case 1:
		ui->statusBar->showMessage(statusString_ + "...");
		++statusTip_;
		break;
	case 2:
		ui->statusBar->showMessage(statusString_ + "......");
		statusTip_ = 0;
		break;
	default:
		break;
	}
}

void SparseTracer::on_actionRun_triggered()
{
	if (!paramPack->activeTree) {
		QMessageBox::warning(this, "no active tree", "please set active tree.");
		return;
	}
	ui->actionTreeChecker->setEnabled(true);
	paramPack->isLocalTrace_ = false;
	//initial
	BOUNDINFOLIST tmpCurBoundInfo;
	if (!largeFilter) {
		largeFilter = LargeSparseTraceFilter::New();
		largeFilter->SetParam(paramPack);
		if (!paramPack->separateTree) paramPack->separateTree = std::make_shared<NeuronPopulation>();
		if (paramPack->separateTree) largeFilter->SetInputOldTree(paramPack->separateTree);
		largeFilter->SetInputMOSTD(mostdReader);
		ui->actionStop->setEnabled(true);
		ui->actionBackTrack->setEnabled(true);
	}
	 {
		 tmpCurBoundInfo = glbox->CurBoundInfo();
		 if (tmpCurBoundInfo.empty()){
			 if (QMessageBox::No == QQUESTION("End current trace", "There is no new trace initial point, or you force to set curve to trace manually.\nwill you end and continue trace another branch?")) {
				 //QMessageBox::question(this, "End current trace", 
				 //"There is no new trace initial point, or you force to set curve to trace manually.\nwill you end and continue trace another branch?",
				 //    QMessageBox::Yes | QMessageBox::No, QMessageBox::No)
				 QMessageBox::information(this, "Information", "Please draw the trace line to add trace initial points.");
				 return;
			 }
		 }
	 }
	 if (!timer->isActive()) timer->start(timerForSaveTime_);
	 //record
	 LOG(INFO) << "sum steps is " << sumSteps_ << ", parameter is " << " treeBinThrev: " << paramPack->binThreshold_ << " axonBinThrev: "
		 << paramPack->axonBinaryThreshold_ << " diffuse value: " << paramPack->diffuseValue_
		 << " axon diffuse value: " << paramPack->axonDiffuseValue_ << " trace value: " << paramPack->traceValue_
		 << " axon trace value: " << paramPack->axonDiffuseValue_
		 << " max bound num: " << paramPack->maxBoundNum_;
	 ++sumSteps_;
	 std::copy(tmpCurBoundInfo.begin(), tmpCurBoundInfo.end(), std::back_inserter(paramPack->activeTree->m_traceInitInfo));
	 /*start running*/
	 glbox->ToggledForbidden(true);
	 ui->actionStop->setEnabled(false);
	 ui->actionBackTrack->setEnabled(false);
	 ui->actionRun->setEnabled(false);
	 ui->statusBar->showMessage("running...");
	 //worker = new BinaryThread(this);
	 //worker->setParam(largeFilter);
	 //connect(worker, SIGNAL(finished()), this, SLOT(EndRunThread()));
	 //connect(worker, SIGNAL(finished()), worker, SLOT(deleteLater()));
	 //worker->start();
	 EndRunThread();
}

void SparseTracer::EndRunThread()
{
	glbox->ToggledForbidden(false);
	ui->actionStop->setEnabled(true);
	ui->actionBackTrack->setEnabled(true);
	ui->actionRun->setEnabled(true);
	auto lres = largeFilter->Update();
	bool retVal = lres->success();//worker->RetVal();
	if (!retVal){
		ui->actionStop->setEnabled(false);
		ui->actionBackTrack->setEnabled(false);
		glbox->RemoveInitPtActor();
		glbox->RebuildActiveTreeActor();
		glbox->ResetCurForbiddenBoundInfo();
		ui->statusBar->showMessage("Error.", 10000);
		QMessageBox::warning(this, "Trace failed", "Trace failed. Maybe the signal is too weak or you set wrong value.");
	}
	else{
		ui->statusBar->showMessage("current block trace complete.", 10000);
	}
	//get image
	if (!largeFilter->ReleaseCurrentImage()) {
		QMessageBox::information(this, "No new image", "Cannot read new image data.");
		return;
	}
	paramPack->OrigImage = largeFilter->ReleaseCurrentImage();
	//get offset
	settingDock->UpdateImageROILabel();
	paramPack->separateTree = largeFilter->GetOutputTree();
	//display
	settingDock->SetMaxThickness(paramPack->zMax_ - paramPack->zMin_);
	glbox->SetParam(paramPack);
	//
	auto &layerDisplay_ = glbox->GetLayerDisplay();
	auto &branchDisplay_ = glbox->GetBranchDisplay();
	if (layerDisplay_.empty() || layerDisplay_.front() == 0) {
		if (paramPack->OrigImage) 
			glbox->SetInputVolume(paramPack->OrigImage, paramPack->OrigImage->DataType());
		glbox->RebuildActiveTreeActor();
		//release previous auxiliary image
		if (paramPack->LayerImage) paramPack->LayerImage.reset();
	}
	else{// construct new image, layerDisplay_ increment level
		//find correspond sub tree
		auto &topo = paramPack->activeTree->m_Topo;
		size_t id = paramPack->curTraceTopoID_; //root node is invalid
		auto curIt = std::find(topo.begin() + 1, topo.end(), LineID(id));
		if (curIt.node == NULL) {
			NG_ERROR_MESSAGE("debug!!");
			system("pause");
		}
		int curDepth = topo.depth(curIt);
		auto subTree = topo.subtree(curIt, topo.next_sibling(curIt));
		int maxDepth = subTree.max_depth();
		if (maxDepth > 0) {
			layerDisplay_.resize(maxDepth);
			int depthID = curDepth - 1;
			std::for_each(layerDisplay_.begin(), layerDisplay_.end(), [&](int &arg){arg = ++depthID; });
			int branch = 0;
			for (auto iter = topo.begin() + 1; iter != topo.end(); ++iter) {
				if (topo.depth(iter) == curDepth) {
					if (iter->id != id)
						++branch;
					else if (iter->id == id)
						break;
				}
			}
			branchDisplay_.clear(); branchDisplay_.push_back(branch);
			DendriteLevelImageMaker maker;
			maker.SetParam(paramPack);
			maker.SetLayer(layerDisplay_); maker.SetBranch(branchDisplay_);
			maker.SetLayerRadius(layerRadiusDisplay_);
			maker.SetHasSoma(hasSoma_);
			auto res = maker.Update();
			if (!res->success()) {
				NG_ERROR_MESSAGE("Layer Image Construct failed.");
				QMessageBox::warning(this, QString::fromStdString(res->ErrorClassName()), QString::fromStdString(res->ErrorInfo()));
				return;
			}
			glbox->SetInputVolume(paramPack->LayerImage);
			glbox->RebuildActiveTreeActor();
		}
		else{
			QMessageBox::warning(this, "Selective Mode", "Selective parameters are not Updated. Nothing new need to update in 3D viewer.");
		}
	}
	//glbox->RebuildActiveTreeActor();
	glbox->ResetCurForbiddenBoundInfo();
	glbox->UpdateBoxWidget_Slot();
	glbox->UpdateVolumeOpac_Slot();
	ui->action2DView->setChecked(false);
}

void SparseTracer::on_actionStop_triggered()
{
	if (QMessageBox::No == QQUESTION("End current trace", "Would you stop the trace work?")){
		//QMessageBox::question(this, "End current trace",
		//"Would you stop the trace work?",
		//QMessageBox::Yes | QMessageBox::No, QMessageBox::No)){
		return;
	}
	//
	LOG(INFO) << " largeSparseTracer Stop. sum Steps is " << sumSteps_;
	sumSteps_ = 0;
	//
	largeFilter->Stop();
	largeFilter.reset();
	glbox->DeactiveAllTrees();
	ui->actionStop->setEnabled(false);
	ui->actionBackTrack->setEnabled(false);
}

void SparseTracer::on_actionChoose_Line_toggled(bool)
{

}

void SparseTracer::on_actionChoose_Vertex_toggled(bool)
{
}

void SparseTracer::on_actionPick_Soma_toggled(bool)
{
}

void SparseTracer::on_actionSelect_Tree_toggled(bool)
{
}

void SparseTracer::SetROIByBoxWidget_Slot()
{//do not change the paramPack, just display
	glbox->GetBoxWidgetBound(boxWidgetBound);
	boxWidgetBound[0] /= paramPack->xRes_;
	boxWidgetBound[1] /= paramPack->xRes_;
	boxWidgetBound[2] /= paramPack->yRes_;
	boxWidgetBound[3] /= paramPack->yRes_;
	boxWidgetBound[4] /= paramPack->zRes_;
	boxWidgetBound[5] /= paramPack->zRes_;
	/*paramPack->xMin_ = boxWidgetBound[0] / paramPack->xRes_;
	paramPack->xMax_ = boxWidgetBound[1] / paramPack->xRes_;
	paramPack->yMin_ = boxWidgetBound[2] / paramPack->yRes_;
	paramPack->yMax_ = boxWidgetBound[3] / paramPack->yRes_;
	paramPack->zMin_ = boxWidgetBound[4] / paramPack->zRes_;
	paramPack->zMax_ = boxWidgetBound[5] / paramPack->zRes_;*/
	settingDock->UpdateImageROILabel(boxWidgetBound);
}

void SparseTracer::on_actionDelete_Line_triggered()
{
}

void SparseTracer::on_actionCut_Vertex_triggered()
{
	//glbox->DivideLine();//sule
}

void SparseTracer::on_actionDraw_Line_toggled(bool)
{
}

void SparseTracer::on_action2DView_toggled(bool arg)
{
	glbox->p2DView->setChecked(arg);
}

void SparseTracer::on_actionZoom_in_triggered()
{
	if (glbox->Is2DView())
		glbox->TwoDViewZoom(true);
}

void SparseTracer::on_actionZoom_out_triggered()
{
	if (glbox->Is2DView())
		glbox->TwoDViewZoom(false);
}

void SparseTracer::on_actionVisible_toggled(bool arg)
{
	//ui->actionChoose_Line->setChecked(false);
	//ui->actionChoose_Vertex->setChecked(false);
	glbox->ToggleTreesVisible(arg);
}

void SparseTracer::ChangeMoveCursor_Slot(bool arg)
{
	if (arg) {
		QApplication::setOverrideCursor(Qt::SizeAllCursor);
	}
	else {
		QApplication::restoreOverrideCursor();
	}
}

void SparseTracer::TimingSave_Slot()
{
	//return;
	if (!paramPack->separateTree){
		printf("TimingSave_Slot : no separate tree or not running.\n");
		return;
	}
	if (!CheckReadImageRangeValid()) {
		printf("image range is out of true range.\n");
		return;
	}
	QString savePath(paramPack->defaultDir);
	QString fileName = imageFileName_.section('/', -1).section('.', 0, -2);
	savePath += '/' + fileName + (QDateTime::currentDateTime()).toString("yyyy-MM-dd_hh_mm") + ".swc";
	if (!SaveSeperateTreeAsOne(savePath)){
		printf("auto save failed.\n");
		return;
	}
	//delete previous data
	/*if (previousSaveTime_.size() > 2) {
		QDir dir(paramPack->defaultDir);
		QStringList dirFileList = dir.entryList(QStringList() << "*.swc");
		for (int k = 0; k < dirFileList.size(); ++k) {
		QString prefix = paramPack->defaultDir + '/' + dirFileList[k].section('.', 0, -2).section('_', 0, -2);
		bool delFlag = false;
		for (int kk = 0; kk < previousSaveTime_.size() - 2; ++kk) {
		if (prefix.compare(previousSaveTime_[kk]) == 0) {
		delFlag = true;
		break;
		}
		}
		if (!delFlag) continue;
		QFile::remove(paramPack->defaultDir + '/' + dirFileList[k]);
		QString staFileName = paramPack->defaultDir + '/' + dirFileList[k].section('.', 0, -2) + ".sta";
		if (QFile::exists(staFileName)) {
		QFile::remove(staFileName);
		}
		}
		previousSaveTime_.erase(previousSaveTime_.begin(), previousSaveTime_.end() - 2);
		}
		previousSaveTime_.push_back(savePath.section('.',0,-2));*/
	printf("auto save as %s.\n", savePath.toStdString().c_str());
}

void SparseTracer::on_actionNGTree_toggled(bool arg)
{
	ui->toolBar_2->setVisible(arg);
	glbox->SetNormalMode();
}

void SparseTracer::on_actionSBWT_toggled(bool arg)
{
	ui->toolBar->setVisible(arg);
	glbox->SetNormalMode();
}

void SparseTracer::on_actionNGTree_Trace_triggered()
{
	//check whether all soma are in image range.
	{
		NG_CREATE_DYNAMIC_CAST(Soma, tmpSoma, paramPack->SomaList);
		if (!tmpSoma) paramPack->SomaList = std::make_shared<Soma>(); {
			for (size_t i = 0; i < tmpSoma->size(); ++i) {
				auto &cell = tmpSoma->GetCell(i);
				if (cell.x < paramPack->xMin_ || cell.x > paramPack->xMax_ ||
					cell.y < paramPack->yMin_ || cell.y > paramPack->yMax_ ||
					cell.z < paramPack->zMin_ || cell.z > paramPack->zMax_) {
					QMessageBox::warning(this, "Invalid Soma", "There are soma out of current image range, please delete them.");
					return;
				}
			}
		}
		//FILE* fp = fopen("F:/quan/TDatcrop.swc", "r");
		//char line[256];
		//double x, y, z, r; int id, type, pid;
		//while (0 != fgets(line, 256, fp)){
		//    sscanf_s(line, "%d %d %lf %lf %lf %lf %d", &id, &type, &x, &y, &z, &r, &pid);
		//    //tmpSoma->push_back(Soma::MakeCell(x / paramPack->xRes_, y / paramPack->yRes_, z));
		//    tmpSoma->push_back(Soma::MakeCell(x , y , z));
		//}
		//fclose(fp);
	}
	//
	ui->statusBar->showMessage("NeuroGPS-Tree tracing.");
	printf("NeuroGPS-Tree.\n");
	NGNeuroGPSTreeFilter ngt = NeuroGPSTreeFilter::New();
	ngt->SetParam(paramPack);
	ngt->SetSoma(paramPack->SomaList);
	auto traceres = ngt->Update();
	if (!traceres->success()) {
		ui->statusBar->showMessage("NeuroGPS-Tree trace error.", 10000);
		QMessageBox::warning(this, QString::fromStdString(traceres->ErrorClassName()), QString::fromStdString(traceres->ErrorInfo()));
		return;
	}
	IDataPointer tmpNewSeparateTree = ngt->ReleaseData();
	NG_CREATE_DYNAMIC_CAST(NeuronPopulation, tmpTree, paramPack->separateTree);
	NG_CREATE_DYNAMIC_CAST(NeuronPopulation, tmpNewTree, tmpNewSeparateTree);
	tmpTree->m_pop.clear();
	for (size_t i = 0; i < tmpNewTree->m_pop.size(); ++i) {
		tmpTree->m_pop.push_back(tmpNewTree->m_pop[i]);
	}
	std::cout << tmpTree->m_pop.size() << std::endl;
	glbox->SetInputNeuronPopulation(paramPack->separateTree);
	ui->statusBar->showMessage("NeuroGPS-Tree trace complete.", 10000);
	paramPack->SomaList = std::make_shared<Soma>();//delete soma
	ui->actionTreeChecker->setEnabled(true);
	settingDock->UpdateTreeWidget();
}

void SparseTracer::on_actionDelete_Soma_triggered()
{
}

void SparseTracer::on_actionDelete_Tree_triggered()
{
}

void SparseTracer::on_actionTest_triggered()
{
	//printf("test");
	//imageFileName_ = "X:/TDI11107005/TDI11107005_mostd/TDifor211scale.mostd";
	//NGTraceStackReader stackReader;
	//mostdReader = MOSTDReader::New();
	//stackReader = TraceStackReader::New();
	//stackReader->SetInput(seperateTree);
	//mostdReader->SetParam(paramPack);
	//if (!mostdReader->SetInputFileName(imageFileName_.toStdString())) {
	//    printf("mostd file is invalid.\n");
	//    return;
	//}
	//mostdReader->GetMaxRange(xRangeMax_, yRangeMax_, zRangeMax_);
	//settingDock->SetImageRange(xRangeMax_, yRangeMax_, zRangeMax_);
	//QString stackFileName = "Z:/hzhou/Trace_Mission/zhouhang/2017-02-28_15_09_33.sta";// "X:/TDI11107005/TDI11107005_mostd/new.sta";
	//settingDock->SetApplyReadImageButton(true);
	//paramPack->highOpac_ = paramPack->dataType_ == IMAGE8 ? 255 : 1000;
	//settingDock->SetDisplayOpacValue(0, paramPack->highOpac_);
	//glbox->SetParam(paramPack);
	////
	//ui->action2DView->setEnabled(true);
	//ui->actionBackTrack->setEnabled(true);
	//ui->actionChoose_Line->setEnabled(true);
	//ui->actionChoose_Vertex->setEnabled(true);
	//ui->actionNGTree_Trace->setEnabled(true);
	//ui->actionPick_Soma->setEnabled(true);
	//ui->actionRun->setEnabled(true);
	//ui->actionSelect_Tree->setEnabled(true);
	//ui->actionStep->setEnabled(true);
	//ui->actionStop->setEnabled(true);
	//ui->actionTrain->setEnabled(true);
	//ui->actionVisible->setEnabled(true);
	////
	//if (!stackFileName.isEmpty()){
	//    stackReader->SetParam(paramPack);
	//    stackReader->SetInputFileName(stackFileName.toStdString());
	//    if (!stackReader->Update()) return;
	//}
	//else{
	//    printf("Image Read, please set the trace start points.\n");
	//    ui->actionRun->setEnabled(true);
	//    return;
	//}
	//
	//printf("Image and Stack information Reading complete.\n");
	////-------------------display-----------//
	//if (paramPack->xMax_ - paramPack->xMin_ > 19 && paramPack->yMax_ - paramPack->yMin_ > 19 && paramPack->zMax_ - paramPack->zMin_ > 19
	//    && paramPack->xMin_ > -1 && paramPack->yMin_ > -1 && paramPack->zMin_ > -1){
	//    mostdReader->SetParam(paramPack);
	//    mostdReader->UpdateROI();
	//    mostdReader->Update();
	//    printf("Image read complete.\n");
	//    paramPack->OrigImage = mostdReader->ReleaseData();
	//    glbox->UpdateOffset();
	//    glbox->SetInputVolume(paramPack->OrigImage, paramPack->dataType_);
	//    glbox->UpdateBoxWidget_Slot();
	//}

	//glbox->SetInputNeuronPopulation(seperateTree);
	//glbox->SetClippingDistance_Slot(paramPack->thicknesss_);
	//settingDock->UpdateGUI();
	//update();
	//test
	NGImageReader imgReader = ImageReader::New();
	imgReader->SetInputFileName("E:/MyProject/SparseTracer2016/shiwei6.tif");
	imgReader->SetParam(paramPack);
	imgReader->Update();
	paramPack->OrigImage = imgReader->ReleaseData();
	NGBinaryFilter filter = BinaryFilter::New();
	filter->SetInput(paramPack->OrigImage);
	filter->SetParam(paramPack);
	filter->SetThreshold(paramPack->axonBinaryThreshold_);//binaryThreshold_
	filter->Update();
	paramPack->BinImage = filter->ReleaseData();
	paramPack->BackImage = filter->ReleaseBackNoiseImage();
	NGSparseTraceAxonFilter sparseAxonFilter = SparseTraceAxonFilter::New();
	sparseAxonFilter->SetInput(paramPack->OrigImage);
	sparseAxonFilter->SetInputBack(paramPack->BackImage);
	sparseAxonFilter->SetInputBin(paramPack->BinImage);
	NG_SMART_POINTER_DEFINE(Soma, tmpSoma);
	NG_SMART_POINTER_NEW_DEFAULT(Soma, tmpSoma);
	tmpSoma->push_back(Cell(0, 41.068401, 2.430548, 31.001518, 1, 1, 1, 1, 1, 1));
	sparseAxonFilter->SetSoma(tmpSoma);
	sparseAxonFilter->SetParam(paramPack);
	sparseAxonFilter->SetInitialDirection(Vec3d(0.0, 0.6, 0.0));
	sparseAxonFilter->Update();
	IDataPointer tmpRelease = sparseAxonFilter->ReleaseData();
	NG_CREATE_DYNAMIC_CAST(TreeCurve, tmpTree, tmpRelease);
	const std::vector<VectorVec5d>& tmpCurve5d = tmpTree->GetCurve();
	//IDataPointer m_Source;
	NG_SMART_POINTER_NEW_DEFAULT(NeuronPopulation, m_Source);
	NG_CREATE_DYNAMIC_CAST(NeuronPopulation, curTree, m_Source);
	curTree->m_pop.push_back(std::make_shared<NeuronTree>());
	std::vector<VectorVec5d> traceLines;
	Vec3d root = tmpCurve5d.front().front().block(0, 0, 3, 1);
	//size_t traceFirstID = 0;
	tree<LineID> traceTopo; traceTopo.set_head(LineID(std::numeric_limits<size_t>::max()));
	for (size_t k = 0; k < tmpCurve5d.size(); ++k){
		if (tmpCurve5d[k].size() < 2) continue;
		VectorVec5d tmpCurve;
		Vec5d tmp;
		for (size_t l = 0; l < tmpCurve5d[k].size(); ++l){
			tmp = tmpCurve5d[k][l];
			tmp.block(0, 0, 3, 1) = tmp.block(0, 0, 3, 1);
			tmpCurve.push_back(tmp);
		}
		traceLines.push_back(VectorVec5d());
		traceLines.back().swap(tmpCurve);
		traceTopo.append_child(traceTopo.begin(), LineID(k));
	}
	curTree->m_pop[0]->m_curveList = traceLines;
	curTree->m_pop[0]->m_Topo = traceTopo;
	KPTREE::SetTreeRoot(*(curTree->m_pop[0]), root);
	//curTree->m_poppush_back(tmpCurve5d);
	glbox->SetParam(paramPack);
	//glbox->UpdateOffset();
	glbox->SetInputVolume(paramPack->OrigImage, paramPack->OrigImage->DataType());
	glbox->SetInputNeuronPopulation(m_Source);
	glbox->UpdateBoxWidget_Slot();
	glbox->UpdateVolumeOpac_Slot();
	//glbox->BuildCurInitTracePt(Vec3d(90.544342, 5.770138, 80.000000));
	//ui->actionChoose_Line->setEnabled(true);
	//ui->actionChoose_Vertex->setEnabled(true);
	ui->action2DView->setEnabled(true);
}

bool SparseTracer::GetDrawLineIsChecked()
{
	return ui->actionDraw_Line->isChecked();
}

void SparseTracer::EditActionGroup_triggered(QAction * arg)
{
	if (oldEditToolbarTag == newEditToolbarTag) {
		arg->setChecked(false);
		oldEditToolbarTag = newEditToolbarTag = EDITTOOLBARTAG::NONE;
		return;
	}
	oldEditToolbarTag = newEditToolbarTag;
}

void SparseTracer::SetProjectMaxThickness_Slot(int arg)
{
	switch (arg){
	case 0:
		settingDock->SetMaxThickness(paramPack->zMax_ - paramPack->zMin_);
		break;
	case 1:
		settingDock->SetMaxThickness(paramPack->xMax_ - paramPack->xMin_);
		break;
	case 2:
		settingDock->SetMaxThickness(paramPack->yMax_ - paramPack->yMin_);
		break;
	default:
		break;
	}
}

void SparseTracer::MaxBoundNumChanged_Slot(int)
{
	glbox->BuildCurrentBoundaryCollection();
}

void SparseTracer::Set2DViewStatus_Slot(bool arg)
{
	ui->actionZoom_in->setEnabled(arg);
	ui->actionZoom_out->setEnabled(arg);
}

void SparseTracer::SetDrawStatus_Slot(bool arg)
{
	ui->actionRun->setEnabled(!arg);
	ui->actionStop->setEnabled(!arg);
}

void SparseTracer::on_actionTrain_triggered()
{
	if (!largeFilter) {
		printf("trace not begin\n");
		return;
	}
	largeFilter->Train();
}

void SparseTracer::closeEvent(QCloseEvent * event)
{
	if (QMessageBox::No == QQUESTION("Exit GTree", "GTree will exit. Do you save all your data?")){
		//QMessageBox::question(this, "Exit NeuroGPS-Tree2",  "NeuroGPS-Tree2 will exit. Do you save all your data?",
		//QMessageBox::Yes | QMessageBox::No, QMessageBox::No)) {
		event->ignore();
	}
	else event->accept();
}

void SparseTracer::on_actionTreeChecker_triggered()
{
	//get layerDisplay
	auto &layerDisplay_ = glbox->GetLayerDisplay();
	auto &branchDisplay_ = glbox->GetBranchDisplay();
	//build dialog
	auto tmp = new TreeChecker(this);
	tmp->SetParam(paramPack);
	tmp->SetRadius(layerRadiusDisplay_);
	tmp->SetLayerID(layerDisplay_);
	tmp->SetBranchID(branchDisplay_);
	tmp->SetHasSoma(hasSoma_);
	if (tmp->exec() == QDialog::Accepted){
		std::vector<int> tmpLayer;
		std::vector<int> tmpBranch;
		if (!tmp->GetLayerID(tmpLayer)) return;
		if (!tmp->GetBranchID(tmpBranch)) return;
		if (tmpLayer.empty() || tmpBranch.empty()) return;
		std::sort(tmpLayer.begin(), tmpLayer.end());//ensure 0 is the first
		tmpLayer.erase(std::unique(tmpLayer.begin(), tmpLayer.end()), tmpLayer.end());
		std::sort(tmpBranch.begin(), tmpBranch.end()); tmpBranch.erase(std::unique(tmpBranch.begin(), tmpBranch.end()), tmpBranch.end());
		std::for_each(tmpBranch.begin(), tmpBranch.end(), [](int &arg){ --arg; });
		if (tmpLayer[0] < 0 || tmpBranch[0] < 0) return;
		if ((tmpLayer[0] == 0 && layerDisplay_[0] == 0) //both have 0
			|| (NGUtility::CompareVector(layerDisplay_, tmpLayer) && NGUtility::CompareVector(branchDisplay_, tmpBranch)
			&& tmp->GetRadius() == layerRadiusDisplay_)) {//all the same
			printf("You do not change any parameters to display.\n");
			return;
		}
		layerDisplay_.swap(tmpLayer);
		branchDisplay_.swap(tmpBranch);
		layerRadiusDisplay_ = tmp->GetRadius();
		hasSoma_ = tmp->HasSoma();
	}
	else return;
	delete tmp;
	//process
	if (layerDisplay_.empty() || layerDisplay_.front() == 0) {
		if (paramPack->OrigImage) glbox->SetInputVolume(paramPack->OrigImage);
		glbox->ToggleActiveTreeLayerDisplay(false);
		//glbox->RebuildActiveTreeActor();
	}
	else{// construct new image
		DendriteLevelImageMaker maker;
		maker.SetParam(paramPack);
		maker.SetLayer(layerDisplay_); maker.SetBranch(branchDisplay_);
		maker.SetLayerRadius(layerRadiusDisplay_);
		maker.SetHasSoma(hasSoma_);
		if (!maker.Update()) {
			NG_ERROR_MESSAGE("Layer Image Construct failed.");
			return;
		}
		glbox->SetInputVolume(paramPack->LayerImage, paramPack->LayerImage->DataType());
		glbox->ToggleActiveTreeLayerDisplay(true);
		//glbox->RebuildActiveTreeActor();
	}
}

void SparseTracer::UpdateGLBoxWidget_Slot()
{
	settingDock->GetROIValues(boxWidgetBound);
	glbox->UpdateBoxWidgetWithoutReset_Slot(boxWidgetBound);
}

void SparseTracer::on_actionLocal_Run_triggered()
{
	if (!paramPack->activeTree) {
		QMessageBox::warning(this, "no active tree", "please set active tree.");
		return;
	}
	//initial
	paramPack->isLocalTrace_ = true;
	int oldTreeDepth = paramPack->activeTree->m_Topo.max_depth();
	BOUNDINFOLIST tmpCurBoundInfo;
	IDataPointer mostImageBackup;
	if (paramPack->dataMode_ == READMODE::MOSTD)
		mostImageBackup = paramPack->OrigImage;
	if (!largeFilter) {
		largeFilter = LargeSparseTraceFilter::New();
		largeFilter->SetParam(paramPack);
		if (!paramPack->separateTree) paramPack->separateTree = std::make_shared<NeuronPopulation>();
		if (paramPack->separateTree) largeFilter->SetInputOldTree(paramPack->separateTree);
		largeFilter->SetInputMOSTD(mostdReader);
	}
	 {
		 tmpCurBoundInfo = glbox->CurBoundInfo();
		 if (tmpCurBoundInfo.empty()){
			 if (QMessageBox::No == QQUESTION("End current trace",
				 "There is no new trace initial point, or you force to set curve to trace manually.\nwill you end and continue trace another branch?")) {
				 QMessageBox::information(this, "Information", "Please draw the trace line to add trace initial points.");
				 return;
			 }
		 }
	 }
	 std::copy(tmpCurBoundInfo.begin(), tmpCurBoundInfo.end(), std::back_inserter(paramPack->activeTree->m_traceInitInfo));
	 /*start running, block*/
	 ui->statusBar->showMessage("running...");
	 //
	 auto lres = largeFilter->Update();
	 bool retVal = lres->success();
	 if (paramPack->dataMode_ == READMODE::MOSTD)
		 paramPack->OrigImage = mostImageBackup;
	 if (!retVal){
		 glbox->RemoveInitPtActor();
		 glbox->RebuildActiveTreeActor();
		 glbox->ResetCurForbiddenBoundInfo();
		 ui->statusBar->showMessage("Error.", 10000);
		 //QMessageBox::warning(this, "Trace failed", "Trace failed. Maybe the signal is too weak or the initial point is near bound.");
		 QMessageBox::warning(this, QString::fromStdString(lres->ErrorClassName()), QString::fromStdString(lres->ErrorInfo()));
		 return;
	 }
	 ui->statusBar->showMessage("current block trace complete.", 10000);
	 paramPack->separateTree = largeFilter->GetOutputTree();
	 glbox->RebuildActiveTreeActor();
	 glbox->ResetCurForbiddenBoundInfo();
	 //process
	 auto &layerDisplay_ = glbox->GetLayerDisplay();
	 auto &branchDisplay_ = glbox->GetBranchDisplay();
	 if (layerDisplay_.empty() || layerDisplay_.front() == 0) {
		 if (paramPack->OrigImage) glbox->SetInputVolume(paramPack->OrigImage);
		 glbox->RebuildActiveTreeActor();
	 }
	 else{// construct new image, layerDisplay_ increment level
		 int newTreeDepth = paramPack->activeTree->m_Topo.max_depth();
		 int newDepth = newTreeDepth - oldTreeDepth;
		 if (newDepth > 0) {
			 int last = layerDisplay_.back();
			 for (int kk = 0; kk < newDepth; ++kk) {
				 layerDisplay_.push_back(++last);
			 }
		 }
		 DendriteLevelImageMaker maker;
		 maker.SetParam(paramPack);
		 maker.SetLayer(layerDisplay_); maker.SetBranch(branchDisplay_);
		 maker.SetLayerRadius(layerRadiusDisplay_);
		 maker.SetHasSoma(hasSoma_);
		 if (!maker.Update()) {
			 NG_ERROR_MESSAGE("Layer Image Construct failed.");
			 return;
		 }
		 glbox->SetInputVolume(paramPack->LayerImage);
		 glbox->RebuildActiveTreeActor();
	 }
	 paramPack->isLocalTrace_ = false;
}

bool SparseTracer::ReadImage()
{
	if (paramPack->dataMode_ == READMODE::MOSTD || paramPack->dataMode_ == READMODE::HDF5) {
		auto res = mostdReader->Update();
		if (!res->success()){
			QMessageBox::warning(this, QString::fromStdString(res->ErrorClassName()), QString::fromStdString(res->ErrorInfo()));
			return false;
		}
		paramPack->OrigImage = mostdReader->ReleaseData();
	}
	else {
		NGImageReader reader = ImageReader::New();
		reader->SetInputFileName(imageFileName_.toStdString());
		reader->SetParam(paramPack);
		auto res = reader->Update();
		if (!res->success()) {
			QMessageBox::warning(this, QString::fromStdString(res->ErrorClassName()), QString::fromStdString(res->ErrorInfo()));
			return false;
		}
		paramPack->OrigImage = reader->ReleaseData();
	}

	return true;
}

void SparseTracer::on_actionSaveLayerImage_triggered()
{
	NG_CREATE_DYNAMIC_CAST(SVolume, tmpImg, paramPack->LayerImage);
	if (tmpImg) {
		QString saveFile = QFileDialog::getSaveFileName(this, "save layer image", paramPack->defaultDir, "tif (*.tif)");
		if (saveFile.isEmpty())
			return;
		NGImageWriter writer = ImageWriter::New();
		writer->SetInput(tmpImg);
		writer->SetOutputFileName(saveFile.toStdString());
		auto res = writer->Update();
		if (!res->success()) {
			QMessageBox::warning(this, QString::fromStdString(res->ErrorClassName()), QString::fromStdString(res->ErrorInfo()));
			return;
		}
	}
}

void SparseTracer::on_actionSaveLayerSwc_triggered()
{
	//get layerDisplay
	auto &layerDisplay_ = glbox->GetLayerDisplay();
	auto &branchDisplay_ = glbox->GetBranchDisplay();
	QString saveFile = QFileDialog::getSaveFileName(this, "save layer swc", paramPack->defaultDir, "SWC (*.swc)");
	if (saveFile.isEmpty())
		return;
	NGTreeWriter writer = TreeWriter::New();
	std::vector<std::string> saveFileList; saveFileList.push_back(saveFile.toStdString());
	writer->SetInput(paramPack->separateTree);
	writer->SetOutputFileName(saveFileList);
	writer->SetParam(paramPack);
	if (!writer->UpdateLayer(layerDisplay_, branchDisplay_))
		NG_ERROR_MESSAGE("cannot save layer swc");
}

void SparseTracer::ClearMOSTDCache_Slot()
{
	if (mostdReader && paramPack->dataMode_ == READMODE::MOSTD) {
		NG_CREATE_DYNAMIC_CAST(MOSTDReader, tmpMost, mostdReader);
		if (tmpMost) {
			tmpMost->ClearCache();
			printf("MOSTD Reader have cleared cache.\n");
		}
		printf("not mostd mode.");
	}
}

void SparseTracer::ActivateTree_Slot(int arg)
{
	if (traverseMode_) return;
	NG_CREATE_DYNAMIC_CAST(NeuronPopulation, tmpPop, paramPack->separateTree);
	if (!tmpPop) {
		printf("there is no tree.\n");
		return;
	}
	//paramPack->activeTree = tmpPop->m_pop[arg];
	glbox->RebuildActiveTreeChangedPopulation(arg);
	settingDock->UpdateTreeWidget();
	settingDock->UpdateCheckWidget();//clear check result
}

void SparseTracer::CompareTree_Slot(int arg)
{
	if (traverseMode_) return;
	NG_CREATE_DYNAMIC_CAST(NeuronPopulation, tmpPop, paramPack->separateTree);
	if (!tmpPop) {
		printf("there is no tree.\n");
		return;
	}
	paramPack->activeTree->BuildHash();
	tmpPop->m_pop[arg]->BuildHash();
	int rad, thev;
	QDialog *dia = new QDialog(this);
	QFormLayout *form = new QFormLayout(dia);
	QString str1("Radius: ");
	QSpinBox *sbox1 = new QSpinBox(dia);
	sbox1->setMinimum(1);
	sbox1->setMaximum(10);
	sbox1->setValue(2);
	form->addRow(str1, sbox1);
	QString str2("Threshold: ");
	QSpinBox *sbox2 = new QSpinBox(dia);
	sbox2->setMinimum(1);
	sbox2->setMaximum(1000);
	sbox2->setValue(5);
	form->addRow(str2, sbox2);
	QDialogButtonBox *bBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, Qt::Horizontal, dia);
	form->addRow(bBox);
	connect(bBox, SIGNAL(accepted()), dia, SLOT(accept()));
	connect(bBox, SIGNAL(rejected()), dia, SLOT(reject()));
	if (dia->exec() == QDialog::Accepted){
		rad = sbox1->value();
		thev = sbox2->value();
		qDebug() << rad << " " << thev;
	}
	else{
		return;
	}
	delete dia;
	//KPTREE::FindSymmericDiff(paramPack->activeTree->m_curveList, tmpPop->m_pop[arg]->m_curveList, paramPack->activeTree->m_hashList, tmpPop->m_pop[arg]->m_hashList, rad, thev, paramPack->activeTree->m_suspiciousPositions_);
	decltype(paramPack->activeTree->m_suspiciousPositions_) tmpSusPos;
	decltype(paramPack->activeTree->m_suspiciousCheckState_) tmpSusState;
	KPTREE::FindSymmericDiff(*(paramPack->activeTree), *(tmpPop->m_pop[arg]), rad, thev, tmpSusPos, tmpSusState);
	paramPack->activeTree->m_suspiciousPositions_.swap(tmpSusPos);
	paramPack->activeTree->m_suspiciousCheckState_.swap(tmpSusState);
	if (!paramPack->activeTree->m_suspiciousPositions_.empty()) {
		settingDock->UpdateCheckWidget();
	}
	//paramPack->activeTree->m_hashList.Clear();
	//tmpPop->m_pop[arg]->m_hashList.Clear();
}

void SparseTracer::ToggleTreeVisible_Slot(int arg)
{
	glbox->ToggleTargetTreeVisible(arg);
}

void SparseTracer::GotoDiff_Slot(int arg)
{
	if (!paramPack->activeTree) return;
	if (paramPack->activeTree->m_suspiciousPositions_.size() <= arg) {
		NG_ERROR_MESSAGE("Please debug");
		return;
	}

	if (paramPack->dataMode_ == READMODE::MOSTD) {
		auto &diffXYZ = paramPack->activeTree->m_suspiciousPositions_[arg];
		paramPack->xMin_ = diffXYZ[0] - paramPack->xExtract_;
		paramPack->yMin_ = diffXYZ[1] - paramPack->yExtract_;
		paramPack->zMin_ = diffXYZ[2] - paramPack->zExtract_;
		paramPack->xMax_ = diffXYZ[0] + paramPack->xExtract_;
		paramPack->yMax_ = diffXYZ[1] + paramPack->yExtract_;
		paramPack->zMax_ = diffXYZ[2] + paramPack->zExtract_;
		settingDock->UpdateGUI();
		ApplyReadNewImage_Slot();
	}
	glbox->LightDiffArrow(arg);
}

void SparseTracer::ToggleTraverse_Slot(bool arg)
{
	if (paramPack->dataMode_ == READMODE::TIFF3D) {
		QMessageBox::information(this, "cannot traverse", "If import tiff, traverse mode is invalid.");
		return;
	}
	if (!mostdTraverseReader) return;
	if (!preStore) {
		preStore = Prestore::New();
		preStore->SetParam(mostdTraverseReader, paramPack);
		connect(preStore.get(), SIGNAL(CacheComplete_Signal()), this, SLOT(CacheComplete_Slot()));
	}
	traverseMode_ = arg;
	if (arg) {
		if (!preStore->Initial()) {
			traverseMode_ = false;
			return;
		}
		//mostdReader->ClearCache();
		if (mostdReader && paramPack->dataMode_ == READMODE::MOSTD) {
			NG_CREATE_DYNAMIC_CAST(MOSTDReader, tmpMost, mostdReader);
			if (tmpMost) {
				tmpMost->ClearCache();
			}
		}
	}
	else{
		//mostdTraverseReader->ClearCache();
		if (!preStore->Finish()){
			traverseMode_ = true;
			return;
		}
	}
	glbox->ToggleTraverseMode(arg);
}

//void SparseTracer::Annotate_Slot()
//{
//
//}

void SparseTracer::BackTraverse_Slot()
{
	if (!traverseMode_) return;
	IDataPointer SaveSvolume;
	if (preStore->GetPreviousImage(SaveSvolume) == false){
		QMessageBox::warning(this, "waiting!", "waiting to read image.\n");
		return;
	}
	//ApplyReadNewImage_Slot();
	settingDock->UpdateImageROILabel();
	if (mostdReader){
		printf("INSE %d %d %d %d %d %d\n", paramPack->xMin_, paramPack->xMax_, paramPack->yMin_, paramPack->yMax_, paramPack->zMin_, paramPack->zMax_);
		paramPack->OrigImage = SaveSvolume;
		printf("Image Reading complete.\n");
		//-------------------display-----------//
		settingDock->SetMaxThickness(paramPack->zMax_ - paramPack->zMin_);
		glbox->SetParam(paramPack);
		glbox->SetInputVolume(paramPack->OrigImage, paramPack->OrigImage->DataType());
		glbox->BuildCurrentBoundaryCollection();
		glbox->UpdateBoxWidget_Slot();
		glbox->UpdateVolumeOpac_Slot();
		glbox->UpdateTraverseFlagActor();
		glbox->update();
		QApplication::restoreOverrideCursor();
		ui->action2DView->setEnabled(true);
	}
	else{
		QMessageBox::warning(this, QString("Mostd Warining"), QString("Please Load Mostd"));
	}
}

void SparseTracer::NextTraverse_Slot()
{
	if (!traverseMode_) return;
	IDataPointer SaveSvolume;
	if (preStore->GetNextImage(SaveSvolume) == false){
		QMessageBox::warning(this, "waiting!", "waiting to read image.\n");
		return;
	}
	settingDock->UpdateImageROILabel();
	if (mostdReader){
		printf("INSE %d %d %d %d %d %d\n", paramPack->xMin_, paramPack->xMax_, paramPack->yMin_, paramPack->yMax_, paramPack->zMin_, paramPack->zMax_);
		paramPack->OrigImage = SaveSvolume;
		printf("Image Reading complete.\n");
		//-------------------display-----------//
		settingDock->SetMaxThickness(paramPack->zMax_ - paramPack->zMin_);
		glbox->SetParam(paramPack);
		glbox->SetInputVolume(paramPack->OrigImage, paramPack->OrigImage->DataType());
		glbox->BuildCurrentBoundaryCollection();
		glbox->UpdateBoxWidget_Slot();
		glbox->UpdateVolumeOpac_Slot();
		glbox->UpdateTraverseFlagActor();
		glbox->update();
		QApplication::restoreOverrideCursor();
		ui->action2DView->setEnabled(true);
	}
	else{
		QMessageBox::warning(this, QString("Mostd warning"), QString("Please Load Mostd"));
	}
}

void SparseTracer::ToggleShowDiffArrow_Slot(bool arg)
{
	glbox->ToggleDiffArrows(arg);
}

void SparseTracer::StartTraverseFromHere_Slot(int lineid, int beg)
{
	if (!traverseMode_) return;
	preStore->SetTraversePosition(lineid, beg);
	printf("traverse start here.");
}

void SparseTracer::ResetTraverse_Slot()
{
	if (!traverseMode_) return;
	ToggleTraverse_Slot(false);
	preStore->Reset();
	ToggleTraverse_Slot(true);
}

void SparseTracer::on_actionStart_triggered()
{
	runningTimer_->start(120000);
	ui->actionStart->setEnabled(false);
	ui->actionHalt->setEnabled(true);
}

void SparseTracer::on_actionHalt_triggered()
{
	runningTimer_->stop();
	ui->actionStart->setEnabled(true);
	ui->actionHalt->setEnabled(false);
}

void SparseTracer::on_actionReset_triggered()
{
	if (QMessageBox::No == QMessageBox::warning(this, "timer reset", "Are you sure to reset running timer ?")) return;
	runningTimer_->stop();
	paramPack->runningMinutes_ = 0;
	ui->actionStart->setEnabled(true);
	ui->actionHalt->setEnabled(false);
}

void SparseTracer::AddRunningTime_Slot()
{
	++paramPack->runningMinutes_;
}

void SparseTracer::CacheComplete_Slot()
{
	ui->statusBar->showMessage("cache read complete.", 3000);
}

void SparseTracer::on_actionNeuroGPS_triggered()
{
	if (!paramPack->OrigImage) {
		QMessageBox::warning(this, "No image data", "Please import image data first.");
		return;
	}
	QDialog *tmpDialog = new QDialog(this);
	QLabel *rLabel = new QLabel("radius", tmpDialog);
	QDoubleSpinBox *radiusSpin = new QDoubleSpinBox(tmpDialog);
	QPushButton *okBtn = new QPushButton(tr("OK"), tmpDialog);
	QPushButton *cancelBtn = new QPushButton(tr("Cancel"), tmpDialog);
	QGridLayout *gridLayout = new QGridLayout(tmpDialog);
	gridLayout->addWidget(rLabel, 0, 0);
	gridLayout->addWidget(radiusSpin, 0, 1);
	gridLayout->addWidget(okBtn, 1, 0);
	gridLayout->addWidget(cancelBtn, 1, 1);
	tmpDialog->setLayout(gridLayout);
	connect(okBtn, SIGNAL(clicked()), tmpDialog, SLOT(accept()));
	connect(cancelBtn, SIGNAL(clicked()), tmpDialog, SLOT(reject()));
	double minRadius = 3.0;
	radiusSpin->setValue(minRadius);
	tmpDialog->show();
	if (tmpDialog->exec() == QDialog::Accepted)
		minRadius = radiusSpin->value();
	else return;
	ui->statusBar->showMessage("Probing...");
	IDataPointer scaleImg;
	NGImageScaleFiter scaler = ImageScaleFiter::New();
	scaler->SetInput(paramPack->OrigImage);
	scaler->SetParam(paramPack);
	auto scres = scaler->Update();
	if (!scres->success()){
		return;
	}
	scaleImg = scaler->ReleaseData();
	NGBinaryFilter binaryFilter = BinaryFilter::New();
	binaryFilter->SetInput(scaleImg);
	//binaryFilter->SetInput(scaleImg_);//param->OrigImage
	binaryFilter->SetParam(paramPack);
	binaryFilter->SetThreshold(paramPack->binThreshold_);
	auto bres = binaryFilter->Update();
	if (!bres->success()) {
		QMessageBox::warning(this, QString::fromStdString(bres->ErrorClassName()), QString::fromStdString(bres->ErrorInfo()));
		return;
	}
	paramPack->BinImage = binaryFilter->ReleaseData();
	paramPack->BackImage = binaryFilter->ReleaseBackNoiseImage();
	NG_SMART_POINTER_DEFINE_NEW_DEFAULT(VectorVec3i, binPtSet);
	binPtSet->swap(binaryFilter->GetBinPtSet());
	NGProbeFilter probeFilter = ProbeFilter::New();
	//probeFilter->SetInput(paramPack->OrigImage);
	probeFilter->SetInput(scaleImg);
	probeFilter->SetInputBin(paramPack->BinImage);
	probeFilter->SetMinRadius(minRadius);
	probeFilter->SetThreadNum(paramPack->threadNum_);
	probeFilter->SetInputBinPtSet(binPtSet);
	auto pres = probeFilter->Update();
	if (!pres->success()) {
		QMessageBox::warning(this, QString::fromStdString(pres->ErrorClassName()), QString::fromStdString(pres->ErrorInfo()));
		return;
	}
	paramPack->SomaList = probeFilter->ReleaseData();
	NG_CREATE_DYNAMIC_CAST(Soma, tmpSomaList, paramPack->SomaList);
	NG_CREATE_DYNAMIC_CAST(SVolume, tmpScale, scaleImg);
	for (auto &it : tmpSomaList->GetAllCell()) {
		it.x *= paramPack->xScale_;
		it.y *= paramPack->yScale_;
		it.x += paramPack->xMin_;
		it.y += paramPack->yMin_;
		it.z += paramPack->zMin_;
		//it.x *= paramPack->xRes_;
		//it.y *= paramPack->yRes_;
		//it.z *= paramPack->zRes_;
	}
	printf("Soma Locate Complete.\n");
	ui->statusBar->showMessage("Soma Locate Complete", 10000);
	paramPack->BinImage.reset(); paramPack->BackImage.reset();
	//clock_t end = clock();
	//printf("Trace use %lf hour\n", double(end - beg) / double(3600000) );
	glbox->UpdateSomaList();
	glbox->update();
}

void SparseTracer::on_actionSaveSoma_triggered()
{
	if (!paramPack->SomaList) {
		QMessageBox::information(this, "No soma data", "There is no soma data, no need to save.");
		return;
	}
	QString somaFile = QFileDialog::getSaveFileName(this, "Save Soma", paramPack->defaultDir, "PBO(*.pbo)");
	if (somaFile.isEmpty()) return;
	FILE *fp = fopen(somaFile.toStdString().c_str(), "w");
	if (!fp) {
		QMessageBox::information(this, "Cannot write soma file", "cannot create soma file to write.");
		return;
	}
	NG_CREATE_DYNAMIC_CAST(Soma, tmpSomaList, paramPack->SomaList);
	int id = 0;
	for (auto &it : tmpSomaList->GetAllCell())
		fprintf(fp, "%d 1 %lf %lf %lf 1.0 -1\n", ++id, it.x * paramPack->xRes_, it.y * paramPack->yRes_, it.z * paramPack->zRes_);
	fclose(fp);
	printf("soma list has been write into %s\n", somaFile.toStdString().c_str());
}

void SparseTracer::on_actionOpenSoma_triggered()
{
	if (paramPack->SomaList) {
		if (QMessageBox::No == QQUESTION("Soma read warning", "If you read new soma list, the existed soma will be removed. Will you continue?"))
			return;
	}
	QString somaFile = QFileDialog::getOpenFileName(this, "Open soma", paramPack->defaultDir, "PBO(*.pbo)");
	FILE * fp = fopen(somaFile.toStdString().c_str(), "r");
	if (!fp) {
		QMessageBox::information(this, "Cannot read soma list", "Please ensure the soma file is existed.");
		return;
	}
	paramPack->SomaList.reset();
	NG_SMART_POINTER_NEW_DEFAULT(Soma, paramPack->SomaList);
	NG_CREATE_DYNAMIC_CAST(Soma, tmpSomaList, paramPack->SomaList);
	Cell tmpCell;
	char line[256];
	int offset = 0, id, type, pid;
	while (!feof(fp)) {
		if (fgets(line, 256, fp) == 0) continue;
		if (line[0] == '#' || line[0] == '\n' || line[0] == '/' || line[0] == '\\') continue;
		offset = 0;
		if (line[0] == ' ') offset = 1;
		sscanf(line + offset, "%d %d %lf %lf %lf %lf %d", &id, &type, &tmpCell.x, &tmpCell.y, &tmpCell.z, &tmpCell.r, &pid);
		tmpCell.x /= paramPack->xRes_;
		tmpCell.y /= paramPack->yRes_;
		tmpCell.z /= paramPack->zRes_;
		tmpSomaList->push_back(tmpCell);
	}
	fclose(fp);
	glbox->UpdateSomaList();
}

int SparseTracer::SetScale(double res)
{
	int scale;
	if (res > 0.25 && res <= 0.5) {
		scale = 2;
	}
	else if (res <= 0.25) {
		scale = 3;
	}
	else{
		scale = 1;
	}
	return scale;
}

void SparseTracer::OpacAdjustApply_Slot()
{
	if (paramPack->highDestOpac_ != 0 && paramPack->highAdjustOpac_ != 0) {
		if (paramPack->OrigImage) {
			NG_CREATE_DYNAMIC_CAST(SVolume, tmpOrig, paramPack->OrigImage);
			if (tmpOrig) {
				tmpOrig->Adjust(paramPack->lowAdjustOpac_, paramPack->highAdjustOpac_, paramPack->lowDestOpac_, paramPack->highDestOpac_);
				glbox->SetInputVolume(tmpOrig, tmpOrig->DataType());
				glbox->update();
			}
		}
	}
}

void SparseTracer::PreviewOpacAdjust_Slot(bool arg)
{
	if (arg) {
		if (paramPack->OrigImage) {
			if (paramPack->previewImage) paramPack->previewImage.reset();
			if (!paramPack->previewImage) NG_SMART_POINTER_NEW_DEFAULT(SVolume, paramPack->previewImage);
			NG_CREATE_DYNAMIC_CAST(SVolume, tmpPreview, paramPack->previewImage);
			NG_CREATE_DYNAMIC_CAST(SVolume, tmpOrig, paramPack->OrigImage);
			if (tmpPreview && tmpOrig &&paramPack->highAdjustOpac_ > 0 && paramPack->highDestOpac_ > 0) {
				tmpPreview->QuickCopy(*tmpOrig);
				tmpPreview->SetResolution(tmpOrig->XResolution(), tmpOrig->YResolution(), tmpOrig->ZResolution());
				tmpPreview->Adjust(paramPack->lowAdjustOpac_, paramPack->highAdjustOpac_, paramPack->lowDestOpac_, paramPack->highDestOpac_);
				glbox->SetInputVolume(tmpPreview, tmpPreview->DataType());
				ui->statusBar->showMessage("Preview opacity adjust.", 10000);
			}
			else if (!tmpPreview){
				QMessageBox::warning(this, "There are other review ", "Please deactivate other preview.");
			}
			else if (!tmpOrig)
				QMessageBox::warning(this, "The type of OrigImage is not  right", "Please Debug.");
			else
				QMessageBox::warning(this, "Wrong Opac Value ", "Please input correct opacity value.");
		}
		else
			QMessageBox::warning(this, "No original image", "Please import image data first.");
	}
	else{
		if (paramPack->previewImage) {
			paramPack->previewImage.reset();
			if (paramPack->OrigImage) glbox->SetInputVolume(paramPack->OrigImage, paramPack->OrigImage->DataType());
			else{
				QMessageBox::warning(this, "No original image", "You may delete original image data? 3D viewer will not update.");
			}
		}
	}
	glbox->update();
}

void SparseTracer::PreviewBinary_Slot(bool arg)
{
	if (arg) {
		if (paramPack->OrigImage) {
			if (paramPack->previewImage) paramPack->previewImage.reset();
			//NG_CREATE_DYNAMIC_CAST(SVolume, tmpPreview, paramPack->previewImage);
			NG_CREATE_DYNAMIC_CAST(SVolume, tmpOrig, paramPack->OrigImage);
			if (tmpOrig) {
				IDataPointer scaleImg;
				NGImageScaleFiter scaler = ImageScaleFiter::New();
				scaler->SetInput(paramPack->OrigImage);
				scaler->SetParam(paramPack);
				scaler->Update();
				scaleImg = scaler->ReleaseData();
				NGBinaryFilter binfilter = BinaryFilter::New();
				binfilter->SetInput(scaleImg);
				binfilter->SetParam(paramPack);
				binfilter->SetThreshold(paramPack->binThreshold_);
				auto res = binfilter->Update();
				if (!res->success()) {
					QMessageBox::warning(this, QString::fromStdString(res->ErrorClassName()), QString::fromStdString(res->ErrorInfo()));
					return;
				}
				paramPack->previewImage = binfilter->ReleaseData();
				glbox->SetInputVolume(paramPack->previewImage, paramPack->previewImage->DataType());
				ui->statusBar->showMessage("Preview dendrite binarization.", 10000);
			}
		}
		else{
			QMessageBox::warning(this, "No original image", "Please import image data first.");
		}
	}
	else{
		if (paramPack->previewImage) {
			paramPack->previewImage.reset();
			if (paramPack->OrigImage) glbox->SetInputVolume(paramPack->OrigImage, paramPack->OrigImage->DataType());
			else{
				QMessageBox::warning(this, "No original image", "You may delete original image data? 3D viewer will not update.");
			}
		}
	}
	glbox->update();
}

void SparseTracer::PreviewAxonBinary_Slot(bool arg)
{
	if (arg) {
		if (paramPack->OrigImage) {
			if (paramPack->previewImage) paramPack->previewImage.reset();
			//NG_CREATE_DYNAMIC_CAST(SVolume, tmpPreview, paramPack->previewImage);
			NG_CREATE_DYNAMIC_CAST(SVolume, tmpOrig, paramPack->OrigImage);
			if (tmpOrig) {
				IDataPointer scaleImg;
				NGImageScaleFiter scaler = ImageScaleFiter::New();
				scaler->SetInput(paramPack->OrigImage);
				scaler->SetParam(paramPack);
				scaler->Update();
				scaleImg = scaler->ReleaseData();
				NGBinaryFilter binfilter = BinaryFilter::New();
				binfilter->SetInput(scaleImg);
				binfilter->SetParam(paramPack);
				binfilter->SetThreshold(paramPack->axonBinaryThreshold_);
				auto res = binfilter->Update();
				if (!res->success()) {
					QMessageBox::warning(this, QString::fromStdString(res->ErrorClassName()), QString::fromStdString(res->ErrorInfo()));
					return;
				}
				paramPack->previewImage = binfilter->ReleaseData();
				glbox->SetInputVolume(paramPack->previewImage, paramPack->previewImage->DataType());
				ui->statusBar->showMessage("Preview axon binary result.", 10000);
			}
		}
		else{
			QMessageBox::warning(this, "No original image", "Please import image data first.");
		}
	}
	else{
		if (paramPack->previewImage) {
			paramPack->previewImage.reset();
			if (paramPack->OrigImage) glbox->SetInputVolume(paramPack->OrigImage, paramPack->OrigImage->DataType());
			else{
				QMessageBox::warning(this, "No original image", "You may delete original image data? 3D viewer will not update.");
			}
		}
	}
	glbox->update();
}

void SparseTracer::on_actionSnapshot_triggered()
{
	QString path = QFileDialog::getSaveFileName(this, "Snapshot", paramPack->defaultDir, "PNG (*.png)");
	snapDialog = new QDialog(this);
	QFormLayout *form = new QFormLayout(snapDialog);
	QString pathLabel("File path");
	QLineEdit *pathEdit = new QLineEdit(snapDialog);
	pathEdit->setText(path);
	form->addRow(pathLabel, pathEdit);
	QString str1("Scale: ");
	QSpinBox *sbox1 = new QSpinBox(snapDialog);
	sbox1->setMinimum(1);
	sbox1->setMaximum(10000);
	sbox1->setValue(2);
	form->addRow(str1, sbox1);
	QDialogButtonBox *bBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, Qt::Horizontal, snapDialog);
	form->addRow(bBox);
	connect(bBox, SIGNAL(accepted()), snapDialog, SLOT(accept()));
	connect(bBox, SIGNAL(rejected()), snapDialog, SLOT(reject()));
	int vx;
	if (snapDialog->exec() == QDialog::Accepted){
		vx = sbox1->value();
		path = pathEdit->text();
	}
	else{
		return;
	}
	delete snapDialog;
	snapDialog = NULL;
	QSize oldSz;
	glbox->Snapshot(path, vx);
	glbox->update();
	//glbox->RestoreSnapShot(oldSz);
}

void SparseTracer::ClearAll()
{
	mostdReader.reset();
	mostdTraverseReader.reset();
	glbox->GetLayerDisplay().clear();
	glbox->GetLayerDisplay().push_back(0);
	paramPack->OrigImage.reset();
	paramPack->BinImage.reset();
	paramPack->BackImage.reset();
	paramPack->SomaList.reset();
	paramPack->activeTree.reset();
	binPtSet.clear();
	soma.reset();
	paramPack->separateTree.reset();
	paramPack->separateTree = std::make_shared<NeuronPopulation>();
	imageFileName_.clear();
	glbox->RemoveAllActor();
	glbox->update();
	settingDock->UpdateTreeWidget();
	settingDock->UpdateCheckWidget();
	printf("Clear all data.\n");
}

void SparseTracer::on_actionCreate_HDF5_triggered()
{
	/***************************************
	QString saveH5 = QFileDialog::getSaveFileName(this, "Create HDF5", paramPack->defaultDir, tr("HDF5 (*.h5)"));
	if (saveH5.isEmpty())
	return;
	QStringList qtiffList = QFileDialog::getOpenFileNames(this, "Choose TIFF Sequence", paramPack->defaultDir, tr("TIFF (*.tif)"));
	if (qtiffList.isEmpty())
	return;
	std::vector<std::string> tiffList;
	std::for_each(qtiffList.begin(), qtiffList.end(), [&](const QString &file){tiffList.push_back(file.toStdString()); });
	//input resolution
	QStringList resList;
	while (true){
	QString resStr = QInputDialog::getText(this, "Input resolution", "input resolution text",QLineEdit::Normal, "0.3 0.3 1.0");
	resList = resStr.split(' ');
	if (resList.size() != 3 ){
	if (QMessageBox::No == QQUESTION("Wrong value", "you input wrong value, retry?"))
	return;
	else
	continue;
	}
	//
	#define RESQUESTION if (!ok) if(QMessageBox::No == QQUESTION("Wrong value", "you input wrong value, retry?")) \
	return;  else continue;
	bool ok = true;
	double xres = resList[0].toDouble(&ok);
	RESQUESTION
	double yres = resList[1].toDouble(&ok);
	RESQUESTION
	double zres = resList[2].toDouble(&ok);
	RESQUESTION
	paramPack->xRes_ = xres;
	paramPack->yRes_ = yres;
	paramPack->zRes_ = zres;
	break;
	}**************************************/

	creathdf5dialog creHdf5diag;//////////////////////////////////////////////////////////////
    creHdf5diag.SetParam(paramPack);
    creHdf5diag.exec();
	paramPack->xRes_ = creHdf5diag.xResolution;
	paramPack->yRes_ = creHdf5diag.yResolution;
	paramPack->zRes_ = creHdf5diag.zResolution;

    if (!creHdf5diag.closeFlag && creHdf5diag.tiffFlag && !creHdf5diag.tiffList.empty())
    {
        //create
        writer = HDF5Writer::New();
        writer->SetOutputFileName(creHdf5diag.saveH5.toStdString());
        writer->SetParam(paramPack);
        writer->SetTIFFFileList(creHdf5diag.tiffList);
        QApplication::setOverrideCursor(Qt::WaitCursor);
        ui->statusBar->showMessage("Creating HDF5...");
        progressDialog = new QProgressDialog(this);
        progressDialog->setWindowTitle("Creating HDF5");
        progressDialog->setLabelText("Creating...");
        progressDialog->setRange(0, int(creHdf5diag.tiffList.size()));
        progressDialog->setModal(true);
        progressDialog->setCancelButtonText("Cancel");
        //QApplication::processEvents();
        progressDialog->show();
        //auto res = writer->Update();
        //QMessageBox::warning(this, QString::fromStdString(res->ErrorClassName()), QString::fromStdString(res->ErrorInfo()));
        worker = new BinaryThread(this);
        worker->setParam(writer);
        connect(progressDialog, SIGNAL(canceled()), this, SLOT(StopCreateHDF5_Slot()));
        connect(&(*writer), SIGNAL(Progress_Signal(int)), this, SLOT(UpdateHDF5ProgressBar_Slot(int)));
        connect(worker, SIGNAL(finished()), this, SLOT(EndCreateHDF5_Slot()));
        connect(worker, SIGNAL(finished()), worker, SLOT(deleteLater()));
        worker->start();
    }
}

void SparseTracer::EndCreateHDF5_Slot()
{
    QApplication::restoreOverrideCursor();
    progressDialog->hide();
    delete progressDialog;
    if (!worker->RetVal()->success()) {
        ui->statusBar->showMessage("HDF5 Error.", 5000);
    }
    else
        ui->statusBar->showMessage("HDF5 Complete.", 5000);
}

void SparseTracer::UpdateHDF5ProgressBar_Slot(int arg)
{
    progressDialog->setValue(arg);
}

void SparseTracer::StopCreateHDF5_Slot()
{
    writer->stopFlag_ = true;
}

void SparseTracer::on_actionSaveCaliber_triggered()
{
	if (!paramPack->caliberBulidData.get())
	{
		printf("ERROR! CaliberBuildaer data is null!\n");
		return;
	}
	//get path
	QString pathStr = QFileDialog::getSaveFileName(this, "Save Caliber Data", paramPack->defaultDir, "cbd(*.cbd)");

	//write data
	ofstream Ofile1(pathStr.toStdString());
	if (Ofile1.is_open()){
		int id = 0, indexNum = 0;
		for (auto &it : *(paramPack->caliberBulidData))
		{
			size_t sz = it.size();
			//customer format
			//ID,class,x,y,z,index,father
			for (size_t i = 0; i < sz; ++i){
				Ofile1 << id << " " << 1 << " " << it[i][0] * paramPack->xRes_ << " " << it[i][1] * paramPack->yRes_ << " " << it[i][2] * paramPack->zRes_ << " " << indexNum << " " << -1 << "\n";
				++id;
			}
			++indexNum;
		}
	}
	//paramPack->defaultDir = pathStr.section('/', 0, -2);
	std::cout << "Save CaliberBulidData Completed!" << std::endl << pathStr.toStdString() << std::endl;
}
