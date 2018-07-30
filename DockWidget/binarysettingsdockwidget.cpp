#include "binarysettingsdockwidget.h"
#include "ui_binarysettingsdockwidget.h"
#include <omp.h>
#include <QSpinBox>
#include <QMessageBox>
#include <QPainter>
#include <QLabel>
#include "../ngtypes/tree.h"
#include "../ngtypes/volume.h"

#include <time.h>
#include <QDebug>//for test

BinarySettingsDockWidget::BinarySettingsDockWidget(QWidget *parent, NGParamPack p):
    QDockWidget(parent),
	paramPack(p),
    ui(new Ui::BinarySettingsDockWidget)
{
    ui->setupUi(this);
    ui->applyImageReadButton->setEnabled(false);
    connect(ui->lowOpacSpinBox, SIGNAL(valueChanged(int)),this,SLOT(OpacValueChanged(int)));
    connect(ui->highOpacSpinBox, SIGNAL(valueChanged(int)),this,SLOT(OpacValueChanged(int)));
    xyzSpinBoxSignalControl_ = false;
    auto *tab2Layout = ui->tab_4->layout();

	//grayRange = 255;//2018-2-10
	
    treeWidget = new NGTreeWidget(ui->tab_4);
    //treeWidget->SetParam(paramPack);
    checkWidget = new NGTreeWidget(ui->tab_4, NGTreeWidget::CHECKWIDGETMODE);
    //checkWidget->SetParam(paramPack);
    buttonWidget_1 = new QWidget(this);
    tab2Layout->addWidget(buttonWidget_1);
    QGridLayout *hLayout1 = new QGridLayout(buttonWidget_1);
    //tab2Layout->addItem(hLayout);
    startTraverse_ = new QPushButton("start traverse",this);
    startTraverse_->setCheckable(true);
    annotateLabel_ = new QPushButton("annotate", this);
    backButton_ = new QPushButton("back", this);
    nextButton_ = new QPushButton("next", this);
    hLayout1->addWidget(annotateLabel_, 0,1,1,1);
    hLayout1->addWidget(startTraverse_, 0, 0, 1, 1);
    hLayout1->addWidget(backButton_ , 1, 0, 1, 1);
    hLayout1->addWidget(nextButton_, 1, 1, 1, 1);
    tab2Layout->addWidget(treeWidget);
    tab2Layout->addWidget(checkWidget);

    QTreeWidgetItem *___qtreewidgetitem = treeWidget->headerItem();
    //___qtreewidgetitem->setText(2, QApplication::translate("BinarySettingsDockWidget", "active", 0));
    ___qtreewidgetitem->setText(0, QApplication::translate("BinarySettingsDockWidget", "ID", 0));
    ___qtreewidgetitem->setText(1, QApplication::translate("BinarySettingsDockWidget", "active", 0));
    QTreeWidgetItem *___qtreewidgetitem1 = checkWidget->headerItem();
    ___qtreewidgetitem1->setText(3, QApplication::translate("BinarySettingsDockWidget", "z", 0));
    ___qtreewidgetitem1->setText(2, QApplication::translate("BinarySettingsDockWidget", "y", 0));
    ___qtreewidgetitem1->setText(1, QApplication::translate("BinarySettingsDockWidget", "x", 0));
    ___qtreewidgetitem1->setText(0, QApplication::translate("BinarySettingsDockWidget", "checked", 0));

    treeWidget->setColumnWidth(0, 50);
    treeWidget->setColumnWidth(1, 20);
    //treeWidget->setColumnWidth(2, 20);
    checkWidget->setColumnWidth(1, 20);
    checkWidget->setColumnWidth(2, 20);
    checkWidget->setColumnWidth(3, 20);

	//108 03 29
	ui->SmoothSlider->setMinimum(1);
	ui->SmoothSlider->setMaximum(100);
	ui->SmoothSlider->setValue(1);
	ui->SmoothSlider->setSingleStep(1);

    connect(treeWidget, SIGNAL(ActiveTree_Signal(int)), this, SLOT(ActivateTree_Slot(int)));
    connect(treeWidget, SIGNAL(CompareTree_Signal(int)), this, SLOT(CompareTree_Slot(int)));
    connect(treeWidget, SIGNAL(ToggleTreeVisible_Signal(int)), this, SLOT(ToggleTreeVisible_Slot(int)));
    connect(treeWidget, SIGNAL(ResetTraverse_Signal()), this, SLOT(ResetTraverse_Slot()));
    connect(checkWidget, SIGNAL(GotoDiff_Signal(int)), this, SLOT(GotoDiff_Slot(int)));
    connect(checkWidget, SIGNAL(ToggleShowDiffArrow_Signal(bool)), this, SLOT(ToggleShowDiffArrow_Slot(bool)));

    connect(startTraverse_, SIGNAL(toggled(bool)), this, SLOT(ToggleTraverse_Slot(bool)));
    connect(nextButton_, SIGNAL(clicked()), this, SLOT(NextTraverse_Slot()));
    connect(backButton_, SIGNAL(clicked()), this, SLOT(BackTraverse_Slot()));
    connect(annotateLabel_, SIGNAL(clicked()), this, SLOT(Annotate_Slot()));

    connect(treeWidget, SIGNAL(ToggleShowTraverseFlag_Signal()), this, SLOT(ToggleShowTraverseFlag_Slot()));
}

BinarySettingsDockWidget::~BinarySettingsDockWidget()
{
    delete ui;
}

DEFINE_SETFUNC(BinarySettingsDockWidget, SetAxonDiffuseValue, axonDiffuseValue_, axonDiffuseValueSpinBox, double);
DEFINE_SETFUNC(BinarySettingsDockWidget, SetAxonTraceValue, axonTraceValue_, axonTraceValueSpinBox, double);
DEFINE_SETFUNC(BinarySettingsDockWidget, SetDiffuseValue, diffuseValue_, diffuseValueSpinBox, double);
DEFINE_SETFUNC(BinarySettingsDockWidget, SetThickness_Slot, thicknesss_, thicknessspinBox, int);
DEFINE_SETFUNC(BinarySettingsDockWidget, SetTraceValue, traceValue_, traceValueSpinBox, double);
DEFINE_SETFUNC(BinarySettingsDockWidget, SetMaxBoundNum, maxBoundNum_, maxBoundNumSpinBox, int);
DEFINE_SETFUNC(BinarySettingsDockWidget, SetAxonSemiautoTraceThreshold, axonSemiautoTraceThreshold, axonSemiautoTraceTresholdSpinBox, double);//2018-1-22
DEFINE_SETFUNC(BinarySettingsDockWidget, SetAxonSemiautoConnectResample, axonSemiautoConnectResampleThreshold, axonSemiautoConnectResampleSpinBox, double);//2018-2-9
///////////////////////////////////////////////////////

DEFINE_EDITINGFINISHED_SIGNAL(BinarySettingsDockWidget, xMinspinBox, xMin_, ROIChanged_Signal());
DEFINE_EDITINGFINISHED_SIGNAL(BinarySettingsDockWidget, xMaxspinBox, xMax_, ROIChanged_Signal());
DEFINE_EDITINGFINISHED_SIGNAL(BinarySettingsDockWidget, yMinspinBox, yMin_, ROIChanged_Signal());
DEFINE_EDITINGFINISHED_SIGNAL(BinarySettingsDockWidget, yMaxspinBox, yMax_, ROIChanged_Signal());
DEFINE_EDITINGFINISHED_SIGNAL(BinarySettingsDockWidget, zMinspinBox, zMin_, ROIChanged_Signal());
DEFINE_EDITINGFINISHED_SIGNAL(BinarySettingsDockWidget, zMaxspinBox, zMax_, ROIChanged_Signal());

DEFINE_VALUECHANGED(BinarySettingsDockWidget, axonDiffuseValueSpinBox, axonDiffuseValue_, double);
DEFINE_VALUECHANGED(BinarySettingsDockWidget, axonThresholdSpinBox, axonBinaryThreshold_, double);
DEFINE_VALUECHANGED(BinarySettingsDockWidget, axonTraceValueSpinBox, axonTraceValue_, double);
DEFINE_VALUECHANGED(BinarySettingsDockWidget, axonSemiautoTraceTresholdSpinBox, axonSemiautoTraceThreshold, double);//2018-1-22
DEFINE_VALUECHANGED(BinarySettingsDockWidget, axonSemiautoConnectResampleSpinBox, axonSemiautoConnectResampleThreshold, double);//2018-2-9
DEFINE_VALUECHANGED(BinarySettingsDockWidget, diffuseValueSpinBox, diffuseValue_, double);
DEFINE_VALUECHANGED(BinarySettingsDockWidget, highOpacSpinBox, highOpac_, int);
DEFINE_VALUECHANGED(BinarySettingsDockWidget, lowOpacSpinBox, lowOpac_, int);
DEFINE_VALUECHANGED(BinarySettingsDockWidget, thresholdSpinBox, binThreshold_, double);
DEFINE_VALUECHANGED(BinarySettingsDockWidget, traceValueSpinBox, traceValue_, double);
DEFINE_VALUECHANGED(BinarySettingsDockWidget, xBlockSizespinBox, xExtract_, int);
DEFINE_VALUECHANGED(BinarySettingsDockWidget, yBlockSizespinBox, yExtract_, int);
DEFINE_VALUECHANGED(BinarySettingsDockWidget, zBlockSizespinBox, zExtract_, int);
DEFINE_VALUECHANGED(BinarySettingsDockWidget, xScaleSpinBox, xScale_, int);
DEFINE_VALUECHANGED(BinarySettingsDockWidget, yScaleSpinBox, yScale_, int);
DEFINE_VALUECHANGED(BinarySettingsDockWidget, zScaleSpinBox, zScale_, int);
DEFINE_VALUECHANGED(BinarySettingsDockWidget, origLowOpacSpinBox, lowAdjustOpac_, int);
DEFINE_VALUECHANGED(BinarySettingsDockWidget, origHighOpacSpinBox, highAdjustOpac_, int);
DEFINE_VALUECHANGED(BinarySettingsDockWidget, destLowOpacSpinBox, lowDestOpac_, int);
DEFINE_VALUECHANGED(BinarySettingsDockWidget, destHighOpacSpinBox, highDestOpac_, int);


void BinarySettingsDockWidget::on_maxBoundNumSpinBox_editingFinished()
{
    if (paramPack->maxBoundNum_ != ui->maxBoundNumSpinBox->value()) {
        paramPack->maxBoundNum_ = ui->maxBoundNumSpinBox->value();
        emit MaxBoundNumChanged_Signal(paramPack->maxBoundNum_);
    }
}

void BinarySettingsDockWidget::SetMaxThickness(int arg)
{
    ui->thicknessspinBox->setMaximum(arg);
}

void BinarySettingsDockWidget::closeEvent(QCloseEvent * event)
{
    event->ignore();
}

void BinarySettingsDockWidget::on_applyImageReadButton_clicked()
{
    paramPack->mostdLevel_ = ui->levelSpinBox->value();
    emit ApplyReadImage_Signal();
}

//hehe
void BinarySettingsDockWidget::OpacValueChanged(int)
{
    paramPack->lowOpac_ = ui->lowOpacSpinBox->value();
    paramPack->highOpac_ = ui->highOpacSpinBox->value();
    emit DisplayNewOpac_Signal();
}

//
void BinarySettingsDockWidget::UpdateImageROILabel()
{
    xyzSpinBoxSignalControl_ = true;
    ui->xMinspinBox->setValue(paramPack->xMin_);
    ui->xMaxspinBox->setValue(paramPack->xMax_);
    ui->yMinspinBox->setValue(paramPack->yMin_);
    ui->yMaxspinBox->setValue(paramPack->yMax_);
    ui->zMinspinBox->setValue(paramPack->zMin_);
    ui->zMaxspinBox->setValue(paramPack->zMax_);
    xyzSpinBoxSignalControl_ = false;
}

void BinarySettingsDockWidget::UpdateImageROILabel(double *arg)
{
    if (!arg) return;
    xyzSpinBoxSignalControl_ = true;
    ui->xMinspinBox->setValue(arg[0]);
    ui->xMaxspinBox->setValue(arg[1]);
    ui->yMinspinBox->setValue(arg[2]);
    ui->yMaxspinBox->setValue(arg[3]);
    ui->zMinspinBox->setValue(arg[4]);
    ui->zMaxspinBox->setValue(arg[5]);
    xyzSpinBoxSignalControl_ = false;
}

void BinarySettingsDockWidget::on_thicknessspinBox_valueChanged(int arg)
{
    paramPack->thicknesss_ = arg;
    emit ThicknessChanged_Signal(arg);
}

void BinarySettingsDockWidget::SetDisplayOpacValue(int low, int high)
{
    paramPack->lowOpac_ = low;
    paramPack->highOpac_ = high;
    ui->lowOpacSpinBox->setValue(low);
    ui->highOpacSpinBox->setValue(high);
}

void BinarySettingsDockWidget::SetApplyReadImageButton(bool arg)
{
    ui->applyImageReadButton->setEnabled(arg);
}

void BinarySettingsDockWidget::SetImageRange(int x, int y, int z)
{
    ui->imageRangeLabel->setText( tr("Image Size: \n%1 x %2 x %3").arg(x).arg(y).arg(z));
    ui->xMinspinBox->setMaximum(x);
    ui->yMinspinBox->setMaximum(y);
    ui->zMinspinBox->setMaximum(z);
}

void BinarySettingsDockWidget::SetParam( NGParamPack arg )
{
    paramPack = arg;
    treeWidget->SetParam(paramPack);
    checkWidget->SetParam(paramPack);
}

void BinarySettingsDockWidget::UpdateImageExtractLabel()
{
    ui->xBlockSizespinBox->setValue(paramPack->xExtract_);
    ui->yBlockSizespinBox->setValue(paramPack->yExtract_);
    ui->zBlockSizespinBox->setValue(paramPack->zExtract_);
}

void BinarySettingsDockWidget::UpdateParam()
{
    if (!paramPack) return;
    paramPack->xMin_ = ui->xMinspinBox->value();
    paramPack->xMax_ = ui->xMaxspinBox->value();
    paramPack->yMin_ = ui->yMinspinBox->value();
    paramPack->yMax_ = ui->yMaxspinBox->value();
    paramPack->zMin_ = ui->zMinspinBox->value();
    paramPack->zMax_ = ui->zMaxspinBox->value();
    paramPack->binThreshold_ = ui->thresholdSpinBox->value();
    paramPack->axonBinaryThreshold_ = ui->axonThresholdSpinBox->value();
    paramPack->diffuseValue_ = ui->diffuseValueSpinBox->value();
    paramPack->axonDiffuseValue_ = ui->axonDiffuseValueSpinBox->value();
    paramPack->traceValue_ = ui->traceValueSpinBox->value();
    paramPack->axonTraceValue_ = ui->axonTraceValueSpinBox->value();
    paramPack->maxBoundNum_ = ui->maxBoundNumSpinBox->value();
	paramPack->axonSemiautoTraceThreshold = ui->axonSemiautoTraceTresholdSpinBox->value();//2018-1-22
	paramPack->axonSemiautoConnectResampleThreshold = ui->axonSemiautoConnectResampleSpinBox->value();//2018-2-9
    paramPack->thicknesss_ = ui->thicknessspinBox->value();
    paramPack->lowOpac_ = ui->lowOpacSpinBox->value();
    paramPack->highOpac_ = ui->highOpacSpinBox->value();
    paramPack->mostdLevel_ = ui->levelSpinBox->value();
    paramPack->enableSVM_ = ui->enableSVMBox->isChecked();
    paramPack->strongSignalMode_ = ui->strongSignalCheckBox->isChecked();
	paramPack->smooth_level = (double)(ui->SmoothSlider->value())/100;//2018 03 29
}

void BinarySettingsDockWidget::UpdateGUI()
{
    ui->highOpacSpinBox->setValue(paramPack->highOpac_);
    ui->lowOpacSpinBox->setValue(paramPack->lowOpac_);
    UpdateImageROILabel();
    ui->thresholdSpinBox->setValue(paramPack->binThreshold_);
    ui->axonThresholdSpinBox->setValue(paramPack->axonBinaryThreshold_);
	ui->axonSemiautoTraceTresholdSpinBox->setValue(paramPack->axonSemiautoTraceThreshold);//2018-1-22
	ui->axonSemiautoConnectResampleSpinBox->setValue(paramPack->axonSemiautoConnectResampleThreshold);//2018-2-9
    ui->traceValueSpinBox->setValue(paramPack->traceValue_);
    ui->diffuseValueSpinBox->setValue(paramPack->diffuseValue_);
    ui->axonDiffuseValueSpinBox->setValue(paramPack->axonDiffuseValue_);
    ui->axonTraceValueSpinBox->setValue(paramPack->axonTraceValue_);
    ui->maxBoundNumSpinBox->setValue(paramPack->maxBoundNum_);
    ui->levelSpinBox->setValue(paramPack->mostdLevel_);
    SetThickness_Slot(paramPack->thicknesss_);
    ui->xBlockSizespinBox->setValue(paramPack->xExtract_);
    ui->yBlockSizespinBox->setValue(paramPack->yExtract_);
    ui->zBlockSizespinBox->setValue(paramPack->zExtract_);
    ui->highOpacSpinBox->setValue(paramPack->highOpac_);
    ui->lowOpacSpinBox->setValue(paramPack->lowOpac_);
    ui->xScaleSpinBox->setValue(paramPack->xScale_);
    ui->yScaleSpinBox->setValue(paramPack->yScale_);
    ui->destLowOpacSpinBox->setValue(paramPack->lowDestOpac_);
    ui->destHighOpacSpinBox->setValue(paramPack->highDestOpac_);
    ui->origLowOpacSpinBox->setValue(paramPack->lowAdjustOpac_);
    ui->origHighOpacSpinBox->setValue(paramPack->highAdjustOpac_);
	ui->SmoothSlider->setValue(100 * paramPack->smooth_level);

    UpdateTreeWidget();
}

void BinarySettingsDockWidget::GetROIValues(double * arg)
{
    if (!arg) return;
    arg[0] = ui->xMinspinBox->value();
    arg[1] = ui->xMaxspinBox->value();
    arg[2] = ui->yMinspinBox->value();
    arg[3] = ui->yMaxspinBox->value();
    arg[4] = ui->zMinspinBox->value();
    arg[5] = ui->zMaxspinBox->value();
}

void BinarySettingsDockWidget::on_action_ClearCache_Button_clicked()
{
    emit ClearCache_Signal();
}

void BinarySettingsDockWidget::UpdateTreeWidget()
{
    treeWidget->clear();
    NG_CREATE_DYNAMIC_CAST(NeuronPopulation, tmpPop, paramPack->separateTree);
    if (!tmpPop) {
        printf("there is no tree.\n");
        return;
    }
    QStringList tmpitemstr;
    for (size_t k = 0; k < tmpPop->m_pop.size(); ++k) {
        tmpitemstr.clear();
        if (paramPack->activeTree && tmpPop->m_pop[k].get() == paramPack->activeTree.get()) {
            tmpitemstr << QString(tr("tree%1").arg(k)) << "a";
        }
        else tmpitemstr << QString(tr("tree%1").arg(k)) << "d";
        treeWidget->addRoot(tmpitemstr);
    }
}

void BinarySettingsDockWidget::ActivateTree_Slot(int arg)
{
    if (IsActiveTree(arg)) 
        return;
    emit ActiveTree_Signal(arg);
}

void BinarySettingsDockWidget::CompareTree_Slot(int arg)
{
    if (IsActiveTree(arg))
        return;
    emit CompareTree_Signal(arg);
}

void BinarySettingsDockWidget::ToggleTreeVisible_Slot(int id)
{
    emit ToggleTreeVisible_Signal(id);
}

void BinarySettingsDockWidget::GotoDiff_Slot(int arg)
{
    emit GotoDiff_Signal(arg);
}

bool BinarySettingsDockWidget::IsActiveTree(int arg)
{
    NG_CREATE_DYNAMIC_CAST(NeuronPopulation, tmpPop, paramPack->separateTree);
    if (!tmpPop) {
        printf("there is no tree.\n");
        return false;
    }
    for (size_t k = 0; k < tmpPop->m_pop.size(); ++k) {
        if (tmpPop->m_pop[k].get() == paramPack->activeTree.get()) {
            if (arg == int(k)) {
                return true;
            }
        }
    }
    return false;
}

void BinarySettingsDockWidget::UpdateCheckWidget()
{
    if (!paramPack->activeTree) return;
    auto &arg = paramPack->activeTree->m_suspiciousPositions_;
    auto &state = paramPack->activeTree->m_suspiciousCheckState_;
    if (state.size() != arg.size()) {
        NG_ERROR_MESSAGE("");
        return;
    }
    checkWidget->clear();
    QStringList tmpitemstr;
    for (size_t k = 0; k <arg.size(); ++k) {
        tmpitemstr.clear();
        if (state[k] == 0) 
            tmpitemstr << QString(tr("no")) << QString(tr("%1").arg(arg[k](0))) << QString(tr("%1").arg(arg[k](1))) << QString(tr("%1").arg(arg[k](2)));//have not checked
        else if (state[k] == 1)
            tmpitemstr << QString(tr("yes")) << QString(tr("%1").arg(arg[k](0))) << QString(tr("%1").arg(arg[k](1))) << QString(tr("%1").arg(arg[k](2)));//no need to check next
        else//2
            tmpitemstr << QString(tr("tol")) << QString(tr("%1").arg(arg[k](0))) << QString(tr("%1").arg(arg[k](1))) << QString(tr("%1").arg(arg[k](2)));//this wrong is tolerable
        checkWidget->addRoot(tmpitemstr);
    }
}

void BinarySettingsDockWidget::ToggleTraverse_Slot(bool arg)
{
    if (!arg) {
        startTraverse_->setText("start traverse");
    }
    else{
        startTraverse_->setText("stop traverse");
    }
    emit ToggleTraverse_Signal(arg);
}

void BinarySettingsDockWidget::Annotate_Slot()
{
    if (!paramPack->activeTree) return;
    Vec3d susPos;
    susPos(0) = (paramPack->xMin_ + paramPack->xMax_) / 2;
    susPos(1) = (paramPack->yMin_ + paramPack->yMax_) / 2;
    susPos(2) = (paramPack->zMin_ + paramPack->zMax_) / 2;
    auto &allLine = paramPack->activeTree->m_curveList;
    double minDist = 50, tmpDist = 20;
    Vec5d *pt = NULL;
    for (size_t i = 0; i < allLine.size(); ++i) {
        for (size_t j = 0; j < allLine[i].size(); ++j) {
            auto &curPos = allLine[i][j];
            if (std::abs(curPos(0) - susPos(0)) < minDist &&std::abs(curPos(1) - susPos(1)) < minDist &&std::abs(curPos(2) - susPos(2)) < minDist){
                tmpDist = std::abs(curPos(0) - susPos(0)) + std::abs(curPos(1) - susPos(1)) + std::abs(curPos(2) - susPos(2));
                if ( tmpDist< minDist) {
                    minDist = tmpDist;
                    pt = &curPos;
                }
            }
        }
    }
    if (pt == NULL) {
        QMessageBox::warning(this, "Wrong!", "There is no point in current area.");
        return;
    }
    paramPack->activeTree->m_suspiciousPositions_.push_back(Vec3d((*pt)(0), (*pt)(1), (*pt)(2)));
    paramPack->activeTree->m_suspiciousCheckState_.push_back(0);
    UpdateCheckWidget();
}

void BinarySettingsDockWidget::BackTraverse_Slot()
{
    emit BackTraverse_Signal();
}

void BinarySettingsDockWidget::NextTraverse_Slot()
{
    emit NextTraverse_Signal();
}
//
//void BinarySettingsDockWidget::UpdateSuspiciousWidget()
//{
//    checkWidget->clear();
//    auto &sp = paramPack->activeTree->m_suspiciousPositions_;
//    QStringList tmpitemstr;
//    for (size_t k = 0; k < sp.size(); ++k) {
//        tmpitemstr.clear();
//        tmpitemstr << QString(tr("%1").arg(k)) << QString(tr("%1").arg(sp[k](0))) << QString(tr("%1").arg(sp[k](1))) << QString(tr("%1").arg(sp[k](2)));
//        checkWidget->addRoot(tmpitemstr);
//    }
//}

void BinarySettingsDockWidget::ToggleShowDiffArrow_Slot(bool arg)
{
    emit ToggleShowDiffArrow_Signal(arg);
}

void BinarySettingsDockWidget::ResetTraverse_Slot()
{
    if (!paramPack->separateTree) {
        NG_ERROR_MESSAGE("");
        return;
    }
    NG_CREATE_DYNAMIC_CAST(NeuronPopulation, tmpPop, paramPack->separateTree);
    if (!tmpPop) {
        printf("there is no tree.\n");
        return;
    }
    emit ResetTraverse_Signal();
    //tmpPop->m_pop[arg]->traverseFlag_.clear();
}

void BinarySettingsDockWidget::on_SmoothSlider_valueChanged(int arg1){
	paramPack->smooth_level = arg1*0.01;
}


void BinarySettingsDockWidget::UpdateTreeWidget_Slot()
{
    UpdateTreeWidget();
}

void BinarySettingsDockWidget::UpdateCheckWidget_Slot()
{
    UpdateCheckWidget();
}

void BinarySettingsDockWidget::ToggleShowTraverseFlag_Slot()
{
    emit ToggleShowTraverseFlag_Signal();
}

void BinarySettingsDockWidget::on_enableSVMBox_toggled(bool arg)
{
    if (arg){
        paramPack->enableSVM_ = true;
    }
    else{
        paramPack->enableSVM_ = false;
    }
}

void BinarySettingsDockWidget::on_GPSSVMcheckBox_toggled(bool arg)
{
	if (arg){
		paramPack->enableGPSSVM_ = true;
	}
	else{
		paramPack->enableGPSSVM_ = false;
	}
}

void BinarySettingsDockWidget::on_strongSignalCheckBox_toggled(bool arg)
{
    if (arg){
        paramPack->strongSignalMode_ = true;
    }
    else{
        paramPack->strongSignalMode_ = false;
    }
}

void BinarySettingsDockWidget::UpdateImageScale()
{
    ui->xScaleSpinBox->setValue(paramPack->xScale_);
    ui->yScaleSpinBox->setValue(paramPack->yScale_);
    ui->zScaleSpinBox->setValue(paramPack->zScale_);
}

void BinarySettingsDockWidget::on_imageOpacInfoButton_clicked()
{
    if (paramPack->OrigImage) {
        NG_CREATE_DYNAMIC_CAST(SVolume, tmpOrig, paramPack->OrigImage);
        if (tmpOrig) {
            tmpOrig->ComputeVolumeMinMax();
            ui->origLowOpacSpinBox->setValue(tmpOrig->MinVal());
            ui->origHighOpacSpinBox->setValue(tmpOrig->MaxVal());
        }
    }
}

void BinarySettingsDockWidget::on_opacApplyButton_clicked()
{
    emit OpacAdjustApply_Signal();
}

void BinarySettingsDockWidget::on_previewIlluminationMapButton_clicked(bool arg)
{
    emit PreviewOpacAdjust_Signal(arg);
}

void BinarySettingsDockWidget::on_previewBinaryButton_clicked(bool arg)
{
    emit PreviewBinary_Signal(arg);
}

void BinarySettingsDockWidget::on_previewAxonBinaryButton_clicked(bool arg)
{
    emit PreviewAxonBinary_Signal(arg);
}

//void BinarySettingsDockWidget::UpdateHistogram()
//{
//	clock_t _start = clock(), _end;//for test
//	qDebug() << "Start to draw.";
//
//	std::shared_ptr<SVolume> oriImg = std::dynamic_pointer_cast<SVolume>(paramPack->OrigImage);
//
//	const unsigned short *temGray1 = (*oriImg).GetPointer();
//	unsigned short grayMin = *temGray1;
//	unsigned short grayMax = *temGray1;
//
//	for (int i = 0; i < oriImg->x()*oriImg->y()*oriImg->z(); ++i)//get max and min gray value
//	{
//		if (*temGray1 > grayMax)
//			grayMax = *temGray1;
//		else if (*temGray1 < grayMin)
//			grayMin = *temGray1;
//
//		++temGray1;
//	}
//
//	grayRange = grayMax - grayMin + 1;
//	VecXd grayArry((int)(grayRange));//the first value of grayArry is the count number of min gray
//	grayArry.setZero();//last is max gray
//
//	const unsigned short *temGray2 = (*oriImg).GetPointer();
//	for (int i = 0; i < oriImg->x()*oriImg->y()*oriImg->z(); ++i)//get the corresponding number of each gray value
//	{
//		++grayArry[(int)(*temGray2 - grayMin)];//has been sorted //no need -1
//		++temGray2;
//	}
//
//	grayArry /= grayArry.maxCoeff();
//	grayArry *= 100;
//
//	//ready to draw
//	_histimage = QImage::QImage(140, 40 + grayRange,QImage::Format_RGB32);//real size
//	_histimage.fill(qRgb(255, 255, 255));
//
//	QPainter _painter(&_histimage);
//	_painter.setBrush(QBrush(QColor(255, 255, 255)));
//	_painter.drawRect(0, 0, _histimage.width(), _histimage.height());
//	_painter.drawLine(20, 20, 20, 20 + grayRange);//y axis
//
//	int yStep = 21;
//	for (int i = 1; i < grayArry.size(); ++i)//i=1,no min
//	{
//		int tem = 20 + (int)(grayArry[i]);
//		_painter.drawLine(20, yStep, tem, yStep);
//		++yStep;
//	}
//	_painter.end();//end to draw histogram
//	_histimage = _histimage.scaled(140, 255, Qt::IgnoreAspectRatio);
//
//	QPainter _painter2(&_histimage);
//	_painter2.drawLine(20, 5, 120, 5);//x axis
//	_painter2.drawText(QPointF(116, 19), QString::number(grayArry.maxCoeff()) + "%");//draw text
//	_painter2.drawText(QPointF(2, 22), QString::number(grayMin));
//	_painter2.drawText(QPointF(2, 245), QString::number(grayMax));
//	_painter2.end();
//
//	ui->histogramLabel->setScaledContents(true);
//	ui->histogramLabel->setPixmap(QPixmap::fromImage(_histimage));
//
//	_end = clock();
//	qDebug() << "This is draw histogram time:" << difftime(_end, _start);
//}
//
//void BinarySettingsDockWidget::mouseMoveEvent(QMouseEvent *event)
//{
//	//need to be modified
//	if (ui->tabWidget->currentIndex() == 1)
//	{
//		if ((paramPack->OrigImage.get()))
//		{
//			//qDebug() << "This is geometry::" << ui->histogramLabel->geometry();
//			qDebug() << event->pos();
//			if (!( event->x() < 9 || event->x() > 149 || event->y() < 240 || event->y() > 495 ))//not out of range
//			{
//				double scaledRate = (40 + grayRange) / 255;//star to translate point
//				double normValue = ((event->y() - 230) * scaledRate - 20) / grayRange;
//				qDebug() << "norm value:" << normValue * 100;
//			}
//		}
//	}
//}
