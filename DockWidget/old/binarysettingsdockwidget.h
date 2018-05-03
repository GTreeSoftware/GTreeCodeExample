/* SparseTracer
*	function: trace single neuron
*	Author : zhouhang
*	2015-10-28
*/
#ifndef BINARYSETTINGSDOCKWIDGET_H
#define BINARYSETTINGSDOCKWIDGET_H

#include <QDockWidget>
#include <QCloseEvent>
#include <QTreeWidget>
#include <QTreeWidgetItem>
#include <QGridLayout>
#include <QPushButton>
#include "../ngtypes/ParamPack.h"
#include "./NGTreeWidget.h"

namespace Ui {
class BinarySettingsDockWidget;
}

#define DEFINE_VALUECHANGED(supername, uiname, varname, type) void supername::on_##uiname##_valueChanged(type arg) \
{ \
    paramPack->varname = arg; \
}

#define DEFINE_EDITINGFINISHED_SIGNAL(supername, uiname, varname, signal) void supername::on_##uiname##_editingFinished() \
{ \
    if(paramPack->varname != ui->uiname->value()){\
        paramPack->varname = ui->uiname->value(); \
        if(!xyzSpinBoxSignalControl_) emit signal; }\
}

#define IMPL_SETFUNC(funcname, type) void funcname(type arg)

#define DEFINE_SETFUNC(supername, funcname, varname, uiname, type) void supername::funcname(type arg) \
{ \
    paramPack->varname = arg; \
    ui->uiname->setValue(arg); \
}

class BinarySettingsDockWidget : public QDockWidget
{
    Q_OBJECT

public:
    explicit BinarySettingsDockWidget(QWidget *parent = 0);
    ~BinarySettingsDockWidget();
    void SetApplyReadImageButton(bool arg);
    void SetAxonDiffuseValue(double);
    void SetAxonTraceValue(double);
    void SetBinaryApplyButtonEnable(bool arg);
    void SetDiffuseValue(double);
    void SetDisplayOpacValue(int, int);
    void SetImageRange(int x, int y, int z);
    void SetMaxThickness(int arg);
    void SetMaxBoundNum(int arg);
    void SetParam(NGParamPack arg);
    void SetTraceValue(double);
    void GetROIValues(double *);
    void UpdateGUI();
    void UpdateImageROILabel();
    void UpdateImageROILabel(double *arg);
    void UpdateParam();
    void UpdateImageExtractLabel();
    void UpdateTreeWidget();
    void UpdateCheckWidget();
signals:
    void ApplyReadImage_Signal();
    void DisplayNewOpac_Signal();
    void MaxBoundNumChanged_Signal(int);
    void NeedBinary_Signal();
    void NeedTrace_Signal();
    void ROIChanged_Signal();
    void ThicknessChanged_Signal(int);
    void ClearCache_Signal();

    void ActiveTree_Signal(int);
    void CompareTree_Signal(int);
    void ToggleTreeVisible_Signal(int);
    void GotoDiff_Signal(int);
    void ToggleTraverse_Signal(bool);
    void Annotate_Signal();
    void BackTraverse_Signal();
    void NextTraverse_Signal();
    void ToggleShowDiffArrow_Signal(bool);
    void ResetTraverse_Signal();
    void ToggleShowTraverseFlag_Signal();

protected:
    void closeEvent(QCloseEvent * event);
    //QTreeWidgetItem* AddTreeNode(QTreeWidgetItem* parent, QString, double x, double y, double z );
    //void CreateActions();
    bool IsActiveTree(int arg);
    //void UpdateSuspiciousWidget();

private slots:
    void on_xMinspinBox_editingFinished();
    void on_xMaxspinBox_editingFinished();
    void on_yMinspinBox_editingFinished();
    void on_yMaxspinBox_editingFinished();
    void on_zMinspinBox_editingFinished();
    void on_zMaxspinBox_editingFinished();

    void on_applyImageReadButton_clicked();
    void on_axonDiffuseValueSpinBox_valueChanged(double);
    void on_axonThresholdSpinBox_valueChanged(double);
    void on_axonTraceValueSpinBox_valueChanged(double);
    void on_diffuseValueSpinBox_valueChanged(double);
    void on_highOpacSpinBox_valueChanged(int);
    void on_lowOpacSpinBox_valueChanged(int);
    void on_maxBoundNumSpinBox_editingFinished();
    void on_thicknessspinBox_valueChanged(int);
    void on_thresholdSpinBox_valueChanged(double);
    void on_traceValueSpinBox_valueChanged(double);
    void OpacValueChanged(int );
    void SetThickness_Slot(int);
	void on_xBlockSizespinBox_valueChanged(int);//radius
	void on_yBlockSizespinBox_valueChanged(int);
	void on_zBlockSizespinBox_valueChanged(int);

    void on_action_ClearCache_Button_clicked();//action_ClearCache_Button

    void ActivateTree_Slot(int);
    void CompareTree_Slot(int);
    void ToggleTreeVisible_Slot(int);
    void GotoDiff_Slot(int);
    void ToggleTraverse_Slot(bool);
    void Annotate_Slot();
    void BackTraverse_Slot();
    void NextTraverse_Slot();

    void ToggleShowDiffArrow_Slot(bool);

    void ResetTraverse_Slot();
    
    void on_SmoothCombox_currentTextChanged(const QString &arg1);

    void on_KernelCombox_currentTextChanged(const QString &arg1);

    void UpdateTreeWidget_Slot();
    void UpdateCheckWidget_Slot();
    void ToggleShowTraverseFlag_Slot();

    void on_enableSVMBox_toggled(bool);
    void on_strongSignalCheckBox_toggled(bool);

private:
    Ui::BinarySettingsDockWidget *ui;
    NGParamPack paramPack;
    bool xyzSpinBoxSignalControl_;
    QPushButton *startTraverse_, *backButton_, *nextButton_, *annotateLabel_;
    QWidget *buttonWidget_1, *buttonWidget_2;
    NGTreeWidget *treeWidget, *checkWidget;
    //QGridLayout *gridLayout;
};

#endif // BINARYSETTINGSDOCKWIDGET_H
