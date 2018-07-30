/********************************************************************************
** Form generated from reading UI file 'sparsetracer.ui'
**
** Created by: Qt User Interface Compiler version 5.5.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_SPARSETRACER_H
#define UI_SPARSETRACER_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_SparseTracer
{
public:
    QAction *actionOpenImage;
    QAction *actionSaveTree;
    QAction *actionAuto;
    QAction *actionClear;
    QAction *actionOpenInitSwc;
    QAction *actionRun;
    QAction *actionStop;
    QAction *actionBackTrack;
    QAction *actionSaveImage;
    QAction *actionChoose_Line;
    QAction *actionChoose_Vertex;
    QAction *actionDelete_Line;
    QAction *actionCut_Vertex;
    QAction *actionDraw_Line;
    QAction *action2DView;
    QAction *actionZoom_in;
    QAction *actionZoom_out;
    QAction *actionMove;
    QAction *actionVisible;
    QAction *actionSaveSVM;
    QAction *actionSaveBack;
    QAction *actionNGTree;
    QAction *actionSBWT;
    QAction *actionPick_Soma;
    QAction *actionDelete_Soma;
    QAction *actionSelect_Tree;
    QAction *actionNGTree_Trace;
    QAction *actionDelete_Tree;
    QAction *actionTest;
    QAction *actionTrain;
    QAction *actionTreeChecker;
    QAction *actionLocal_Run;
    QAction *actionSaveLayerImage;
    QAction *actionSaveLayerSwc;
    QAction *actionStart;
    QAction *actionHalt;
    QAction *actionReset;
    QAction *actionNeuroGPS;
    QAction *actionSaveSoma;
    QAction *actionOpenSoma;
    QAction *actionSnapshot;
    QAction *actionCreate_HDF5;
    QAction *actionSaveCaliber;
    QWidget *centralWidget;
    QGridLayout *gridLayout;
    QHBoxLayout *mainLayout;
    QMenuBar *menuBar;
    QMenu *menu_File;
    QMenu *menuTrace;
    QMenu *menuEdit;
    QMenu *menuMode;
    QMenu *menuDebug;
    QMenu *menuTimer_2;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;
    QToolBar *toolBar;
    QToolBar *toolBar_2;

    void setupUi(QMainWindow *SparseTracer)
    {
        if (SparseTracer->objectName().isEmpty())
            SparseTracer->setObjectName(QStringLiteral("SparseTracer"));
        SparseTracer->resize(706, 584);
        QIcon icon;
        icon.addFile(QStringLiteral(":/new/tracer/Resource/Science.png"), QSize(), QIcon::Normal, QIcon::Off);
        SparseTracer->setWindowIcon(icon);
        actionOpenImage = new QAction(SparseTracer);
        actionOpenImage->setObjectName(QStringLiteral("actionOpenImage"));
        QIcon icon1;
        icon1.addFile(QStringLiteral(":/new/tracer/Resource/open.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionOpenImage->setIcon(icon1);
        actionSaveTree = new QAction(SparseTracer);
        actionSaveTree->setObjectName(QStringLiteral("actionSaveTree"));
        QIcon icon2;
        icon2.addFile(QStringLiteral(":/new/tracer/Resource/save.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionSaveTree->setIcon(icon2);
        actionAuto = new QAction(SparseTracer);
        actionAuto->setObjectName(QStringLiteral("actionAuto"));
        QIcon icon3;
        icon3.addFile(QStringLiteral(":/new/tracer/Resource/go-next.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionAuto->setIcon(icon3);
        actionClear = new QAction(SparseTracer);
        actionClear->setObjectName(QStringLiteral("actionClear"));
        QIcon icon4;
        icon4.addFile(QStringLiteral(":/new/tracer/Resource/sweep.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionClear->setIcon(icon4);
        actionOpenInitSwc = new QAction(SparseTracer);
        actionOpenInitSwc->setObjectName(QStringLiteral("actionOpenInitSwc"));
        QIcon icon5;
        icon5.addFile(QStringLiteral(":/new/tracer/Resource/seed.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionOpenInitSwc->setIcon(icon5);
        actionRun = new QAction(SparseTracer);
        actionRun->setObjectName(QStringLiteral("actionRun"));
        QIcon icon6;
        icon6.addFile(QStringLiteral(":/new/tracer/Resource/run.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionRun->setIcon(icon6);
        actionStop = new QAction(SparseTracer);
        actionStop->setObjectName(QStringLiteral("actionStop"));
        QIcon icon7;
        icon7.addFile(QStringLiteral(":/new/tracer/Resource/stop.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionStop->setIcon(icon7);
        actionBackTrack = new QAction(SparseTracer);
        actionBackTrack->setObjectName(QStringLiteral("actionBackTrack"));
        QIcon icon8;
        icon8.addFile(QStringLiteral(":/new/tracer/Resource/reload.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionBackTrack->setIcon(icon8);
        actionSaveImage = new QAction(SparseTracer);
        actionSaveImage->setObjectName(QStringLiteral("actionSaveImage"));
        actionChoose_Line = new QAction(SparseTracer);
        actionChoose_Line->setObjectName(QStringLiteral("actionChoose_Line"));
        actionChoose_Line->setCheckable(true);
        QIcon icon9;
        icon9.addFile(QStringLiteral(":/new/tracer/Resource/chooseline.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionChoose_Line->setIcon(icon9);
        actionChoose_Vertex = new QAction(SparseTracer);
        actionChoose_Vertex->setObjectName(QStringLiteral("actionChoose_Vertex"));
        actionChoose_Vertex->setCheckable(true);
        QIcon icon10;
        icon10.addFile(QStringLiteral(":/new/tracer/Resource/choosepoint.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionChoose_Vertex->setIcon(icon10);
        actionDelete_Line = new QAction(SparseTracer);
        actionDelete_Line->setObjectName(QStringLiteral("actionDelete_Line"));
        QIcon icon11;
        icon11.addFile(QStringLiteral(":/new/tracer/Resource/deleteline.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionDelete_Line->setIcon(icon11);
        actionCut_Vertex = new QAction(SparseTracer);
        actionCut_Vertex->setObjectName(QStringLiteral("actionCut_Vertex"));
        QIcon icon12;
        icon12.addFile(QStringLiteral(":/new/tracer/Resource/cut.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionCut_Vertex->setIcon(icon12);
        actionDraw_Line = new QAction(SparseTracer);
        actionDraw_Line->setObjectName(QStringLiteral("actionDraw_Line"));
        actionDraw_Line->setCheckable(true);
        QIcon icon13;
        icon13.addFile(QStringLiteral(":/new/tracer/Resource/drawline.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionDraw_Line->setIcon(icon13);
        action2DView = new QAction(SparseTracer);
        action2DView->setObjectName(QStringLiteral("action2DView"));
        action2DView->setCheckable(true);
        QIcon icon14;
        icon14.addFile(QStringLiteral(":/new/tracer/Resource/view.png"), QSize(), QIcon::Normal, QIcon::Off);
        action2DView->setIcon(icon14);
        actionZoom_in = new QAction(SparseTracer);
        actionZoom_in->setObjectName(QStringLiteral("actionZoom_in"));
        QIcon icon15;
        icon15.addFile(QStringLiteral(":/new/tracer/Resource/zoomin.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionZoom_in->setIcon(icon15);
        actionZoom_out = new QAction(SparseTracer);
        actionZoom_out->setObjectName(QStringLiteral("actionZoom_out"));
        QIcon icon16;
        icon16.addFile(QStringLiteral(":/new/tracer/Resource/zoomout.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionZoom_out->setIcon(icon16);
        actionMove = new QAction(SparseTracer);
        actionMove->setObjectName(QStringLiteral("actionMove"));
        actionMove->setCheckable(true);
        QIcon icon17;
        icon17.addFile(QStringLiteral(":/new/tracer/Resource/move.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionMove->setIcon(icon17);
        actionVisible = new QAction(SparseTracer);
        actionVisible->setObjectName(QStringLiteral("actionVisible"));
        actionVisible->setCheckable(true);
        actionVisible->setChecked(true);
        QIcon icon18;
        icon18.addFile(QStringLiteral(":/new/tracer/Resource/visible.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionVisible->setIcon(icon18);
        actionSaveSVM = new QAction(SparseTracer);
        actionSaveSVM->setObjectName(QStringLiteral("actionSaveSVM"));
        actionSaveBack = new QAction(SparseTracer);
        actionSaveBack->setObjectName(QStringLiteral("actionSaveBack"));
        actionNGTree = new QAction(SparseTracer);
        actionNGTree->setObjectName(QStringLiteral("actionNGTree"));
        actionNGTree->setCheckable(true);
        actionNGTree->setChecked(false);
        actionSBWT = new QAction(SparseTracer);
        actionSBWT->setObjectName(QStringLiteral("actionSBWT"));
        actionSBWT->setCheckable(true);
        actionSBWT->setChecked(true);
        actionPick_Soma = new QAction(SparseTracer);
        actionPick_Soma->setObjectName(QStringLiteral("actionPick_Soma"));
        actionPick_Soma->setCheckable(true);
        QIcon icon19;
        icon19.addFile(QStringLiteral(":/new/tracer/Resource/choosesoma.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionPick_Soma->setIcon(icon19);
        actionDelete_Soma = new QAction(SparseTracer);
        actionDelete_Soma->setObjectName(QStringLiteral("actionDelete_Soma"));
        actionDelete_Soma->setCheckable(false);
        QIcon icon20;
        icon20.addFile(QStringLiteral(":/new/tracer/Resource/deletesoma.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionDelete_Soma->setIcon(icon20);
        actionSelect_Tree = new QAction(SparseTracer);
        actionSelect_Tree->setObjectName(QStringLiteral("actionSelect_Tree"));
        actionSelect_Tree->setCheckable(true);
        QIcon icon21;
        icon21.addFile(QStringLiteral(":/new/tracer/Resource/choosetree.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionSelect_Tree->setIcon(icon21);
        actionNGTree_Trace = new QAction(SparseTracer);
        actionNGTree_Trace->setObjectName(QStringLiteral("actionNGTree_Trace"));
        QIcon icon22;
        icon22.addFile(QStringLiteral(":/new/tracer/Resource/NGTree.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionNGTree_Trace->setIcon(icon22);
        actionDelete_Tree = new QAction(SparseTracer);
        actionDelete_Tree->setObjectName(QStringLiteral("actionDelete_Tree"));
        QIcon icon23;
        icon23.addFile(QStringLiteral(":/new/tracer/Resource/deletetree.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionDelete_Tree->setIcon(icon23);
        actionTest = new QAction(SparseTracer);
        actionTest->setObjectName(QStringLiteral("actionTest"));
        actionTrain = new QAction(SparseTracer);
        actionTrain->setObjectName(QStringLiteral("actionTrain"));
        actionTreeChecker = new QAction(SparseTracer);
        actionTreeChecker->setObjectName(QStringLiteral("actionTreeChecker"));
        QIcon icon24;
        icon24.addFile(QStringLiteral(":/new/tracer/Resource/tree.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionTreeChecker->setIcon(icon24);
        actionLocal_Run = new QAction(SparseTracer);
        actionLocal_Run->setObjectName(QStringLiteral("actionLocal_Run"));
        QIcon icon25;
        icon25.addFile(QStringLiteral(":/new/tracer/Resource/step.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionLocal_Run->setIcon(icon25);
        actionSaveLayerImage = new QAction(SparseTracer);
        actionSaveLayerImage->setObjectName(QStringLiteral("actionSaveLayerImage"));
        actionSaveLayerSwc = new QAction(SparseTracer);
        actionSaveLayerSwc->setObjectName(QStringLiteral("actionSaveLayerSwc"));
        actionStart = new QAction(SparseTracer);
        actionStart->setObjectName(QStringLiteral("actionStart"));
        actionHalt = new QAction(SparseTracer);
        actionHalt->setObjectName(QStringLiteral("actionHalt"));
        actionReset = new QAction(SparseTracer);
        actionReset->setObjectName(QStringLiteral("actionReset"));
        actionNeuroGPS = new QAction(SparseTracer);
        actionNeuroGPS->setObjectName(QStringLiteral("actionNeuroGPS"));
        actionNeuroGPS->setIcon(icon5);
        actionSaveSoma = new QAction(SparseTracer);
        actionSaveSoma->setObjectName(QStringLiteral("actionSaveSoma"));
        actionOpenSoma = new QAction(SparseTracer);
        actionOpenSoma->setObjectName(QStringLiteral("actionOpenSoma"));
        actionSnapshot = new QAction(SparseTracer);
        actionSnapshot->setObjectName(QStringLiteral("actionSnapshot"));
        QIcon icon26;
        icon26.addFile(QStringLiteral(":/new/tracer/Resource/snap.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionSnapshot->setIcon(icon26);
        actionCreate_HDF5 = new QAction(SparseTracer);
        actionCreate_HDF5->setObjectName(QStringLiteral("actionCreate_HDF5"));
        actionSaveCaliber = new QAction(SparseTracer);
        actionSaveCaliber->setObjectName(QStringLiteral("actionSaveCaliber"));
        centralWidget = new QWidget(SparseTracer);
        centralWidget->setObjectName(QStringLiteral("centralWidget"));
        gridLayout = new QGridLayout(centralWidget);
        gridLayout->setSpacing(6);
        gridLayout->setContentsMargins(11, 11, 11, 11);
        gridLayout->setObjectName(QStringLiteral("gridLayout"));
        mainLayout = new QHBoxLayout();
        mainLayout->setSpacing(6);
        mainLayout->setObjectName(QStringLiteral("mainLayout"));

        gridLayout->addLayout(mainLayout, 0, 0, 1, 1);

        SparseTracer->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(SparseTracer);
        menuBar->setObjectName(QStringLiteral("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 706, 23));
        menu_File = new QMenu(menuBar);
        menu_File->setObjectName(QStringLiteral("menu_File"));
        menuTrace = new QMenu(menuBar);
        menuTrace->setObjectName(QStringLiteral("menuTrace"));
        menuEdit = new QMenu(menuBar);
        menuEdit->setObjectName(QStringLiteral("menuEdit"));
        menuMode = new QMenu(menuBar);
        menuMode->setObjectName(QStringLiteral("menuMode"));
        menuDebug = new QMenu(menuBar);
        menuDebug->setObjectName(QStringLiteral("menuDebug"));
        menuTimer_2 = new QMenu(menuDebug);
        menuTimer_2->setObjectName(QStringLiteral("menuTimer_2"));
        SparseTracer->setMenuBar(menuBar);
        mainToolBar = new QToolBar(SparseTracer);
        mainToolBar->setObjectName(QStringLiteral("mainToolBar"));
        SparseTracer->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(SparseTracer);
        statusBar->setObjectName(QStringLiteral("statusBar"));
        SparseTracer->setStatusBar(statusBar);
        toolBar = new QToolBar(SparseTracer);
        toolBar->setObjectName(QStringLiteral("toolBar"));
        SparseTracer->addToolBar(Qt::TopToolBarArea, toolBar);
        toolBar_2 = new QToolBar(SparseTracer);
        toolBar_2->setObjectName(QStringLiteral("toolBar_2"));
        SparseTracer->addToolBar(Qt::TopToolBarArea, toolBar_2);

        menuBar->addAction(menu_File->menuAction());
        menuBar->addAction(menuEdit->menuAction());
        menuBar->addAction(menuTrace->menuAction());
        menuBar->addAction(menuMode->menuAction());
        menuBar->addAction(menuDebug->menuAction());
        menu_File->addAction(actionOpenImage);
        menu_File->addAction(actionOpenSoma);
        menu_File->addAction(actionSaveSoma);
        menu_File->addAction(actionSaveTree);
        menu_File->addAction(actionSaveImage);
        menu_File->addAction(actionSaveCaliber);
        menu_File->addAction(actionClear);
        menu_File->addAction(actionSnapshot);
        menuTrace->addAction(actionNeuroGPS);
        menuTrace->addAction(actionNGTree_Trace);
        menuTrace->addAction(actionLocal_Run);
        menuTrace->addAction(actionTreeChecker);
        menuTrace->addAction(actionRun);
        menuTrace->addAction(actionTrain);
        menuTrace->addAction(actionStop);
        menuEdit->addAction(action2DView);
        menuEdit->addAction(actionZoom_in);
        menuEdit->addAction(actionZoom_out);
        menuEdit->addAction(actionVisible);
        menuEdit->addSeparator();
        menuEdit->addSeparator();
        menuEdit->addAction(actionCreate_HDF5);
        menuMode->addAction(actionNGTree);
        menuMode->addAction(actionSBWT);
        menuDebug->addAction(actionSaveSVM);
        menuDebug->addAction(actionSaveBack);
        menuDebug->addAction(actionSaveLayerImage);
        menuDebug->addAction(actionSaveLayerSwc);
        menuDebug->addAction(actionChoose_Line);
        menuDebug->addAction(actionChoose_Vertex);
        menuDebug->addAction(actionDelete_Line);
        menuDebug->addAction(actionCut_Vertex);
        menuDebug->addAction(actionDraw_Line);
        menuDebug->addAction(actionPick_Soma);
        menuDebug->addAction(actionDelete_Soma);
        menuDebug->addAction(actionSelect_Tree);
        menuDebug->addAction(actionDelete_Tree);
        menuDebug->addAction(actionTest);
        menuDebug->addAction(menuTimer_2->menuAction());
        menuTimer_2->addAction(actionStart);
        menuTimer_2->addAction(actionHalt);
        menuTimer_2->addAction(actionReset);
        mainToolBar->addAction(actionOpenImage);
        mainToolBar->addAction(actionSaveTree);
        mainToolBar->addAction(actionSnapshot);
        mainToolBar->addAction(actionClear);
        mainToolBar->addSeparator();
        mainToolBar->addAction(actionNGTree);
        mainToolBar->addAction(actionSBWT);
        mainToolBar->addSeparator();
        mainToolBar->addAction(actionZoom_in);
        mainToolBar->addAction(actionZoom_out);
        mainToolBar->addAction(actionVisible);
        toolBar->addSeparator();
        toolBar->addAction(actionRun);
        toolBar->addAction(actionStop);
        toolBar->addAction(actionTrain);
        toolBar_2->addAction(actionNGTree_Trace);
        toolBar_2->addAction(actionLocal_Run);
        toolBar_2->addAction(actionTreeChecker);
        toolBar_2->addAction(actionNeuroGPS);

        retranslateUi(SparseTracer);

        QMetaObject::connectSlotsByName(SparseTracer);
    } // setupUi

    void retranslateUi(QMainWindow *SparseTracer)
    {
        SparseTracer->setWindowTitle(QApplication::translate("SparseTracer", "GTree", 0));
        actionOpenImage->setText(QApplication::translate("SparseTracer", "OpenImage", 0));
        actionOpenImage->setShortcut(QApplication::translate("SparseTracer", "Ctrl+O", 0));
        actionSaveTree->setText(QApplication::translate("SparseTracer", "SaveTree", 0));
        actionSaveTree->setShortcut(QApplication::translate("SparseTracer", "Ctrl+S", 0));
        actionAuto->setText(QApplication::translate("SparseTracer", "Auto", 0));
        actionClear->setText(QApplication::translate("SparseTracer", "Clear", 0));
        actionOpenInitSwc->setText(QApplication::translate("SparseTracer", "OpenInitSwc", 0));
        actionRun->setText(QApplication::translate("SparseTracer", "Run", 0));
        actionRun->setShortcut(QApplication::translate("SparseTracer", "T", 0));
        actionStop->setText(QApplication::translate("SparseTracer", "Stop", 0));
        actionBackTrack->setText(QApplication::translate("SparseTracer", "BackTrack", 0));
        actionSaveImage->setText(QApplication::translate("SparseTracer", "SaveImage", 0));
        actionChoose_Line->setText(QApplication::translate("SparseTracer", "Choose Line", 0));
        actionChoose_Vertex->setText(QApplication::translate("SparseTracer", "Choose Vertex", 0));
        actionDelete_Line->setText(QApplication::translate("SparseTracer", "Delete Line", 0));
        actionCut_Vertex->setText(QApplication::translate("SparseTracer", "Cut Vertex", 0));
        actionDraw_Line->setText(QApplication::translate("SparseTracer", "Draw Line", 0));
        action2DView->setText(QApplication::translate("SparseTracer", "2DView", 0));
        action2DView->setShortcut(QApplication::translate("SparseTracer", "2", 0));
        actionZoom_in->setText(QApplication::translate("SparseTracer", "Zoom in", 0));
        actionZoom_out->setText(QApplication::translate("SparseTracer", "Zoom out", 0));
        actionMove->setText(QApplication::translate("SparseTracer", "Move", 0));
        actionVisible->setText(QApplication::translate("SparseTracer", "Visible", 0));
        actionVisible->setShortcut(QApplication::translate("SparseTracer", "A", 0));
        actionSaveSVM->setText(QApplication::translate("SparseTracer", "SaveSVM", 0));
        actionSaveBack->setText(QApplication::translate("SparseTracer", "SaveBack", 0));
        actionNGTree->setText(QApplication::translate("SparseTracer", "NGTree", 0));
        actionSBWT->setText(QApplication::translate("SparseTracer", "SBWT", 0));
        actionPick_Soma->setText(QApplication::translate("SparseTracer", "Pick Soma", 0));
        actionDelete_Soma->setText(QApplication::translate("SparseTracer", "Delete Soma", 0));
        actionSelect_Tree->setText(QApplication::translate("SparseTracer", "Select Tree", 0));
        actionNGTree_Trace->setText(QApplication::translate("SparseTracer", "NGTree Trace", 0));
        actionDelete_Tree->setText(QApplication::translate("SparseTracer", "Delete Tree", 0));
        actionTest->setText(QApplication::translate("SparseTracer", "Test", 0));
        actionTrain->setText(QApplication::translate("SparseTracer", "Train", 0));
        actionTreeChecker->setText(QApplication::translate("SparseTracer", "TreeChecker", 0));
        actionLocal_Run->setText(QApplication::translate("SparseTracer", "Local Run", 0));
        actionSaveLayerImage->setText(QApplication::translate("SparseTracer", "SaveLayerImage", 0));
        actionSaveLayerSwc->setText(QApplication::translate("SparseTracer", "SaveLayerSwc", 0));
        actionStart->setText(QApplication::translate("SparseTracer", "Start", 0));
        actionHalt->setText(QApplication::translate("SparseTracer", "Halt", 0));
        actionReset->setText(QApplication::translate("SparseTracer", "Reset", 0));
        actionNeuroGPS->setText(QApplication::translate("SparseTracer", "NeuroGPS", 0));
        actionSaveSoma->setText(QApplication::translate("SparseTracer", "SaveSoma", 0));
        actionOpenSoma->setText(QApplication::translate("SparseTracer", "OpenSoma", 0));
        actionSnapshot->setText(QApplication::translate("SparseTracer", "Snapshot", 0));
        actionCreate_HDF5->setText(QApplication::translate("SparseTracer", "Create HDF5", 0));
        actionSaveCaliber->setText(QApplication::translate("SparseTracer", "SaveCaliber", 0));
        menu_File->setTitle(QApplication::translate("SparseTracer", "&File", 0));
        menuTrace->setTitle(QApplication::translate("SparseTracer", "Trace", 0));
        menuEdit->setTitle(QApplication::translate("SparseTracer", "Edit", 0));
        menuMode->setTitle(QApplication::translate("SparseTracer", "Mode", 0));
        menuDebug->setTitle(QApplication::translate("SparseTracer", "Debug", 0));
        menuTimer_2->setTitle(QApplication::translate("SparseTracer", "Timer", 0));
        toolBar->setWindowTitle(QApplication::translate("SparseTracer", "toolBar", 0));
        toolBar_2->setWindowTitle(QApplication::translate("SparseTracer", "toolBar_2", 0));
#ifndef QT_NO_TOOLTIP
        toolBar_2->setToolTip(QApplication::translate("SparseTracer", "locate soma", 0));
#endif // QT_NO_TOOLTIP
    } // retranslateUi

};

namespace Ui {
    class SparseTracer: public Ui_SparseTracer {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_SPARSETRACER_H
