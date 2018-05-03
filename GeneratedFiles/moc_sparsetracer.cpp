/****************************************************************************
** Meta object code from reading C++ file 'sparsetracer.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.5.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../sparsetracer.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'sparsetracer.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.5.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_SparseTracer_t {
    QByteArrayData data[76];
    char stringdata0[1792];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_SparseTracer_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_SparseTracer_t qt_meta_stringdata_SparseTracer = {
    {
QT_MOC_LITERAL(0, 0, 12), // "SparseTracer"
QT_MOC_LITERAL(1, 13, 13), // "renderCaliber"
QT_MOC_LITERAL(2, 27, 0), // ""
QT_MOC_LITERAL(3, 28, 28), // "on_actionOpenImage_triggered"
QT_MOC_LITERAL(4, 57, 27), // "on_actionSaveTree_triggered"
QT_MOC_LITERAL(5, 85, 24), // "on_actionClear_triggered"
QT_MOC_LITERAL(6, 110, 22), // "on_actionRun_triggered"
QT_MOC_LITERAL(7, 133, 12), // "EndRunThread"
QT_MOC_LITERAL(8, 146, 23), // "on_actionStop_triggered"
QT_MOC_LITERAL(9, 170, 17), // "UpdateStatus_Slot"
QT_MOC_LITERAL(10, 188, 22), // "ApplyReadNewImage_Slot"
QT_MOC_LITERAL(11, 211, 28), // "on_actionSaveImage_triggered"
QT_MOC_LITERAL(12, 240, 27), // "on_actionSaveBack_triggered"
QT_MOC_LITERAL(13, 268, 22), // "SetROIByBoxWidget_Slot"
QT_MOC_LITERAL(14, 291, 28), // "on_actionChoose_Line_toggled"
QT_MOC_LITERAL(15, 320, 30), // "on_actionChoose_Vertex_toggled"
QT_MOC_LITERAL(16, 351, 30), // "on_actionDelete_Line_triggered"
QT_MOC_LITERAL(17, 382, 29), // "on_actionCut_Vertex_triggered"
QT_MOC_LITERAL(18, 412, 26), // "on_actionDraw_Line_toggled"
QT_MOC_LITERAL(19, 439, 23), // "on_action2DView_toggled"
QT_MOC_LITERAL(20, 463, 26), // "on_actionZoom_in_triggered"
QT_MOC_LITERAL(21, 490, 27), // "on_actionZoom_out_triggered"
QT_MOC_LITERAL(22, 518, 24), // "on_actionVisible_toggled"
QT_MOC_LITERAL(23, 543, 26), // "on_actionSaveSVM_triggered"
QT_MOC_LITERAL(24, 570, 21), // "ChangeMoveCursor_Slot"
QT_MOC_LITERAL(25, 592, 15), // "TimingSave_Slot"
QT_MOC_LITERAL(26, 608, 26), // "on_actionPick_Soma_toggled"
QT_MOC_LITERAL(27, 635, 3), // "arg"
QT_MOC_LITERAL(28, 639, 28), // "on_actionSelect_Tree_toggled"
QT_MOC_LITERAL(29, 668, 23), // "on_actionNGTree_toggled"
QT_MOC_LITERAL(30, 692, 21), // "on_actionSBWT_toggled"
QT_MOC_LITERAL(31, 714, 31), // "on_actionNGTree_Trace_triggered"
QT_MOC_LITERAL(32, 746, 30), // "on_actionDelete_Soma_triggered"
QT_MOC_LITERAL(33, 777, 30), // "on_actionDelete_Tree_triggered"
QT_MOC_LITERAL(34, 808, 23), // "on_actionTest_triggered"
QT_MOC_LITERAL(35, 832, 25), // "EditActionGroup_triggered"
QT_MOC_LITERAL(36, 858, 8), // "QAction*"
QT_MOC_LITERAL(37, 867, 27), // "SetProjectMaxThickness_Slot"
QT_MOC_LITERAL(38, 895, 23), // "MaxBoundNumChanged_Slot"
QT_MOC_LITERAL(39, 919, 20), // "Set2DViewStatus_Slot"
QT_MOC_LITERAL(40, 940, 18), // "SetDrawStatus_Slot"
QT_MOC_LITERAL(41, 959, 24), // "on_actionTrain_triggered"
QT_MOC_LITERAL(42, 984, 30), // "on_actionTreeChecker_triggered"
QT_MOC_LITERAL(43, 1015, 22), // "UpdateGLBoxWidget_Slot"
QT_MOC_LITERAL(44, 1038, 28), // "on_actionLocal_Run_triggered"
QT_MOC_LITERAL(45, 1067, 33), // "on_actionSaveLayerImage_trigg..."
QT_MOC_LITERAL(46, 1101, 31), // "on_actionSaveLayerSwc_triggered"
QT_MOC_LITERAL(47, 1133, 20), // "ClearMOSTDCache_Slot"
QT_MOC_LITERAL(48, 1154, 17), // "ActivateTree_Slot"
QT_MOC_LITERAL(49, 1172, 16), // "CompareTree_Slot"
QT_MOC_LITERAL(50, 1189, 22), // "ToggleTreeVisible_Slot"
QT_MOC_LITERAL(51, 1212, 13), // "GotoDiff_Slot"
QT_MOC_LITERAL(52, 1226, 19), // "ToggleTraverse_Slot"
QT_MOC_LITERAL(53, 1246, 17), // "BackTraverse_Slot"
QT_MOC_LITERAL(54, 1264, 17), // "NextTraverse_Slot"
QT_MOC_LITERAL(55, 1282, 24), // "ToggleShowDiffArrow_Slot"
QT_MOC_LITERAL(56, 1307, 26), // "StartTraverseFromHere_Slot"
QT_MOC_LITERAL(57, 1334, 18), // "ResetTraverse_Slot"
QT_MOC_LITERAL(58, 1353, 24), // "on_actionStart_triggered"
QT_MOC_LITERAL(59, 1378, 23), // "on_actionHalt_triggered"
QT_MOC_LITERAL(60, 1402, 24), // "on_actionReset_triggered"
QT_MOC_LITERAL(61, 1427, 19), // "AddRunningTime_Slot"
QT_MOC_LITERAL(62, 1447, 18), // "CacheComplete_Slot"
QT_MOC_LITERAL(63, 1466, 27), // "on_actionNeuroGPS_triggered"
QT_MOC_LITERAL(64, 1494, 27), // "on_actionSaveSoma_triggered"
QT_MOC_LITERAL(65, 1522, 27), // "on_actionOpenSoma_triggered"
QT_MOC_LITERAL(66, 1550, 30), // "on_actionSaveCaliber_triggered"
QT_MOC_LITERAL(67, 1581, 20), // "OpacAdjustApply_Slot"
QT_MOC_LITERAL(68, 1602, 22), // "PreviewOpacAdjust_Slot"
QT_MOC_LITERAL(69, 1625, 18), // "PreviewBinary_Slot"
QT_MOC_LITERAL(70, 1644, 22), // "PreviewAxonBinary_Slot"
QT_MOC_LITERAL(71, 1667, 27), // "on_actionSnapshot_triggered"
QT_MOC_LITERAL(72, 1695, 30), // "on_actionCreate_HDF5_triggered"
QT_MOC_LITERAL(73, 1726, 18), // "EndCreateHDF5_Slot"
QT_MOC_LITERAL(74, 1745, 26), // "UpdateHDF5ProgressBar_Slot"
QT_MOC_LITERAL(75, 1772, 19) // "StopCreateHDF5_Slot"

    },
    "SparseTracer\0renderCaliber\0\0"
    "on_actionOpenImage_triggered\0"
    "on_actionSaveTree_triggered\0"
    "on_actionClear_triggered\0"
    "on_actionRun_triggered\0EndRunThread\0"
    "on_actionStop_triggered\0UpdateStatus_Slot\0"
    "ApplyReadNewImage_Slot\0"
    "on_actionSaveImage_triggered\0"
    "on_actionSaveBack_triggered\0"
    "SetROIByBoxWidget_Slot\0"
    "on_actionChoose_Line_toggled\0"
    "on_actionChoose_Vertex_toggled\0"
    "on_actionDelete_Line_triggered\0"
    "on_actionCut_Vertex_triggered\0"
    "on_actionDraw_Line_toggled\0"
    "on_action2DView_toggled\0"
    "on_actionZoom_in_triggered\0"
    "on_actionZoom_out_triggered\0"
    "on_actionVisible_toggled\0"
    "on_actionSaveSVM_triggered\0"
    "ChangeMoveCursor_Slot\0TimingSave_Slot\0"
    "on_actionPick_Soma_toggled\0arg\0"
    "on_actionSelect_Tree_toggled\0"
    "on_actionNGTree_toggled\0on_actionSBWT_toggled\0"
    "on_actionNGTree_Trace_triggered\0"
    "on_actionDelete_Soma_triggered\0"
    "on_actionDelete_Tree_triggered\0"
    "on_actionTest_triggered\0"
    "EditActionGroup_triggered\0QAction*\0"
    "SetProjectMaxThickness_Slot\0"
    "MaxBoundNumChanged_Slot\0Set2DViewStatus_Slot\0"
    "SetDrawStatus_Slot\0on_actionTrain_triggered\0"
    "on_actionTreeChecker_triggered\0"
    "UpdateGLBoxWidget_Slot\0"
    "on_actionLocal_Run_triggered\0"
    "on_actionSaveLayerImage_triggered\0"
    "on_actionSaveLayerSwc_triggered\0"
    "ClearMOSTDCache_Slot\0ActivateTree_Slot\0"
    "CompareTree_Slot\0ToggleTreeVisible_Slot\0"
    "GotoDiff_Slot\0ToggleTraverse_Slot\0"
    "BackTraverse_Slot\0NextTraverse_Slot\0"
    "ToggleShowDiffArrow_Slot\0"
    "StartTraverseFromHere_Slot\0"
    "ResetTraverse_Slot\0on_actionStart_triggered\0"
    "on_actionHalt_triggered\0"
    "on_actionReset_triggered\0AddRunningTime_Slot\0"
    "CacheComplete_Slot\0on_actionNeuroGPS_triggered\0"
    "on_actionSaveSoma_triggered\0"
    "on_actionOpenSoma_triggered\0"
    "on_actionSaveCaliber_triggered\0"
    "OpacAdjustApply_Slot\0PreviewOpacAdjust_Slot\0"
    "PreviewBinary_Slot\0PreviewAxonBinary_Slot\0"
    "on_actionSnapshot_triggered\0"
    "on_actionCreate_HDF5_triggered\0"
    "EndCreateHDF5_Slot\0UpdateHDF5ProgressBar_Slot\0"
    "StopCreateHDF5_Slot"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_SparseTracer[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
      72,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    0,  374,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
       3,    0,  375,    2, 0x08 /* Private */,
       4,    0,  376,    2, 0x08 /* Private */,
       5,    0,  377,    2, 0x08 /* Private */,
       6,    0,  378,    2, 0x08 /* Private */,
       7,    0,  379,    2, 0x08 /* Private */,
       8,    0,  380,    2, 0x08 /* Private */,
       9,    0,  381,    2, 0x08 /* Private */,
      10,    0,  382,    2, 0x08 /* Private */,
      11,    0,  383,    2, 0x08 /* Private */,
      12,    0,  384,    2, 0x08 /* Private */,
      13,    0,  385,    2, 0x08 /* Private */,
      14,    1,  386,    2, 0x08 /* Private */,
      15,    1,  389,    2, 0x08 /* Private */,
      16,    0,  392,    2, 0x08 /* Private */,
      17,    0,  393,    2, 0x08 /* Private */,
      18,    1,  394,    2, 0x08 /* Private */,
      19,    1,  397,    2, 0x08 /* Private */,
      20,    0,  400,    2, 0x08 /* Private */,
      21,    0,  401,    2, 0x08 /* Private */,
      22,    1,  402,    2, 0x08 /* Private */,
      23,    0,  405,    2, 0x08 /* Private */,
      24,    1,  406,    2, 0x08 /* Private */,
      25,    0,  409,    2, 0x08 /* Private */,
      26,    1,  410,    2, 0x08 /* Private */,
      28,    1,  413,    2, 0x08 /* Private */,
      29,    1,  416,    2, 0x08 /* Private */,
      30,    1,  419,    2, 0x08 /* Private */,
      31,    0,  422,    2, 0x08 /* Private */,
      32,    0,  423,    2, 0x08 /* Private */,
      33,    0,  424,    2, 0x08 /* Private */,
      34,    0,  425,    2, 0x08 /* Private */,
      35,    1,  426,    2, 0x08 /* Private */,
      37,    1,  429,    2, 0x08 /* Private */,
      38,    1,  432,    2, 0x08 /* Private */,
      39,    1,  435,    2, 0x08 /* Private */,
      40,    1,  438,    2, 0x08 /* Private */,
      41,    0,  441,    2, 0x08 /* Private */,
      42,    0,  442,    2, 0x08 /* Private */,
      43,    0,  443,    2, 0x08 /* Private */,
      44,    0,  444,    2, 0x08 /* Private */,
      45,    0,  445,    2, 0x08 /* Private */,
      46,    0,  446,    2, 0x08 /* Private */,
      47,    0,  447,    2, 0x08 /* Private */,
      48,    1,  448,    2, 0x08 /* Private */,
      49,    1,  451,    2, 0x08 /* Private */,
      50,    1,  454,    2, 0x08 /* Private */,
      51,    1,  457,    2, 0x08 /* Private */,
      52,    1,  460,    2, 0x08 /* Private */,
      53,    0,  463,    2, 0x08 /* Private */,
      54,    0,  464,    2, 0x08 /* Private */,
      55,    1,  465,    2, 0x08 /* Private */,
      56,    2,  468,    2, 0x08 /* Private */,
      57,    0,  473,    2, 0x08 /* Private */,
      58,    0,  474,    2, 0x08 /* Private */,
      59,    0,  475,    2, 0x08 /* Private */,
      60,    0,  476,    2, 0x08 /* Private */,
      61,    0,  477,    2, 0x08 /* Private */,
      62,    0,  478,    2, 0x08 /* Private */,
      63,    0,  479,    2, 0x08 /* Private */,
      64,    0,  480,    2, 0x08 /* Private */,
      65,    0,  481,    2, 0x08 /* Private */,
      66,    0,  482,    2, 0x08 /* Private */,
      67,    0,  483,    2, 0x08 /* Private */,
      68,    1,  484,    2, 0x08 /* Private */,
      69,    1,  487,    2, 0x08 /* Private */,
      70,    1,  490,    2, 0x08 /* Private */,
      71,    0,  493,    2, 0x08 /* Private */,
      72,    0,  494,    2, 0x08 /* Private */,
      73,    0,  495,    2, 0x08 /* Private */,
      74,    1,  496,    2, 0x08 /* Private */,
      75,    0,  499,    2, 0x08 /* Private */,

 // signals: parameters
    QMetaType::Void,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Bool,   27,
    QMetaType::Void, QMetaType::Bool,   27,
    QMetaType::Void, QMetaType::Bool,   27,
    QMetaType::Void, QMetaType::Bool,   27,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, 0x80000000 | 36,    2,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int,   27,
    QMetaType::Void, QMetaType::Int,   27,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void, QMetaType::Int, QMetaType::Int,    2,    2,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void,

       0        // eod
};

void SparseTracer::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        SparseTracer *_t = static_cast<SparseTracer *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->renderCaliber(); break;
        case 1: _t->on_actionOpenImage_triggered(); break;
        case 2: _t->on_actionSaveTree_triggered(); break;
        case 3: _t->on_actionClear_triggered(); break;
        case 4: _t->on_actionRun_triggered(); break;
        case 5: _t->EndRunThread(); break;
        case 6: _t->on_actionStop_triggered(); break;
        case 7: _t->UpdateStatus_Slot(); break;
        case 8: _t->ApplyReadNewImage_Slot(); break;
        case 9: _t->on_actionSaveImage_triggered(); break;
        case 10: _t->on_actionSaveBack_triggered(); break;
        case 11: _t->SetROIByBoxWidget_Slot(); break;
        case 12: _t->on_actionChoose_Line_toggled((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 13: _t->on_actionChoose_Vertex_toggled((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 14: _t->on_actionDelete_Line_triggered(); break;
        case 15: _t->on_actionCut_Vertex_triggered(); break;
        case 16: _t->on_actionDraw_Line_toggled((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 17: _t->on_action2DView_toggled((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 18: _t->on_actionZoom_in_triggered(); break;
        case 19: _t->on_actionZoom_out_triggered(); break;
        case 20: _t->on_actionVisible_toggled((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 21: _t->on_actionSaveSVM_triggered(); break;
        case 22: _t->ChangeMoveCursor_Slot((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 23: _t->TimingSave_Slot(); break;
        case 24: _t->on_actionPick_Soma_toggled((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 25: _t->on_actionSelect_Tree_toggled((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 26: _t->on_actionNGTree_toggled((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 27: _t->on_actionSBWT_toggled((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 28: _t->on_actionNGTree_Trace_triggered(); break;
        case 29: _t->on_actionDelete_Soma_triggered(); break;
        case 30: _t->on_actionDelete_Tree_triggered(); break;
        case 31: _t->on_actionTest_triggered(); break;
        case 32: _t->EditActionGroup_triggered((*reinterpret_cast< QAction*(*)>(_a[1]))); break;
        case 33: _t->SetProjectMaxThickness_Slot((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 34: _t->MaxBoundNumChanged_Slot((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 35: _t->Set2DViewStatus_Slot((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 36: _t->SetDrawStatus_Slot((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 37: _t->on_actionTrain_triggered(); break;
        case 38: _t->on_actionTreeChecker_triggered(); break;
        case 39: _t->UpdateGLBoxWidget_Slot(); break;
        case 40: _t->on_actionLocal_Run_triggered(); break;
        case 41: _t->on_actionSaveLayerImage_triggered(); break;
        case 42: _t->on_actionSaveLayerSwc_triggered(); break;
        case 43: _t->ClearMOSTDCache_Slot(); break;
        case 44: _t->ActivateTree_Slot((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 45: _t->CompareTree_Slot((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 46: _t->ToggleTreeVisible_Slot((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 47: _t->GotoDiff_Slot((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 48: _t->ToggleTraverse_Slot((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 49: _t->BackTraverse_Slot(); break;
        case 50: _t->NextTraverse_Slot(); break;
        case 51: _t->ToggleShowDiffArrow_Slot((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 52: _t->StartTraverseFromHere_Slot((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        case 53: _t->ResetTraverse_Slot(); break;
        case 54: _t->on_actionStart_triggered(); break;
        case 55: _t->on_actionHalt_triggered(); break;
        case 56: _t->on_actionReset_triggered(); break;
        case 57: _t->AddRunningTime_Slot(); break;
        case 58: _t->CacheComplete_Slot(); break;
        case 59: _t->on_actionNeuroGPS_triggered(); break;
        case 60: _t->on_actionSaveSoma_triggered(); break;
        case 61: _t->on_actionOpenSoma_triggered(); break;
        case 62: _t->on_actionSaveCaliber_triggered(); break;
        case 63: _t->OpacAdjustApply_Slot(); break;
        case 64: _t->PreviewOpacAdjust_Slot((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 65: _t->PreviewBinary_Slot((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 66: _t->PreviewAxonBinary_Slot((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 67: _t->on_actionSnapshot_triggered(); break;
        case 68: _t->on_actionCreate_HDF5_triggered(); break;
        case 69: _t->EndCreateHDF5_Slot(); break;
        case 70: _t->UpdateHDF5ProgressBar_Slot((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 71: _t->StopCreateHDF5_Slot(); break;
        default: ;
        }
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        switch (_id) {
        default: *reinterpret_cast<int*>(_a[0]) = -1; break;
        case 32:
            switch (*reinterpret_cast<int*>(_a[1])) {
            default: *reinterpret_cast<int*>(_a[0]) = -1; break;
            case 0:
                *reinterpret_cast<int*>(_a[0]) = qRegisterMetaType< QAction* >(); break;
            }
            break;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        void **func = reinterpret_cast<void **>(_a[1]);
        {
            typedef void (SparseTracer::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&SparseTracer::renderCaliber)) {
                *result = 0;
            }
        }
    }
}

const QMetaObject SparseTracer::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_SparseTracer.data,
      qt_meta_data_SparseTracer,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *SparseTracer::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *SparseTracer::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_SparseTracer.stringdata0))
        return static_cast<void*>(const_cast< SparseTracer*>(this));
    return QMainWindow::qt_metacast(_clname);
}

int SparseTracer::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 72)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 72;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 72)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 72;
    }
    return _id;
}

// SIGNAL 0
void SparseTracer::renderCaliber()
{
    QMetaObject::activate(this, &staticMetaObject, 0, Q_NULLPTR);
}
QT_END_MOC_NAMESPACE
