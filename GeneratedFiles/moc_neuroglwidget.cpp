/****************************************************************************
** Meta object code from reading C++ file 'neuroglwidget.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.5.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../Render/neuroglwidget.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'neuroglwidget.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.5.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_NeuroGLWidget_t {
    QByteArrayData data[73];
    char stringdata0[1364];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_NeuroGLWidget_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_NeuroGLWidget_t qt_meta_stringdata_NeuroGLWidget = {
    {
QT_MOC_LITERAL(0, 0, 13), // "NeuroGLWidget"
QT_MOC_LITERAL(1, 14, 23), // "BoxWidgetChanged_Signal"
QT_MOC_LITERAL(2, 38, 0), // ""
QT_MOC_LITERAL(3, 39, 23), // "ChangeMoveCursor_Signal"
QT_MOC_LITERAL(4, 63, 30), // "clippingDistanceChanged_Signal"
QT_MOC_LITERAL(5, 94, 26), // "InteractModeChanged_Signal"
QT_MOC_LITERAL(6, 121, 19), // "Set2DProject_Signal"
QT_MOC_LITERAL(7, 141, 16), // "Set2DView_Signal"
QT_MOC_LITERAL(8, 158, 14), // "SetDraw_Signal"
QT_MOC_LITERAL(9, 173, 21), // "SetTraversePos_Signal"
QT_MOC_LITERAL(10, 195, 20), // "NeedUpdateTreeWidget"
QT_MOC_LITERAL(11, 216, 21), // "NeedUpdateCheckWidget"
QT_MOC_LITERAL(12, 238, 19), // "NextTraverse_Signal"
QT_MOC_LITERAL(13, 258, 12), // "DeleteBranch"
QT_MOC_LITERAL(14, 271, 14), // "SeparateBranch"
QT_MOC_LITERAL(15, 286, 10), // "DeleteSoma"
QT_MOC_LITERAL(16, 297, 11), // "DeleteTrees"
QT_MOC_LITERAL(17, 309, 10), // "DivideLine"
QT_MOC_LITERAL(18, 320, 18), // "MergeTwoPickedTree"
QT_MOC_LITERAL(19, 339, 17), // "SaveSelectedTrees"
QT_MOC_LITERAL(20, 357, 18), // "SetActiveTree_Slot"
QT_MOC_LITERAL(21, 376, 24), // "SetClippingDistance_Slot"
QT_MOC_LITERAL(22, 401, 3), // "arg"
QT_MOC_LITERAL(23, 405, 17), // "SetHeadAsTreeRoot"
QT_MOC_LITERAL(24, 423, 22), // "SetInteractEnable_Slot"
QT_MOC_LITERAL(25, 446, 17), // "SetTailAsTreeRoot"
QT_MOC_LITERAL(26, 464, 14), // "Toggle2D3DView"
QT_MOC_LITERAL(27, 479, 12), // "Toggle2DDraw"
QT_MOC_LITERAL(28, 492, 14), // "ToggleLineMode"
QT_MOC_LITERAL(29, 507, 15), // "TogglePointMode"
QT_MOC_LITERAL(30, 523, 31), // "ToggleShowCurrentLinesDirection"
QT_MOC_LITERAL(31, 555, 18), // "ToggleShowTreeRoot"
QT_MOC_LITERAL(32, 574, 14), // "ToggleSomaMode"
QT_MOC_LITERAL(33, 589, 14), // "ToggleTreeMode"
QT_MOC_LITERAL(34, 604, 18), // "ToggleTreesVisible"
QT_MOC_LITERAL(35, 623, 27), // "UnDeleteOneActiveTreeBranch"
QT_MOC_LITERAL(36, 651, 8), // "UndoDraw"
QT_MOC_LITERAL(37, 660, 20), // "UpdateBoxWidget_Slot"
QT_MOC_LITERAL(38, 681, 32), // "UpdateBoxWidgetWithoutReset_Slot"
QT_MOC_LITERAL(39, 714, 7), // "double*"
QT_MOC_LITERAL(40, 722, 21), // "UpdateVolumeOpac_Slot"
QT_MOC_LITERAL(41, 744, 12), // "XYReset_Slot"
QT_MOC_LITERAL(42, 757, 12), // "XZReset_Slot"
QT_MOC_LITERAL(43, 770, 12), // "YZReset_Slot"
QT_MOC_LITERAL(44, 783, 17), // "ShowPosition_Slot"
QT_MOC_LITERAL(45, 801, 14), // "ConjChildCurve"
QT_MOC_LITERAL(46, 816, 22), // "ToggleEnableCurveTrace"
QT_MOC_LITERAL(47, 839, 14), // "ToggleShowInit"
QT_MOC_LITERAL(48, 854, 13), // "UndeleteTrees"
QT_MOC_LITERAL(49, 868, 16), // "ToggleBoundCheck"
QT_MOC_LITERAL(50, 885, 18), // "ForceBoundCheckAll"
QT_MOC_LITERAL(51, 904, 21), // "ToggleSmartSemiTracer"
QT_MOC_LITERAL(52, 926, 28), // "ToggleActiveTreeLayerDisplay"
QT_MOC_LITERAL(53, 955, 15), // "Annotation_Slot"
QT_MOC_LITERAL(54, 971, 17), // "Unannotation_Slot"
QT_MOC_LITERAL(55, 989, 17), // "ShowLevelBranchID"
QT_MOC_LITERAL(56, 1007, 26), // "StartTraverseFromHere_Slot"
QT_MOC_LITERAL(57, 1034, 22), // "InsertPointBefore_Slot"
QT_MOC_LITERAL(58, 1057, 16), // "RemovePoint_Slot"
QT_MOC_LITERAL(59, 1074, 20), // "ToggleMovePoint_Slot"
QT_MOC_LITERAL(60, 1095, 22), // "ToggleReconParent_Slot"
QT_MOC_LITERAL(61, 1118, 16), // "RefreshTree_Slot"
QT_MOC_LITERAL(62, 1135, 17), // "BuildCaliber_Slot"
QT_MOC_LITERAL(63, 1153, 23), // "SetBranchAsDenrite_Slot"
QT_MOC_LITERAL(64, 1177, 20), // "SetBranchAsAxon_Slot"
QT_MOC_LITERAL(65, 1198, 17), // "HideDendrite_Slot"
QT_MOC_LITERAL(66, 1216, 13), // "HideAxon_Slot"
QT_MOC_LITERAL(67, 1230, 16), // "HideCaliber_Slot"
QT_MOC_LITERAL(68, 1247, 23), // "SmoothCurrentCurve_Slot"
QT_MOC_LITERAL(69, 1271, 27), // "SaveImageAndSwcForTest_Slot"
QT_MOC_LITERAL(70, 1299, 27), // "ToggleShowTraverseFlag_Slot"
QT_MOC_LITERAL(71, 1327, 18), // "renderCaliber_Slot"
QT_MOC_LITERAL(72, 1346, 17) // "AllCaliberBulider"

    },
    "NeuroGLWidget\0BoxWidgetChanged_Signal\0"
    "\0ChangeMoveCursor_Signal\0"
    "clippingDistanceChanged_Signal\0"
    "InteractModeChanged_Signal\0"
    "Set2DProject_Signal\0Set2DView_Signal\0"
    "SetDraw_Signal\0SetTraversePos_Signal\0"
    "NeedUpdateTreeWidget\0NeedUpdateCheckWidget\0"
    "NextTraverse_Signal\0DeleteBranch\0"
    "SeparateBranch\0DeleteSoma\0DeleteTrees\0"
    "DivideLine\0MergeTwoPickedTree\0"
    "SaveSelectedTrees\0SetActiveTree_Slot\0"
    "SetClippingDistance_Slot\0arg\0"
    "SetHeadAsTreeRoot\0SetInteractEnable_Slot\0"
    "SetTailAsTreeRoot\0Toggle2D3DView\0"
    "Toggle2DDraw\0ToggleLineMode\0TogglePointMode\0"
    "ToggleShowCurrentLinesDirection\0"
    "ToggleShowTreeRoot\0ToggleSomaMode\0"
    "ToggleTreeMode\0ToggleTreesVisible\0"
    "UnDeleteOneActiveTreeBranch\0UndoDraw\0"
    "UpdateBoxWidget_Slot\0"
    "UpdateBoxWidgetWithoutReset_Slot\0"
    "double*\0UpdateVolumeOpac_Slot\0"
    "XYReset_Slot\0XZReset_Slot\0YZReset_Slot\0"
    "ShowPosition_Slot\0ConjChildCurve\0"
    "ToggleEnableCurveTrace\0ToggleShowInit\0"
    "UndeleteTrees\0ToggleBoundCheck\0"
    "ForceBoundCheckAll\0ToggleSmartSemiTracer\0"
    "ToggleActiveTreeLayerDisplay\0"
    "Annotation_Slot\0Unannotation_Slot\0"
    "ShowLevelBranchID\0StartTraverseFromHere_Slot\0"
    "InsertPointBefore_Slot\0RemovePoint_Slot\0"
    "ToggleMovePoint_Slot\0ToggleReconParent_Slot\0"
    "RefreshTree_Slot\0BuildCaliber_Slot\0"
    "SetBranchAsDenrite_Slot\0SetBranchAsAxon_Slot\0"
    "HideDendrite_Slot\0HideAxon_Slot\0"
    "HideCaliber_Slot\0SmoothCurrentCurve_Slot\0"
    "SaveImageAndSwcForTest_Slot\0"
    "ToggleShowTraverseFlag_Slot\0"
    "renderCaliber_Slot\0AllCaliberBulider"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_NeuroGLWidget[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
      69,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
      11,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    0,  359,    2, 0x06 /* Public */,
       3,    1,  360,    2, 0x06 /* Public */,
       4,    1,  363,    2, 0x06 /* Public */,
       5,    1,  366,    2, 0x06 /* Public */,
       6,    1,  369,    2, 0x06 /* Public */,
       7,    1,  372,    2, 0x06 /* Public */,
       8,    1,  375,    2, 0x06 /* Public */,
       9,    2,  378,    2, 0x06 /* Public */,
      10,    0,  383,    2, 0x06 /* Public */,
      11,    0,  384,    2, 0x06 /* Public */,
      12,    0,  385,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
      13,    0,  386,    2, 0x0a /* Public */,
      14,    0,  387,    2, 0x0a /* Public */,
      15,    0,  388,    2, 0x0a /* Public */,
      16,    0,  389,    2, 0x0a /* Public */,
      17,    0,  390,    2, 0x0a /* Public */,
      18,    0,  391,    2, 0x0a /* Public */,
      19,    0,  392,    2, 0x0a /* Public */,
      20,    0,  393,    2, 0x0a /* Public */,
      21,    1,  394,    2, 0x0a /* Public */,
      23,    0,  397,    2, 0x0a /* Public */,
      24,    1,  398,    2, 0x0a /* Public */,
      25,    0,  401,    2, 0x0a /* Public */,
      26,    1,  402,    2, 0x0a /* Public */,
      27,    1,  405,    2, 0x0a /* Public */,
      28,    1,  408,    2, 0x0a /* Public */,
      29,    1,  411,    2, 0x0a /* Public */,
      30,    1,  414,    2, 0x0a /* Public */,
      31,    1,  417,    2, 0x0a /* Public */,
      32,    1,  420,    2, 0x0a /* Public */,
      33,    1,  423,    2, 0x0a /* Public */,
      34,    1,  426,    2, 0x0a /* Public */,
      35,    0,  429,    2, 0x0a /* Public */,
      36,    0,  430,    2, 0x0a /* Public */,
      37,    0,  431,    2, 0x0a /* Public */,
      38,    1,  432,    2, 0x0a /* Public */,
      40,    0,  435,    2, 0x0a /* Public */,
      41,    0,  436,    2, 0x0a /* Public */,
      42,    0,  437,    2, 0x0a /* Public */,
      43,    0,  438,    2, 0x0a /* Public */,
      44,    0,  439,    2, 0x0a /* Public */,
      45,    0,  440,    2, 0x0a /* Public */,
      46,    0,  441,    2, 0x0a /* Public */,
      47,    1,  442,    2, 0x0a /* Public */,
      48,    0,  445,    2, 0x0a /* Public */,
      49,    1,  446,    2, 0x0a /* Public */,
      50,    0,  449,    2, 0x0a /* Public */,
      51,    1,  450,    2, 0x0a /* Public */,
      52,    1,  453,    2, 0x0a /* Public */,
      53,    0,  456,    2, 0x0a /* Public */,
      54,    0,  457,    2, 0x0a /* Public */,
      55,    0,  458,    2, 0x0a /* Public */,
      56,    0,  459,    2, 0x0a /* Public */,
      57,    0,  460,    2, 0x0a /* Public */,
      58,    0,  461,    2, 0x0a /* Public */,
      59,    1,  462,    2, 0x0a /* Public */,
      60,    1,  465,    2, 0x0a /* Public */,
      61,    0,  468,    2, 0x0a /* Public */,
      62,    0,  469,    2, 0x0a /* Public */,
      63,    0,  470,    2, 0x0a /* Public */,
      64,    0,  471,    2, 0x0a /* Public */,
      65,    1,  472,    2, 0x0a /* Public */,
      66,    1,  475,    2, 0x0a /* Public */,
      67,    1,  478,    2, 0x0a /* Public */,
      68,    0,  481,    2, 0x0a /* Public */,
      69,    0,  482,    2, 0x0a /* Public */,
      70,    0,  483,    2, 0x0a /* Public */,
      71,    0,  484,    2, 0x0a /* Public */,
      72,    0,  485,    2, 0x0a /* Public */,

 // signals: parameters
    QMetaType::Void,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void, QMetaType::Int, QMetaType::Int,    2,    2,
    QMetaType::Void,
    QMetaType::Void,
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
    QMetaType::Void, QMetaType::Int,   22,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Bool,   22,
    QMetaType::Void, QMetaType::Bool,   22,
    QMetaType::Void, QMetaType::Bool,   22,
    QMetaType::Void, QMetaType::Bool,   22,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void, QMetaType::Bool,   22,
    QMetaType::Void, QMetaType::Bool,   22,
    QMetaType::Void, QMetaType::Bool,   22,
    QMetaType::Void, QMetaType::Bool,   22,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, 0x80000000 | 39,    2,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Bool,   22,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Bool,   22,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Bool,   22,
    QMetaType::Void, QMetaType::Bool,   22,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Bool, QMetaType::Bool,    2,
    QMetaType::Void, QMetaType::Bool,    2,
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
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void NeuroGLWidget::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        NeuroGLWidget *_t = static_cast<NeuroGLWidget *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->BoxWidgetChanged_Signal(); break;
        case 1: _t->ChangeMoveCursor_Signal((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 2: _t->clippingDistanceChanged_Signal((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 3: _t->InteractModeChanged_Signal((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 4: _t->Set2DProject_Signal((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 5: _t->Set2DView_Signal((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 6: _t->SetDraw_Signal((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 7: _t->SetTraversePos_Signal((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        case 8: _t->NeedUpdateTreeWidget(); break;
        case 9: _t->NeedUpdateCheckWidget(); break;
        case 10: _t->NextTraverse_Signal(); break;
        case 11: _t->DeleteBranch(); break;
        case 12: _t->SeparateBranch(); break;
        case 13: _t->DeleteSoma(); break;
        case 14: _t->DeleteTrees(); break;
        case 15: _t->DivideLine(); break;
        case 16: _t->MergeTwoPickedTree(); break;
        case 17: _t->SaveSelectedTrees(); break;
        case 18: _t->SetActiveTree_Slot(); break;
        case 19: _t->SetClippingDistance_Slot((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 20: _t->SetHeadAsTreeRoot(); break;
        case 21: _t->SetInteractEnable_Slot((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 22: _t->SetTailAsTreeRoot(); break;
        case 23: _t->Toggle2D3DView((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 24: _t->Toggle2DDraw((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 25: _t->ToggleLineMode((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 26: _t->TogglePointMode((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 27: _t->ToggleShowCurrentLinesDirection((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 28: _t->ToggleShowTreeRoot((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 29: _t->ToggleSomaMode((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 30: _t->ToggleTreeMode((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 31: _t->ToggleTreesVisible((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 32: _t->UnDeleteOneActiveTreeBranch(); break;
        case 33: _t->UndoDraw(); break;
        case 34: _t->UpdateBoxWidget_Slot(); break;
        case 35: _t->UpdateBoxWidgetWithoutReset_Slot((*reinterpret_cast< double*(*)>(_a[1]))); break;
        case 36: _t->UpdateVolumeOpac_Slot(); break;
        case 37: _t->XYReset_Slot(); break;
        case 38: _t->XZReset_Slot(); break;
        case 39: _t->YZReset_Slot(); break;
        case 40: _t->ShowPosition_Slot(); break;
        case 41: _t->ConjChildCurve(); break;
        case 42: _t->ToggleEnableCurveTrace(); break;
        case 43: _t->ToggleShowInit((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 44: _t->UndeleteTrees(); break;
        case 45: _t->ToggleBoundCheck((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 46: _t->ForceBoundCheckAll(); break;
        case 47: _t->ToggleSmartSemiTracer((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 48: _t->ToggleActiveTreeLayerDisplay((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 49: _t->Annotation_Slot(); break;
        case 50: _t->Unannotation_Slot(); break;
        case 51: _t->ShowLevelBranchID(); break;
        case 52: _t->StartTraverseFromHere_Slot(); break;
        case 53: _t->InsertPointBefore_Slot(); break;
        case 54: _t->RemovePoint_Slot(); break;
        case 55: { bool _r = _t->ToggleMovePoint_Slot((*reinterpret_cast< bool(*)>(_a[1])));
            if (_a[0]) *reinterpret_cast< bool*>(_a[0]) = _r; }  break;
        case 56: _t->ToggleReconParent_Slot((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 57: _t->RefreshTree_Slot(); break;
        case 58: _t->BuildCaliber_Slot(); break;
        case 59: _t->SetBranchAsDenrite_Slot(); break;
        case 60: _t->SetBranchAsAxon_Slot(); break;
        case 61: _t->HideDendrite_Slot((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 62: _t->HideAxon_Slot((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 63: _t->HideCaliber_Slot((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 64: _t->SmoothCurrentCurve_Slot(); break;
        case 65: _t->SaveImageAndSwcForTest_Slot(); break;
        case 66: _t->ToggleShowTraverseFlag_Slot(); break;
        case 67: _t->renderCaliber_Slot(); break;
        case 68: _t->AllCaliberBulider(); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        void **func = reinterpret_cast<void **>(_a[1]);
        {
            typedef void (NeuroGLWidget::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&NeuroGLWidget::BoxWidgetChanged_Signal)) {
                *result = 0;
            }
        }
        {
            typedef void (NeuroGLWidget::*_t)(bool );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&NeuroGLWidget::ChangeMoveCursor_Signal)) {
                *result = 1;
            }
        }
        {
            typedef void (NeuroGLWidget::*_t)(int );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&NeuroGLWidget::clippingDistanceChanged_Signal)) {
                *result = 2;
            }
        }
        {
            typedef void (NeuroGLWidget::*_t)(bool );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&NeuroGLWidget::InteractModeChanged_Signal)) {
                *result = 3;
            }
        }
        {
            typedef void (NeuroGLWidget::*_t)(int );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&NeuroGLWidget::Set2DProject_Signal)) {
                *result = 4;
            }
        }
        {
            typedef void (NeuroGLWidget::*_t)(bool );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&NeuroGLWidget::Set2DView_Signal)) {
                *result = 5;
            }
        }
        {
            typedef void (NeuroGLWidget::*_t)(bool );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&NeuroGLWidget::SetDraw_Signal)) {
                *result = 6;
            }
        }
        {
            typedef void (NeuroGLWidget::*_t)(int , int );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&NeuroGLWidget::SetTraversePos_Signal)) {
                *result = 7;
            }
        }
        {
            typedef void (NeuroGLWidget::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&NeuroGLWidget::NeedUpdateTreeWidget)) {
                *result = 8;
            }
        }
        {
            typedef void (NeuroGLWidget::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&NeuroGLWidget::NeedUpdateCheckWidget)) {
                *result = 9;
            }
        }
        {
            typedef void (NeuroGLWidget::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&NeuroGLWidget::NextTraverse_Signal)) {
                *result = 10;
            }
        }
    }
}

const QMetaObject NeuroGLWidget::staticMetaObject = {
    { &QVTKOpenGLWidget::staticMetaObject, qt_meta_stringdata_NeuroGLWidget.data,
      qt_meta_data_NeuroGLWidget,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *NeuroGLWidget::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *NeuroGLWidget::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_NeuroGLWidget.stringdata0))
        return static_cast<void*>(const_cast< NeuroGLWidget*>(this));
    return QVTKOpenGLWidget::qt_metacast(_clname);
}

int NeuroGLWidget::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QVTKOpenGLWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 69)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 69;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 69)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 69;
    }
    return _id;
}

// SIGNAL 0
void NeuroGLWidget::BoxWidgetChanged_Signal()
{
    QMetaObject::activate(this, &staticMetaObject, 0, Q_NULLPTR);
}

// SIGNAL 1
void NeuroGLWidget::ChangeMoveCursor_Signal(bool _t1)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}

// SIGNAL 2
void NeuroGLWidget::clippingDistanceChanged_Signal(int _t1)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 2, _a);
}

// SIGNAL 3
void NeuroGLWidget::InteractModeChanged_Signal(bool _t1)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 3, _a);
}

// SIGNAL 4
void NeuroGLWidget::Set2DProject_Signal(int _t1)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 4, _a);
}

// SIGNAL 5
void NeuroGLWidget::Set2DView_Signal(bool _t1)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 5, _a);
}

// SIGNAL 6
void NeuroGLWidget::SetDraw_Signal(bool _t1)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 6, _a);
}

// SIGNAL 7
void NeuroGLWidget::SetTraversePos_Signal(int _t1, int _t2)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 7, _a);
}

// SIGNAL 8
void NeuroGLWidget::NeedUpdateTreeWidget()
{
    QMetaObject::activate(this, &staticMetaObject, 8, Q_NULLPTR);
}

// SIGNAL 9
void NeuroGLWidget::NeedUpdateCheckWidget()
{
    QMetaObject::activate(this, &staticMetaObject, 9, Q_NULLPTR);
}

// SIGNAL 10
void NeuroGLWidget::NextTraverse_Signal()
{
    QMetaObject::activate(this, &staticMetaObject, 10, Q_NULLPTR);
}
QT_END_MOC_NAMESPACE
