/****************************************************************************
** Meta object code from reading C++ file 'binarysettingsdockwidget.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.5.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../DockWidget/binarysettingsdockwidget.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'binarysettingsdockwidget.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.5.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_BinarySettingsDockWidget_t {
    QByteArrayData data[80];
    char stringdata0[2133];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_BinarySettingsDockWidget_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_BinarySettingsDockWidget_t qt_meta_stringdata_BinarySettingsDockWidget = {
    {
QT_MOC_LITERAL(0, 0, 24), // "BinarySettingsDockWidget"
QT_MOC_LITERAL(1, 25, 21), // "ApplyReadImage_Signal"
QT_MOC_LITERAL(2, 47, 0), // ""
QT_MOC_LITERAL(3, 48, 21), // "DisplayNewOpac_Signal"
QT_MOC_LITERAL(4, 70, 25), // "MaxBoundNumChanged_Signal"
QT_MOC_LITERAL(5, 96, 17), // "NeedBinary_Signal"
QT_MOC_LITERAL(6, 114, 16), // "NeedTrace_Signal"
QT_MOC_LITERAL(7, 131, 17), // "ROIChanged_Signal"
QT_MOC_LITERAL(8, 149, 23), // "ThicknessChanged_Signal"
QT_MOC_LITERAL(9, 173, 17), // "ClearCache_Signal"
QT_MOC_LITERAL(10, 191, 17), // "ActiveTree_Signal"
QT_MOC_LITERAL(11, 209, 18), // "CompareTree_Signal"
QT_MOC_LITERAL(12, 228, 24), // "ToggleTreeVisible_Signal"
QT_MOC_LITERAL(13, 253, 15), // "GotoDiff_Signal"
QT_MOC_LITERAL(14, 269, 21), // "ToggleTraverse_Signal"
QT_MOC_LITERAL(15, 291, 15), // "Annotate_Signal"
QT_MOC_LITERAL(16, 307, 19), // "BackTraverse_Signal"
QT_MOC_LITERAL(17, 327, 19), // "NextTraverse_Signal"
QT_MOC_LITERAL(18, 347, 26), // "ToggleShowDiffArrow_Signal"
QT_MOC_LITERAL(19, 374, 20), // "ResetTraverse_Signal"
QT_MOC_LITERAL(20, 395, 29), // "ToggleShowTraverseFlag_Signal"
QT_MOC_LITERAL(21, 425, 22), // "OpacAdjustApply_Signal"
QT_MOC_LITERAL(22, 448, 24), // "PreviewOpacAdjust_Signal"
QT_MOC_LITERAL(23, 473, 20), // "PreviewBinary_Signal"
QT_MOC_LITERAL(24, 494, 24), // "PreviewAxonBinary_Signal"
QT_MOC_LITERAL(25, 519, 30), // "on_xMinspinBox_editingFinished"
QT_MOC_LITERAL(26, 550, 30), // "on_xMaxspinBox_editingFinished"
QT_MOC_LITERAL(27, 581, 30), // "on_yMinspinBox_editingFinished"
QT_MOC_LITERAL(28, 612, 30), // "on_yMaxspinBox_editingFinished"
QT_MOC_LITERAL(29, 643, 30), // "on_zMinspinBox_editingFinished"
QT_MOC_LITERAL(30, 674, 30), // "on_zMaxspinBox_editingFinished"
QT_MOC_LITERAL(31, 705, 31), // "on_applyImageReadButton_clicked"
QT_MOC_LITERAL(32, 737, 39), // "on_axonDiffuseValueSpinBox_va..."
QT_MOC_LITERAL(33, 777, 36), // "on_axonThresholdSpinBox_value..."
QT_MOC_LITERAL(34, 814, 37), // "on_axonTraceValueSpinBox_valu..."
QT_MOC_LITERAL(35, 852, 48), // "on_axonSemiautoTraceTresholdS..."
QT_MOC_LITERAL(36, 901, 50), // "on_axonSemiautoConnectResampl..."
QT_MOC_LITERAL(37, 952, 35), // "on_diffuseValueSpinBox_valueC..."
QT_MOC_LITERAL(38, 988, 31), // "on_highOpacSpinBox_valueChanged"
QT_MOC_LITERAL(39, 1020, 30), // "on_lowOpacSpinBox_valueChanged"
QT_MOC_LITERAL(40, 1051, 37), // "on_maxBoundNumSpinBox_editing..."
QT_MOC_LITERAL(41, 1089, 32), // "on_thicknessspinBox_valueChanged"
QT_MOC_LITERAL(42, 1122, 32), // "on_thresholdSpinBox_valueChanged"
QT_MOC_LITERAL(43, 1155, 33), // "on_traceValueSpinBox_valueCha..."
QT_MOC_LITERAL(44, 1189, 16), // "OpacValueChanged"
QT_MOC_LITERAL(45, 1206, 17), // "SetThickness_Slot"
QT_MOC_LITERAL(46, 1224, 33), // "on_xBlockSizespinBox_valueCha..."
QT_MOC_LITERAL(47, 1258, 33), // "on_yBlockSizespinBox_valueCha..."
QT_MOC_LITERAL(48, 1292, 33), // "on_zBlockSizespinBox_valueCha..."
QT_MOC_LITERAL(49, 1326, 35), // "on_action_ClearCache_Button_c..."
QT_MOC_LITERAL(50, 1362, 17), // "ActivateTree_Slot"
QT_MOC_LITERAL(51, 1380, 16), // "CompareTree_Slot"
QT_MOC_LITERAL(52, 1397, 22), // "ToggleTreeVisible_Slot"
QT_MOC_LITERAL(53, 1420, 13), // "GotoDiff_Slot"
QT_MOC_LITERAL(54, 1434, 19), // "ToggleTraverse_Slot"
QT_MOC_LITERAL(55, 1454, 13), // "Annotate_Slot"
QT_MOC_LITERAL(56, 1468, 17), // "BackTraverse_Slot"
QT_MOC_LITERAL(57, 1486, 17), // "NextTraverse_Slot"
QT_MOC_LITERAL(58, 1504, 24), // "ToggleShowDiffArrow_Slot"
QT_MOC_LITERAL(59, 1529, 18), // "ResetTraverse_Slot"
QT_MOC_LITERAL(60, 1548, 28), // "on_SmoothSlider_valueChanged"
QT_MOC_LITERAL(61, 1577, 4), // "arg1"
QT_MOC_LITERAL(62, 1582, 21), // "UpdateTreeWidget_Slot"
QT_MOC_LITERAL(63, 1604, 22), // "UpdateCheckWidget_Slot"
QT_MOC_LITERAL(64, 1627, 27), // "ToggleShowTraverseFlag_Slot"
QT_MOC_LITERAL(65, 1655, 23), // "on_enableSVMBox_toggled"
QT_MOC_LITERAL(66, 1679, 25), // "on_GPSSVMcheckBox_toggled"
QT_MOC_LITERAL(67, 1705, 31), // "on_strongSignalCheckBox_toggled"
QT_MOC_LITERAL(68, 1737, 29), // "on_xScaleSpinBox_valueChanged"
QT_MOC_LITERAL(69, 1767, 29), // "on_yScaleSpinBox_valueChanged"
QT_MOC_LITERAL(70, 1797, 29), // "on_zScaleSpinBox_valueChanged"
QT_MOC_LITERAL(71, 1827, 34), // "on_origLowOpacSpinBox_valueCh..."
QT_MOC_LITERAL(72, 1862, 35), // "on_origHighOpacSpinBox_valueC..."
QT_MOC_LITERAL(73, 1898, 34), // "on_destLowOpacSpinBox_valueCh..."
QT_MOC_LITERAL(74, 1933, 35), // "on_destHighOpacSpinBox_valueC..."
QT_MOC_LITERAL(75, 1969, 30), // "on_imageOpacInfoButton_clicked"
QT_MOC_LITERAL(76, 2000, 26), // "on_opacApplyButton_clicked"
QT_MOC_LITERAL(77, 2027, 39), // "on_previewIlluminationMapButt..."
QT_MOC_LITERAL(78, 2067, 30), // "on_previewBinaryButton_clicked"
QT_MOC_LITERAL(79, 2098, 34) // "on_previewAxonBinaryButton_cl..."

    },
    "BinarySettingsDockWidget\0ApplyReadImage_Signal\0"
    "\0DisplayNewOpac_Signal\0MaxBoundNumChanged_Signal\0"
    "NeedBinary_Signal\0NeedTrace_Signal\0"
    "ROIChanged_Signal\0ThicknessChanged_Signal\0"
    "ClearCache_Signal\0ActiveTree_Signal\0"
    "CompareTree_Signal\0ToggleTreeVisible_Signal\0"
    "GotoDiff_Signal\0ToggleTraverse_Signal\0"
    "Annotate_Signal\0BackTraverse_Signal\0"
    "NextTraverse_Signal\0ToggleShowDiffArrow_Signal\0"
    "ResetTraverse_Signal\0ToggleShowTraverseFlag_Signal\0"
    "OpacAdjustApply_Signal\0PreviewOpacAdjust_Signal\0"
    "PreviewBinary_Signal\0PreviewAxonBinary_Signal\0"
    "on_xMinspinBox_editingFinished\0"
    "on_xMaxspinBox_editingFinished\0"
    "on_yMinspinBox_editingFinished\0"
    "on_yMaxspinBox_editingFinished\0"
    "on_zMinspinBox_editingFinished\0"
    "on_zMaxspinBox_editingFinished\0"
    "on_applyImageReadButton_clicked\0"
    "on_axonDiffuseValueSpinBox_valueChanged\0"
    "on_axonThresholdSpinBox_valueChanged\0"
    "on_axonTraceValueSpinBox_valueChanged\0"
    "on_axonSemiautoTraceTresholdSpinBox_valueChanged\0"
    "on_axonSemiautoConnectResampleSpinBox_valueChanged\0"
    "on_diffuseValueSpinBox_valueChanged\0"
    "on_highOpacSpinBox_valueChanged\0"
    "on_lowOpacSpinBox_valueChanged\0"
    "on_maxBoundNumSpinBox_editingFinished\0"
    "on_thicknessspinBox_valueChanged\0"
    "on_thresholdSpinBox_valueChanged\0"
    "on_traceValueSpinBox_valueChanged\0"
    "OpacValueChanged\0SetThickness_Slot\0"
    "on_xBlockSizespinBox_valueChanged\0"
    "on_yBlockSizespinBox_valueChanged\0"
    "on_zBlockSizespinBox_valueChanged\0"
    "on_action_ClearCache_Button_clicked\0"
    "ActivateTree_Slot\0CompareTree_Slot\0"
    "ToggleTreeVisible_Slot\0GotoDiff_Slot\0"
    "ToggleTraverse_Slot\0Annotate_Slot\0"
    "BackTraverse_Slot\0NextTraverse_Slot\0"
    "ToggleShowDiffArrow_Slot\0ResetTraverse_Slot\0"
    "on_SmoothSlider_valueChanged\0arg1\0"
    "UpdateTreeWidget_Slot\0UpdateCheckWidget_Slot\0"
    "ToggleShowTraverseFlag_Slot\0"
    "on_enableSVMBox_toggled\0"
    "on_GPSSVMcheckBox_toggled\0"
    "on_strongSignalCheckBox_toggled\0"
    "on_xScaleSpinBox_valueChanged\0"
    "on_yScaleSpinBox_valueChanged\0"
    "on_zScaleSpinBox_valueChanged\0"
    "on_origLowOpacSpinBox_valueChanged\0"
    "on_origHighOpacSpinBox_valueChanged\0"
    "on_destLowOpacSpinBox_valueChanged\0"
    "on_destHighOpacSpinBox_valueChanged\0"
    "on_imageOpacInfoButton_clicked\0"
    "on_opacApplyButton_clicked\0"
    "on_previewIlluminationMapButton_clicked\0"
    "on_previewBinaryButton_clicked\0"
    "on_previewAxonBinaryButton_clicked"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_BinarySettingsDockWidget[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
      77,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
      23,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    0,  399,    2, 0x06 /* Public */,
       3,    0,  400,    2, 0x06 /* Public */,
       4,    1,  401,    2, 0x06 /* Public */,
       5,    0,  404,    2, 0x06 /* Public */,
       6,    0,  405,    2, 0x06 /* Public */,
       7,    0,  406,    2, 0x06 /* Public */,
       8,    1,  407,    2, 0x06 /* Public */,
       9,    0,  410,    2, 0x06 /* Public */,
      10,    1,  411,    2, 0x06 /* Public */,
      11,    1,  414,    2, 0x06 /* Public */,
      12,    1,  417,    2, 0x06 /* Public */,
      13,    1,  420,    2, 0x06 /* Public */,
      14,    1,  423,    2, 0x06 /* Public */,
      15,    0,  426,    2, 0x06 /* Public */,
      16,    0,  427,    2, 0x06 /* Public */,
      17,    0,  428,    2, 0x06 /* Public */,
      18,    1,  429,    2, 0x06 /* Public */,
      19,    0,  432,    2, 0x06 /* Public */,
      20,    0,  433,    2, 0x06 /* Public */,
      21,    0,  434,    2, 0x06 /* Public */,
      22,    1,  435,    2, 0x06 /* Public */,
      23,    1,  438,    2, 0x06 /* Public */,
      24,    1,  441,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
      25,    0,  444,    2, 0x08 /* Private */,
      26,    0,  445,    2, 0x08 /* Private */,
      27,    0,  446,    2, 0x08 /* Private */,
      28,    0,  447,    2, 0x08 /* Private */,
      29,    0,  448,    2, 0x08 /* Private */,
      30,    0,  449,    2, 0x08 /* Private */,
      31,    0,  450,    2, 0x08 /* Private */,
      32,    1,  451,    2, 0x08 /* Private */,
      33,    1,  454,    2, 0x08 /* Private */,
      34,    1,  457,    2, 0x08 /* Private */,
      35,    1,  460,    2, 0x08 /* Private */,
      36,    1,  463,    2, 0x08 /* Private */,
      37,    1,  466,    2, 0x08 /* Private */,
      38,    1,  469,    2, 0x08 /* Private */,
      39,    1,  472,    2, 0x08 /* Private */,
      40,    0,  475,    2, 0x08 /* Private */,
      41,    1,  476,    2, 0x08 /* Private */,
      42,    1,  479,    2, 0x08 /* Private */,
      43,    1,  482,    2, 0x08 /* Private */,
      44,    1,  485,    2, 0x08 /* Private */,
      45,    1,  488,    2, 0x08 /* Private */,
      46,    1,  491,    2, 0x08 /* Private */,
      47,    1,  494,    2, 0x08 /* Private */,
      48,    1,  497,    2, 0x08 /* Private */,
      49,    0,  500,    2, 0x08 /* Private */,
      50,    1,  501,    2, 0x08 /* Private */,
      51,    1,  504,    2, 0x08 /* Private */,
      52,    1,  507,    2, 0x08 /* Private */,
      53,    1,  510,    2, 0x08 /* Private */,
      54,    1,  513,    2, 0x08 /* Private */,
      55,    0,  516,    2, 0x08 /* Private */,
      56,    0,  517,    2, 0x08 /* Private */,
      57,    0,  518,    2, 0x08 /* Private */,
      58,    1,  519,    2, 0x08 /* Private */,
      59,    0,  522,    2, 0x08 /* Private */,
      60,    1,  523,    2, 0x08 /* Private */,
      62,    0,  526,    2, 0x08 /* Private */,
      63,    0,  527,    2, 0x08 /* Private */,
      64,    0,  528,    2, 0x08 /* Private */,
      65,    1,  529,    2, 0x08 /* Private */,
      66,    1,  532,    2, 0x08 /* Private */,
      67,    1,  535,    2, 0x08 /* Private */,
      68,    1,  538,    2, 0x08 /* Private */,
      69,    1,  541,    2, 0x08 /* Private */,
      70,    1,  544,    2, 0x08 /* Private */,
      71,    1,  547,    2, 0x08 /* Private */,
      72,    1,  550,    2, 0x08 /* Private */,
      73,    1,  553,    2, 0x08 /* Private */,
      74,    1,  556,    2, 0x08 /* Private */,
      75,    0,  559,    2, 0x08 /* Private */,
      76,    0,  560,    2, 0x08 /* Private */,
      77,    1,  561,    2, 0x08 /* Private */,
      78,    1,  564,    2, 0x08 /* Private */,
      79,    1,  567,    2, 0x08 /* Private */,

 // signals: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void, QMetaType::Bool,    2,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Double,    2,
    QMetaType::Void, QMetaType::Double,    2,
    QMetaType::Void, QMetaType::Double,    2,
    QMetaType::Void, QMetaType::Double,    2,
    QMetaType::Void, QMetaType::Double,    2,
    QMetaType::Void, QMetaType::Double,    2,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Double,    2,
    QMetaType::Void, QMetaType::Double,    2,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int,   61,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void, QMetaType::Bool,    2,

       0        // eod
};

void BinarySettingsDockWidget::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        BinarySettingsDockWidget *_t = static_cast<BinarySettingsDockWidget *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->ApplyReadImage_Signal(); break;
        case 1: _t->DisplayNewOpac_Signal(); break;
        case 2: _t->MaxBoundNumChanged_Signal((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 3: _t->NeedBinary_Signal(); break;
        case 4: _t->NeedTrace_Signal(); break;
        case 5: _t->ROIChanged_Signal(); break;
        case 6: _t->ThicknessChanged_Signal((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 7: _t->ClearCache_Signal(); break;
        case 8: _t->ActiveTree_Signal((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 9: _t->CompareTree_Signal((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 10: _t->ToggleTreeVisible_Signal((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 11: _t->GotoDiff_Signal((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 12: _t->ToggleTraverse_Signal((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 13: _t->Annotate_Signal(); break;
        case 14: _t->BackTraverse_Signal(); break;
        case 15: _t->NextTraverse_Signal(); break;
        case 16: _t->ToggleShowDiffArrow_Signal((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 17: _t->ResetTraverse_Signal(); break;
        case 18: _t->ToggleShowTraverseFlag_Signal(); break;
        case 19: _t->OpacAdjustApply_Signal(); break;
        case 20: _t->PreviewOpacAdjust_Signal((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 21: _t->PreviewBinary_Signal((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 22: _t->PreviewAxonBinary_Signal((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 23: _t->on_xMinspinBox_editingFinished(); break;
        case 24: _t->on_xMaxspinBox_editingFinished(); break;
        case 25: _t->on_yMinspinBox_editingFinished(); break;
        case 26: _t->on_yMaxspinBox_editingFinished(); break;
        case 27: _t->on_zMinspinBox_editingFinished(); break;
        case 28: _t->on_zMaxspinBox_editingFinished(); break;
        case 29: _t->on_applyImageReadButton_clicked(); break;
        case 30: _t->on_axonDiffuseValueSpinBox_valueChanged((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 31: _t->on_axonThresholdSpinBox_valueChanged((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 32: _t->on_axonTraceValueSpinBox_valueChanged((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 33: _t->on_axonSemiautoTraceTresholdSpinBox_valueChanged((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 34: _t->on_axonSemiautoConnectResampleSpinBox_valueChanged((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 35: _t->on_diffuseValueSpinBox_valueChanged((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 36: _t->on_highOpacSpinBox_valueChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 37: _t->on_lowOpacSpinBox_valueChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 38: _t->on_maxBoundNumSpinBox_editingFinished(); break;
        case 39: _t->on_thicknessspinBox_valueChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 40: _t->on_thresholdSpinBox_valueChanged((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 41: _t->on_traceValueSpinBox_valueChanged((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 42: _t->OpacValueChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 43: _t->SetThickness_Slot((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 44: _t->on_xBlockSizespinBox_valueChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 45: _t->on_yBlockSizespinBox_valueChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 46: _t->on_zBlockSizespinBox_valueChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 47: _t->on_action_ClearCache_Button_clicked(); break;
        case 48: _t->ActivateTree_Slot((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 49: _t->CompareTree_Slot((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 50: _t->ToggleTreeVisible_Slot((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 51: _t->GotoDiff_Slot((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 52: _t->ToggleTraverse_Slot((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 53: _t->Annotate_Slot(); break;
        case 54: _t->BackTraverse_Slot(); break;
        case 55: _t->NextTraverse_Slot(); break;
        case 56: _t->ToggleShowDiffArrow_Slot((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 57: _t->ResetTraverse_Slot(); break;
        case 58: _t->on_SmoothSlider_valueChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 59: _t->UpdateTreeWidget_Slot(); break;
        case 60: _t->UpdateCheckWidget_Slot(); break;
        case 61: _t->ToggleShowTraverseFlag_Slot(); break;
        case 62: _t->on_enableSVMBox_toggled((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 63: _t->on_GPSSVMcheckBox_toggled((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 64: _t->on_strongSignalCheckBox_toggled((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 65: _t->on_xScaleSpinBox_valueChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 66: _t->on_yScaleSpinBox_valueChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 67: _t->on_zScaleSpinBox_valueChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 68: _t->on_origLowOpacSpinBox_valueChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 69: _t->on_origHighOpacSpinBox_valueChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 70: _t->on_destLowOpacSpinBox_valueChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 71: _t->on_destHighOpacSpinBox_valueChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 72: _t->on_imageOpacInfoButton_clicked(); break;
        case 73: _t->on_opacApplyButton_clicked(); break;
        case 74: _t->on_previewIlluminationMapButton_clicked((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 75: _t->on_previewBinaryButton_clicked((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 76: _t->on_previewAxonBinaryButton_clicked((*reinterpret_cast< bool(*)>(_a[1]))); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        void **func = reinterpret_cast<void **>(_a[1]);
        {
            typedef void (BinarySettingsDockWidget::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&BinarySettingsDockWidget::ApplyReadImage_Signal)) {
                *result = 0;
            }
        }
        {
            typedef void (BinarySettingsDockWidget::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&BinarySettingsDockWidget::DisplayNewOpac_Signal)) {
                *result = 1;
            }
        }
        {
            typedef void (BinarySettingsDockWidget::*_t)(int );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&BinarySettingsDockWidget::MaxBoundNumChanged_Signal)) {
                *result = 2;
            }
        }
        {
            typedef void (BinarySettingsDockWidget::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&BinarySettingsDockWidget::NeedBinary_Signal)) {
                *result = 3;
            }
        }
        {
            typedef void (BinarySettingsDockWidget::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&BinarySettingsDockWidget::NeedTrace_Signal)) {
                *result = 4;
            }
        }
        {
            typedef void (BinarySettingsDockWidget::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&BinarySettingsDockWidget::ROIChanged_Signal)) {
                *result = 5;
            }
        }
        {
            typedef void (BinarySettingsDockWidget::*_t)(int );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&BinarySettingsDockWidget::ThicknessChanged_Signal)) {
                *result = 6;
            }
        }
        {
            typedef void (BinarySettingsDockWidget::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&BinarySettingsDockWidget::ClearCache_Signal)) {
                *result = 7;
            }
        }
        {
            typedef void (BinarySettingsDockWidget::*_t)(int );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&BinarySettingsDockWidget::ActiveTree_Signal)) {
                *result = 8;
            }
        }
        {
            typedef void (BinarySettingsDockWidget::*_t)(int );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&BinarySettingsDockWidget::CompareTree_Signal)) {
                *result = 9;
            }
        }
        {
            typedef void (BinarySettingsDockWidget::*_t)(int );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&BinarySettingsDockWidget::ToggleTreeVisible_Signal)) {
                *result = 10;
            }
        }
        {
            typedef void (BinarySettingsDockWidget::*_t)(int );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&BinarySettingsDockWidget::GotoDiff_Signal)) {
                *result = 11;
            }
        }
        {
            typedef void (BinarySettingsDockWidget::*_t)(bool );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&BinarySettingsDockWidget::ToggleTraverse_Signal)) {
                *result = 12;
            }
        }
        {
            typedef void (BinarySettingsDockWidget::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&BinarySettingsDockWidget::Annotate_Signal)) {
                *result = 13;
            }
        }
        {
            typedef void (BinarySettingsDockWidget::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&BinarySettingsDockWidget::BackTraverse_Signal)) {
                *result = 14;
            }
        }
        {
            typedef void (BinarySettingsDockWidget::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&BinarySettingsDockWidget::NextTraverse_Signal)) {
                *result = 15;
            }
        }
        {
            typedef void (BinarySettingsDockWidget::*_t)(bool );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&BinarySettingsDockWidget::ToggleShowDiffArrow_Signal)) {
                *result = 16;
            }
        }
        {
            typedef void (BinarySettingsDockWidget::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&BinarySettingsDockWidget::ResetTraverse_Signal)) {
                *result = 17;
            }
        }
        {
            typedef void (BinarySettingsDockWidget::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&BinarySettingsDockWidget::ToggleShowTraverseFlag_Signal)) {
                *result = 18;
            }
        }
        {
            typedef void (BinarySettingsDockWidget::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&BinarySettingsDockWidget::OpacAdjustApply_Signal)) {
                *result = 19;
            }
        }
        {
            typedef void (BinarySettingsDockWidget::*_t)(bool );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&BinarySettingsDockWidget::PreviewOpacAdjust_Signal)) {
                *result = 20;
            }
        }
        {
            typedef void (BinarySettingsDockWidget::*_t)(bool );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&BinarySettingsDockWidget::PreviewBinary_Signal)) {
                *result = 21;
            }
        }
        {
            typedef void (BinarySettingsDockWidget::*_t)(bool );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&BinarySettingsDockWidget::PreviewAxonBinary_Signal)) {
                *result = 22;
            }
        }
    }
}

const QMetaObject BinarySettingsDockWidget::staticMetaObject = {
    { &QDockWidget::staticMetaObject, qt_meta_stringdata_BinarySettingsDockWidget.data,
      qt_meta_data_BinarySettingsDockWidget,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *BinarySettingsDockWidget::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *BinarySettingsDockWidget::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_BinarySettingsDockWidget.stringdata0))
        return static_cast<void*>(const_cast< BinarySettingsDockWidget*>(this));
    return QDockWidget::qt_metacast(_clname);
}

int BinarySettingsDockWidget::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDockWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 77)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 77;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 77)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 77;
    }
    return _id;
}

// SIGNAL 0
void BinarySettingsDockWidget::ApplyReadImage_Signal()
{
    QMetaObject::activate(this, &staticMetaObject, 0, Q_NULLPTR);
}

// SIGNAL 1
void BinarySettingsDockWidget::DisplayNewOpac_Signal()
{
    QMetaObject::activate(this, &staticMetaObject, 1, Q_NULLPTR);
}

// SIGNAL 2
void BinarySettingsDockWidget::MaxBoundNumChanged_Signal(int _t1)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 2, _a);
}

// SIGNAL 3
void BinarySettingsDockWidget::NeedBinary_Signal()
{
    QMetaObject::activate(this, &staticMetaObject, 3, Q_NULLPTR);
}

// SIGNAL 4
void BinarySettingsDockWidget::NeedTrace_Signal()
{
    QMetaObject::activate(this, &staticMetaObject, 4, Q_NULLPTR);
}

// SIGNAL 5
void BinarySettingsDockWidget::ROIChanged_Signal()
{
    QMetaObject::activate(this, &staticMetaObject, 5, Q_NULLPTR);
}

// SIGNAL 6
void BinarySettingsDockWidget::ThicknessChanged_Signal(int _t1)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 6, _a);
}

// SIGNAL 7
void BinarySettingsDockWidget::ClearCache_Signal()
{
    QMetaObject::activate(this, &staticMetaObject, 7, Q_NULLPTR);
}

// SIGNAL 8
void BinarySettingsDockWidget::ActiveTree_Signal(int _t1)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 8, _a);
}

// SIGNAL 9
void BinarySettingsDockWidget::CompareTree_Signal(int _t1)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 9, _a);
}

// SIGNAL 10
void BinarySettingsDockWidget::ToggleTreeVisible_Signal(int _t1)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 10, _a);
}

// SIGNAL 11
void BinarySettingsDockWidget::GotoDiff_Signal(int _t1)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 11, _a);
}

// SIGNAL 12
void BinarySettingsDockWidget::ToggleTraverse_Signal(bool _t1)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 12, _a);
}

// SIGNAL 13
void BinarySettingsDockWidget::Annotate_Signal()
{
    QMetaObject::activate(this, &staticMetaObject, 13, Q_NULLPTR);
}

// SIGNAL 14
void BinarySettingsDockWidget::BackTraverse_Signal()
{
    QMetaObject::activate(this, &staticMetaObject, 14, Q_NULLPTR);
}

// SIGNAL 15
void BinarySettingsDockWidget::NextTraverse_Signal()
{
    QMetaObject::activate(this, &staticMetaObject, 15, Q_NULLPTR);
}

// SIGNAL 16
void BinarySettingsDockWidget::ToggleShowDiffArrow_Signal(bool _t1)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 16, _a);
}

// SIGNAL 17
void BinarySettingsDockWidget::ResetTraverse_Signal()
{
    QMetaObject::activate(this, &staticMetaObject, 17, Q_NULLPTR);
}

// SIGNAL 18
void BinarySettingsDockWidget::ToggleShowTraverseFlag_Signal()
{
    QMetaObject::activate(this, &staticMetaObject, 18, Q_NULLPTR);
}

// SIGNAL 19
void BinarySettingsDockWidget::OpacAdjustApply_Signal()
{
    QMetaObject::activate(this, &staticMetaObject, 19, Q_NULLPTR);
}

// SIGNAL 20
void BinarySettingsDockWidget::PreviewOpacAdjust_Signal(bool _t1)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 20, _a);
}

// SIGNAL 21
void BinarySettingsDockWidget::PreviewBinary_Signal(bool _t1)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 21, _a);
}

// SIGNAL 22
void BinarySettingsDockWidget::PreviewAxonBinary_Signal(bool _t1)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 22, _a);
}
QT_END_MOC_NAMESPACE
