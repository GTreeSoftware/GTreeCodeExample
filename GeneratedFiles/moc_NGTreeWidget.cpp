/****************************************************************************
** Meta object code from reading C++ file 'NGTreeWidget.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.5.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../DockWidget/NGTreeWidget.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'NGTreeWidget.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.5.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_NGTreeWidget_t {
    QByteArrayData data[18];
    char stringdata0[359];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_NGTreeWidget_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_NGTreeWidget_t qt_meta_stringdata_NGTreeWidget = {
    {
QT_MOC_LITERAL(0, 0, 12), // "NGTreeWidget"
QT_MOC_LITERAL(1, 13, 17), // "ActiveTree_Signal"
QT_MOC_LITERAL(2, 31, 0), // ""
QT_MOC_LITERAL(3, 32, 18), // "CompareTree_Signal"
QT_MOC_LITERAL(4, 51, 24), // "ToggleTreeVisible_Signal"
QT_MOC_LITERAL(5, 76, 15), // "GotoDiff_Signal"
QT_MOC_LITERAL(6, 92, 26), // "ToggleShowDiffArrow_Signal"
QT_MOC_LITERAL(7, 119, 20), // "ResetTraverse_Signal"
QT_MOC_LITERAL(8, 140, 29), // "ToggleShowTraverseFlag_Signal"
QT_MOC_LITERAL(9, 170, 17), // "ActivateTree_Slot"
QT_MOC_LITERAL(10, 188, 16), // "CompareTree_Slot"
QT_MOC_LITERAL(11, 205, 22), // "ToggleTreeVisible_Slot"
QT_MOC_LITERAL(12, 228, 13), // "GotoDiff_Slot"
QT_MOC_LITERAL(13, 242, 24), // "TriggeredCheckState_Slot"
QT_MOC_LITERAL(14, 267, 19), // "DeleteAnnotate_Slot"
QT_MOC_LITERAL(15, 287, 24), // "ToggleShowDiffArrow_Slot"
QT_MOC_LITERAL(16, 312, 18), // "ResetTraverse_Slot"
QT_MOC_LITERAL(17, 331, 27) // "ToggleShowTraverseFlag_Slot"

    },
    "NGTreeWidget\0ActiveTree_Signal\0\0"
    "CompareTree_Signal\0ToggleTreeVisible_Signal\0"
    "GotoDiff_Signal\0ToggleShowDiffArrow_Signal\0"
    "ResetTraverse_Signal\0ToggleShowTraverseFlag_Signal\0"
    "ActivateTree_Slot\0CompareTree_Slot\0"
    "ToggleTreeVisible_Slot\0GotoDiff_Slot\0"
    "TriggeredCheckState_Slot\0DeleteAnnotate_Slot\0"
    "ToggleShowDiffArrow_Slot\0ResetTraverse_Slot\0"
    "ToggleShowTraverseFlag_Slot"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_NGTreeWidget[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
      16,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       7,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    1,   94,    2, 0x06 /* Public */,
       3,    1,   97,    2, 0x06 /* Public */,
       4,    1,  100,    2, 0x06 /* Public */,
       5,    1,  103,    2, 0x06 /* Public */,
       6,    1,  106,    2, 0x06 /* Public */,
       7,    0,  109,    2, 0x06 /* Public */,
       8,    0,  110,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
       9,    0,  111,    2, 0x0a /* Public */,
      10,    0,  112,    2, 0x0a /* Public */,
      11,    0,  113,    2, 0x0a /* Public */,
      12,    0,  114,    2, 0x0a /* Public */,
      13,    0,  115,    2, 0x0a /* Public */,
      14,    0,  116,    2, 0x0a /* Public */,
      15,    1,  117,    2, 0x0a /* Public */,
      16,    0,  120,    2, 0x0a /* Public */,
      17,    0,  121,    2, 0x0a /* Public */,

 // signals: parameters
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void,
    QMetaType::Void,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Bool,    2,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void NGTreeWidget::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        NGTreeWidget *_t = static_cast<NGTreeWidget *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->ActiveTree_Signal((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 1: _t->CompareTree_Signal((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 2: _t->ToggleTreeVisible_Signal((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 3: _t->GotoDiff_Signal((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 4: _t->ToggleShowDiffArrow_Signal((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 5: _t->ResetTraverse_Signal(); break;
        case 6: _t->ToggleShowTraverseFlag_Signal(); break;
        case 7: _t->ActivateTree_Slot(); break;
        case 8: _t->CompareTree_Slot(); break;
        case 9: _t->ToggleTreeVisible_Slot(); break;
        case 10: _t->GotoDiff_Slot(); break;
        case 11: _t->TriggeredCheckState_Slot(); break;
        case 12: _t->DeleteAnnotate_Slot(); break;
        case 13: _t->ToggleShowDiffArrow_Slot((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 14: _t->ResetTraverse_Slot(); break;
        case 15: _t->ToggleShowTraverseFlag_Slot(); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        void **func = reinterpret_cast<void **>(_a[1]);
        {
            typedef void (NGTreeWidget::*_t)(int );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&NGTreeWidget::ActiveTree_Signal)) {
                *result = 0;
            }
        }
        {
            typedef void (NGTreeWidget::*_t)(int );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&NGTreeWidget::CompareTree_Signal)) {
                *result = 1;
            }
        }
        {
            typedef void (NGTreeWidget::*_t)(int );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&NGTreeWidget::ToggleTreeVisible_Signal)) {
                *result = 2;
            }
        }
        {
            typedef void (NGTreeWidget::*_t)(int );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&NGTreeWidget::GotoDiff_Signal)) {
                *result = 3;
            }
        }
        {
            typedef void (NGTreeWidget::*_t)(bool );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&NGTreeWidget::ToggleShowDiffArrow_Signal)) {
                *result = 4;
            }
        }
        {
            typedef void (NGTreeWidget::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&NGTreeWidget::ResetTraverse_Signal)) {
                *result = 5;
            }
        }
        {
            typedef void (NGTreeWidget::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&NGTreeWidget::ToggleShowTraverseFlag_Signal)) {
                *result = 6;
            }
        }
    }
}

const QMetaObject NGTreeWidget::staticMetaObject = {
    { &QTreeWidget::staticMetaObject, qt_meta_stringdata_NGTreeWidget.data,
      qt_meta_data_NGTreeWidget,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *NGTreeWidget::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *NGTreeWidget::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_NGTreeWidget.stringdata0))
        return static_cast<void*>(const_cast< NGTreeWidget*>(this));
    return QTreeWidget::qt_metacast(_clname);
}

int NGTreeWidget::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QTreeWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 16)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 16;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 16)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 16;
    }
    return _id;
}

// SIGNAL 0
void NGTreeWidget::ActiveTree_Signal(int _t1)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void NGTreeWidget::CompareTree_Signal(int _t1)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}

// SIGNAL 2
void NGTreeWidget::ToggleTreeVisible_Signal(int _t1)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 2, _a);
}

// SIGNAL 3
void NGTreeWidget::GotoDiff_Signal(int _t1)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 3, _a);
}

// SIGNAL 4
void NGTreeWidget::ToggleShowDiffArrow_Signal(bool _t1)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 4, _a);
}

// SIGNAL 5
void NGTreeWidget::ResetTraverse_Signal()
{
    QMetaObject::activate(this, &staticMetaObject, 5, Q_NULLPTR);
}

// SIGNAL 6
void NGTreeWidget::ToggleShowTraverseFlag_Signal()
{
    QMetaObject::activate(this, &staticMetaObject, 6, Q_NULLPTR);
}
QT_END_MOC_NAMESPACE
