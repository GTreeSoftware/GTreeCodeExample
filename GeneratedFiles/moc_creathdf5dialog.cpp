/****************************************************************************
** Meta object code from reading C++ file 'creathdf5dialog.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.5.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../Dialog/creathdf5dialog.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'creathdf5dialog.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.5.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_creathdf5dialog_t {
    QByteArrayData data[8];
    char stringdata0[138];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_creathdf5dialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_creathdf5dialog_t qt_meta_stringdata_creathdf5dialog = {
    {
QT_MOC_LITERAL(0, 0, 15), // "creathdf5dialog"
QT_MOC_LITERAL(1, 16, 19), // "on_okButton_clicked"
QT_MOC_LITERAL(2, 36, 0), // ""
QT_MOC_LITERAL(3, 37, 23), // "on_SelectButton_clicked"
QT_MOC_LITERAL(4, 61, 21), // "on_pathButton_clicked"
QT_MOC_LITERAL(5, 83, 23), // "on_cancelButton_clicked"
QT_MOC_LITERAL(6, 107, 15), // "on_enterPressed"
QT_MOC_LITERAL(7, 123, 14) // "getMaxLevelNum"

    },
    "creathdf5dialog\0on_okButton_clicked\0"
    "\0on_SelectButton_clicked\0on_pathButton_clicked\0"
    "on_cancelButton_clicked\0on_enterPressed\0"
    "getMaxLevelNum"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_creathdf5dialog[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       6,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   44,    2, 0x08 /* Private */,
       3,    0,   45,    2, 0x08 /* Private */,
       4,    0,   46,    2, 0x08 /* Private */,
       5,    0,   47,    2, 0x08 /* Private */,
       6,    0,   48,    2, 0x08 /* Private */,
       7,    3,   49,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Int, QMetaType::Int, QMetaType::Int, QMetaType::Int,    2,    2,    2,

       0        // eod
};

void creathdf5dialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        creathdf5dialog *_t = static_cast<creathdf5dialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->on_okButton_clicked(); break;
        case 1: _t->on_SelectButton_clicked(); break;
        case 2: _t->on_pathButton_clicked(); break;
        case 3: _t->on_cancelButton_clicked(); break;
        case 4: _t->on_enterPressed(); break;
        case 5: { int _r = _t->getMaxLevelNum((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2])),(*reinterpret_cast< int(*)>(_a[3])));
            if (_a[0]) *reinterpret_cast< int*>(_a[0]) = _r; }  break;
        default: ;
        }
    }
}

const QMetaObject creathdf5dialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_creathdf5dialog.data,
      qt_meta_data_creathdf5dialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *creathdf5dialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *creathdf5dialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_creathdf5dialog.stringdata0))
        return static_cast<void*>(const_cast< creathdf5dialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int creathdf5dialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 6)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 6;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 6)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 6;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
