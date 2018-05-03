#-------------------------------------------------
#
# Project created by QtCreator 2013-12-04T13:06:10
#
#-------------------------------------------------

QT       += core gui opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = SparseTracer
TEMPLATE = app

CONFIG(debug, debug|release) {
    TARGET = SparseTracer_d
    CONFIG(debug, debug|release)�� -= release
    CONFIG += debug console
    MOC_DIR += ./GeneratedFiles/debug
    OBJECTS_DIR += debug
    DESTDIR = ../x64/Debug
} else {
    TARGET = NeuroTracer
    CONFIG -= debug
    CONFIG += release console
    #QMAKE_CXXFLAGS_RELEASE += -O3
    MOC_DIR += ./GeneratedFiles/release
    OBJECTS_DIR += release
    DESTDIR = ../x64/Release
}

DEFINES +=  QT_DLL QT_OPENGL_LIB
QMAKE_CXXFLAGS += -fopenmp -std=c++0x
INCLUDEPATH += . \
               $$PWD/../../../../usr/lib/ \
              /home/zhouhang/eigen/eigen3.2.2/
LIBS += -lpthread \
        -lgomp \
        -ltiff \
        -L$$PWD/../../../../usr/lib/x86_64-linux-gnu -lGLU
DEPENDPATH += . \
                $$PWD/../../../../usr/lib/

UI_DIR += ./GeneratedFiles
RCC_DIR += ./GeneratedFiles

SOURCES += main.cpp\
        sparsetracer.cpp \
    Dialog/openimagedialog.cpp \
    DockWidget/binarysettingsdockwidget.cpp \
    Function/volumealgo.cpp \
    Function/contourutil.cpp \
    Function/binaryfilter.cpp \
    Function/IO/somawriter.cpp \
    Function/IO/somareader.cpp \
    Function/IO/imagereader.cpp \
    ngtypes/volume.cpp \
    ngtypes/tree.cpp \
    ngtypes/soma.cpp \
    ngtypes/shape.cpp \
    Function/Trace/traceutil.cpp \
    Function/IO/treewriter.cpp \
    Function/Trace/SparseTrace/sparsetracefilter.cpp \
    Render/neuroglwidget.cpp \
    Thread/binarythread.cpp \
    Render/nginteractorstyletrackballcamera.cpp \
    Function/Trace/SparseTrace/optracefilter.cpp \
    Function/IO/imagewriter.cpp \
    Function/Trace/SparseTrace/largesparsetracefilter.cpp

HEADERS  += sparsetracer.h \
    Dialog/openimagedialog.h \
    DockWidget/binarysettingsdockwidget.h \
    Function/volumealgo.h \
    Function/ineuronprocessobject.h \
    Function/ineuronio.h \
    Function/contourutil.h \
    Function/binaryfilter.h \
    Function/IO/somawriter.h \
    Function/IO/somareader.h \
    Function/IO/imagereader.h \
    ngtypes/volumecreator.h \
    ngtypes/volume.h \
    ngtypes/tree.h \
    ngtypes/soma.h \
    ngtypes/shape.h \
    ngtypes/ineurondataobject.h \
    ngtypes/basetypes.h \
    Function/Trace/traceutil.h \
    Function/IO/treewriter.h \
    Function/Trace/SparseTrace/sparsetracefilter.h \
    Render/neuroglwidget.h \
    Thread/binarythread.h \
    Render/nginteractorstyletrackballcamera.h \
    Function/Trace/SparseTrace/optracefilter.h \
    Function/IO/imagewriter.h \
    Function/Trace/SparseTrace/largesparsetracefilter.h

FORMS    += sparsetracer.ui \
    Dialog/openimagedialog.ui \
    DockWidget/binarysettingsdockwidget.ui

INCLUDEPATH += /usr/local/include/vtk-5.8
LIBS += -L/usr/local/lib/vtk-5.8 -lvtkCommon -lvtksys -lQVTK -lvtkQtChart \
 -lvtkViews -lvtkWidgets -lvtkInfovis -lvtkRendering -lvtkGraphics \
 -lvtkImaging -lvtkIO -lvtkFiltering -lvtklibxml2 -lvtkDICOMParser \
 -lvtkpng -lvtkpng -lvtktiff -lvtkzlib -lvtkjpeg -lvtkalglib -lvtkexpat \
 -lvtkverdict -lvtkmetaio -lvtkNetCDF -lvtksqlite -lvtkexoIIc -lvtkftgl \
 -lvtkfreetype -lvtkHybrid -lvtkNetCDF_cxx -lvtkVolumeRendering

RESOURCES += \
    sparsetracer.qrc
