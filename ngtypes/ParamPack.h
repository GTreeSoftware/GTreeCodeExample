#ifndef PARAMPACK_H//#pragma once
#define PARAMPACK_H
#include <QString>
#include <QAction>
#include <vector>
#include <deque>
#include "ngtypes/basetypes.h"
#include "ngtypes/soma.h"
#include "ngtypes/ineurondataobject.h"
#include <tuple>
#include <omp.h>

class NeuronTree;
typedef std::shared_ptr<NeuronTree> NeuronTreePointer;
class NeuronPopulation;
typedef std::shared_ptr<NeuronPopulation> NeuronPopulationPointer;
typedef std::shared_ptr<const NeuronPopulation> CNeuronPopulationPointer;

struct ParamPack 
{
	ParamPack() :isLocalTrace_(false), enableSVM_(false), strongSignalMode_(true), curTraceTopoID_(0), dataMode_(READMODE::INVALID), filterNum_(5), threadNum_(omp_get_max_threads()), xScale_(1), yScale_(1), zScale_(1),
		maxBoundNum_(10), thicknesss_(10.0), lowOpac_(0), highOpac_(1000), lowAdjustOpac_(0), highAdjustOpac_(0), lowDestOpac_(0), highDestOpac_(0),
		xExtract_(60), yExtract_(60), zExtract_(54),
		xMin_(0), yMin_(0), zMin_(0), xMax_(0), yMax_(0), zMax_(0), xRangeMax_(0), yRangeMax_(0), zRangeMax_(0),
		binThreshold_(20.0), axonBinaryThreshold_(20.0),
		xRes_(1.0), yRes_(1.0), zRes_(1.0), diffuseValue_(4096.0), axonDiffuseValue_(50.0), axonSemiautoTraceThreshold(0.08), axonSemiautoConnectResampleThreshold(4),
    boundaryDistanceThreshold_(4.0), traceValue_(1.0), axonTraceValue_(1.0), bifurcationValue_(4.0), crossDistanceThrev_(10.0),
    mostdLevel_(1), maxValue_(4096), b(0.0), runningMinutes_(0),smooth_level(0.01){
        NG_SMART_POINTER_NEW_DEFAULT(Soma, SomaList);
        w = VecXd(1); w.setZero();
        defaultDir = QString("./");
    }

    //
    bool isLocalTrace_;
    bool enableSVM_;
    bool strongSignalMode_;
    //
    size_t curTraceTopoID_;
    //
    READMODE dataMode_;
    //binary
    int filterNum_;
    int threadNum_; 
    //image scale
    int xScale_;
    int yScale_;
    int zScale_;
    //bound check
    int maxBoundNum_;
    //display
    int thicknesss_;
    int lowOpac_;
    int highOpac_;
    //adjust illumination
    int lowAdjustOpac_;
    int highAdjustOpac_;
    int lowDestOpac_;
    int highDestOpac_;
    //image extract range
    int xExtract_;
    int yExtract_;
    int zExtract_;
    //image range
    int xMin_;
    int yMin_;
    int zMin_;
    int xMax_;
    int yMax_;
    int zMax_;
    

    int xRangeMax_;
    int yRangeMax_;
    int zRangeMax_;
    double xRes_;   //Resolution
    double yRes_;
    double zRes_;
    //binary threshold
    double binThreshold_;
    double axonBinaryThreshold_;
    //trace
    double traceValue_;
    double axonTraceValue_;
    double bifurcationValue_;
    double diffuseValue_;
    double axonDiffuseValue_;
    double boundaryDistanceThreshold_;
    double crossDistanceThrev_;
	double axonSemiautoTraceThreshold;//2018-1-22
	double axonSemiautoConnectResampleThreshold;//2018-2-9

	// smooth option 
	double kernel_size;
	double smooth_level;

    //DATA
    //DATATYPE dataType_;
    int mostdLevel_;
    int maxValue_;
    IDataPointer OrigImage;
    IDataPointer BinImage;
    IDataPointer BackImage;
    IDataPointer SomaList;
    IDataPointer LayerImage;
    IDataPointer previewImage;
    IDataPointer separateTree;//true tree data
    QString defaultDir;

    //SVM
    VecXd w; double b;

    //end point of current trace curve
    VectorVec3d allowableBoundCheckPtList_;
    //manual
    VectorVec5d manualLabelCurve_;
    //active
    NeuronTreePointer activeTree;

	//record 
	std::vector<Record> UserRecord;

    //time
    size_t runningMinutes_;

	//CaliberBulider data
	std::shared_ptr<CellVec3d> caliberBulidData;
};

NG_SMART_POINTER_TYPEDEF(ParamPack, NGParamPack);

// expand function
namespace PARAMPACK{
    void CompareSimilarity(const VectorVec5d&, const VectorVec5d&, size_t&);
    void ConcatTracedCurves(VectorVec5d&, VectorVec5d&, int);
}

#endif