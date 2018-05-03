/* SparseTracer
*	function: trace single neuron
*	Author : zhouhang, lishiwei
*	2015-10-28
*/
#ifndef LARGESPARSETRACEFILTER_H
#define LARGESPARSETRACEFILTER_H
#include <stack>
#include <utility>
#include <exception>
#ifdef _WIN32
#include <ctime>
#endif
#include "../../ineuronprocessobject.h"
#include "../../../ngtypes/basetypes.h"
#include "../../../ngtypes/volume.h"
#include "../../../ngtypes/ParamPack.h"
#include "../../../ngtypes/tree.h"
class LargeSparseTraceFilter;
class INeuronBigReader;
class BinaryFilter;
class SparseTraceFilter;
class SparseTraceAxonFilter;
NG_SMART_POINTER_TYPEDEF(INeuronBigReader, NGNeuronBigReader);
NG_SMART_POINTER_TYPEDEF(BinaryFilter, NGBinaryFilter);
NG_SMART_POINTER_TYPEDEF(SparseTraceFilter, NGSparseTraceFilter);
NG_SMART_POINTER_TYPEDEF(SparseTraceAxonFilter, NGSparseTraceAxonFilter);
NG_SMART_POINTER_TYPEDEF(LargeSparseTraceFilter, NGLargeSparseTraceFilter);
NG_SMART_POINTER_TYPEDEF(std::vector<VectorVec5d>, SMARTTREE);

class MyException :public std::exception
{
public:
    const char* what()const throw()  
    {
        return "ERROR! Please call zhouhang.\n";
    }
};


class LargeSparseTraceFilter : public INeuronProcessObject
{
public:
    static NGLargeSparseTraceFilter New(){return NGLargeSparseTraceFilter(new LargeSparseTraceFilter());}
    LargeSparseTraceFilter();
    ~LargeSparseTraceFilter();

    virtual ProcStatPointer Update();
    ConstIDataPointer GetOutput();//modification is allowed
    IDataPointer GetOutputTree();//modification is allowed


    //the prefix of axon means the parameter of axon trace
    void SetInputOldTree(IDataPointer);
    void SetInputMOSTD(NGNeuronBigReader);//mostd reader
    void SetSoma(IDataPointer);
    void SetParam(NGParamPack arg){if(arg == paramPack) return; paramPack = arg;}
    IDataPointer ReleaseCurrentImage();
    
    //user control
    ProcStatPointer Run() throw();
    void Step();
    void Stop();
    void Train();

protected:
    IDataPointer ReleaseData();
    void SetInput();
    void CalcPtCurveDirection(const VectorVec3d &ptCurve, Vec3d &fir);
    //bool CollectionBoundaryPt(std::vector < std::tuple<Vec3d, Vec3d, size_t> >& arg);
    bool IsInCurrentBlock(const Vec3d&);
    bool IsInCurrentBlock(const Vec5d&);
    void errorcheck();
    void CollectAllowableBoundCheckPt(const std::vector<VectorVec5d>&);
    void PopTraceInitInfoCurrentBlock();
    bool IsInOrNearCurrentBlock(const Vec3d& arg);
    bool Binary();
    bool CalcImageRange();
private:
    
    //
    bool isStartTrace_;
    enum TRACESTATUS{SOMA, AXONSTART, AXONCONTINUE};
    TRACESTATUS status_;
    int curStepNum_;
    //int xScale_, yScale_, zScale_;
    //init axon pos and direction
    Vec3d initAxonPosition_;
    Vec3d initAxonDirection_;
    Vec3i largeVolumeBoundary;
    Vec3d tmpBoundaryPtDirection;
    //data for function class
    NGParamPack paramPack;
    NGParamPack paramPackBakup_;
    NG_SMART_POINTER_DEFINE(Soma, tmpSoma);
    //function class
    NGSparseTraceFilter sparseFilter;
    NGSparseTraceAxonFilter sparseAxonFilter;
    NGBinaryFilter binaryFilter;
    //for most-d
    NGNeuronBigReader mostdReader;
    int sumIndex;
    //
    std::pair<Vec3d, Vec3d> nextInit_;
};

#endif // LARGESPARSETRACEFILTER_H
