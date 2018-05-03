#ifndef NEUROGPSTREEFILTER_H
#define NEUROGPSTREEFILTER_H
#include "../ineuronprocessobject.h"
#include "../ngtypes/ineurondataobject.h"
#include "../../ngtypes/basetypes.h"
#include "../../ngtypes/volume.h"
#include "../../ngtypes/ParamPack.h"
#include "../../ngtypes/tree.h"
class NeuroGPSTreeFilter;
class BinaryFilter;
class TraceFilter;
class BridgeBreaker;
class NeuroTreeCluster;
NG_SMART_POINTER_TYPEDEF(NeuroGPSTreeFilter, NGNeuroGPSTreeFilter);
NG_SMART_POINTER_TYPEDEF(BinaryFilter, NGBinaryFilter);
NG_SMART_POINTER_TYPEDEF(TraceFilter, NGTraceFilter);
NG_SMART_POINTER_TYPEDEF(BridgeBreaker, NGBridgeBreaker);
NG_SMART_POINTER_TYPEDEF(NeuroTreeCluster, NGNeuronTreeCluster);

class NeuroGPSTreeFilter: public INeuronProcessObject
{
public:
    static NGNeuroGPSTreeFilter New(){ return NGNeuroGPSTreeFilter(new NeuroGPSTreeFilter()); }
    NeuroGPSTreeFilter();
    ~NeuroGPSTreeFilter();
    INEURONPROCESSOBJECT_DEFINE
    /*virtual bool Update();
    ConstIDataPointer GetOutput();
    IDataPointer ReleaseData();*/
    void SetParam(NGParamPack arg){ param_ = arg; }
    void SetSoma(ConstIDataPointer arg);
   
protected:
    int SetScale(double res);
    bool GenerateDendritePopulation(std::vector<NeuronTreePointer>& neuronTree, std::vector<size_t> &treeIDList, VectorVec3d& somaList);
private:
    NGParamPack param_;//orig, bin, back
    //Data
    IDataPointer scaleImg_;
    NG_SMART_POINTER_DEFINE(VectorVec3i, m_binPtSet);
    ConstIDataPointer soma_;
    IDataPointer tree_;
    IDataPointer treeConInfo_;
    //filter
    NGBinaryFilter binaryFilter;
    NGTraceFilter traceFilter;
    NGBridgeBreaker breaker;
    NGNeuronTreeCluster treeCluster;

};
#endif
