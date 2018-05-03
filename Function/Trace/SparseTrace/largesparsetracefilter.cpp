/* SparseTracer
*	function: trace single neuron
*	Author : zhouhang, lishiwei
*	2015-10-28
*/
#include <QMessageBox>
#include <exception>
#include "largesparsetracefilter.h"
#include "Function/binaryfilter.h"
#include "sparsetracefilter.h"
#include "sparsetraceaxonfilter.h"
#include "../../../ngtypes/tree.h"
#include "../../../ngtypes/soma.h"
#include "Function/IO/imagereader.h"
#include "Function/IO/MOSTDReader.h"
#include "Function/volumealgo.h"
#include "Function/Trace/traceutil.h"
#include "Function/IO/imagewriter.h"
#include <omp.h>
#include "Function/ImageScaleFiter.h"
#define GLOG_NO_ABBREVIATED_SEVERITIES
#include <glog/logging.h>
#include <iostream>
using namespace google;

LargeSparseTraceFilter::LargeSparseTraceFilter()
{
    className_ = std::string("LargeSparseTraceFilter");
    NG_SMART_POINTER_NEW(NeuronPopulation, m_Source, this);
    status_ = AXONCONTINUE;
    isStartTrace_ = false;
    nextInit_.first.setZero();
    nextInit_.second.setZero();
}

LargeSparseTraceFilter::~LargeSparseTraceFilter()
{
}

ProcStatPointer LargeSparseTraceFilter::Update()
{
    return Run();
}


IDataPointer LargeSparseTraceFilter::GetOutputTree()
{
    if (!m_Source)
        NG_SMART_POINTER_DEFINE_NEW(NeuronPopulation, m_Source, this);
    return m_Source;
}


ConstIDataPointer LargeSparseTraceFilter::GetOutput()
{
    if(!m_Source)
        m_Source = std::shared_ptr<NeuronPopulation>(new NeuronPopulation(this));
    return m_Source;
}

IDataPointer LargeSparseTraceFilter::ReleaseData()
{
    m_Source->ReleaseProcessObject();
    IDataPointer tData(m_Source);
    m_Source.reset();
    return tData;
}

void LargeSparseTraceFilter::SetInputMOSTD( NGNeuronBigReader reader )
{
    mostdReader = reader;
}

ProcStatPointer LargeSparseTraceFilter::Run()
{
    if(!m_Source)  NG_SMART_POINTER_NEW(NeuronPopulation, m_Source, this);
    NG_SMART_POINTER_DEFINE(TreeCurve, tmpTree);
    NG_CREATE_DYNAMIC_CAST(NeuronPopulation, curTree, m_Source);
    isStartTrace_ = true;
    if (paramPack->dataMode_ == READMODE::MOSTD) {
        mostdReader->GetMaxRange(largeVolumeBoundary(0), largeVolumeBoundary(1), largeVolumeBoundary(2));
    }
    else{
        NG_CREATE_DYNAMIC_CAST(SVolume, m_OrigImg, paramPack->OrigImage);
        largeVolumeBoundary << m_OrigImg->x(), m_OrigImg->y(), m_OrigImg->z();
    }
    
    /*deep search*/
    IDataPointer tmpRelease;
    switch(status_){
        case AXONSTART:
        case AXONCONTINUE:
            {
                if(status_ != AXONSTART) {
                } else status_ = AXONCONTINUE;
                printf("%u remain boundary seed  \n", paramPack->activeTree->m_traceInitInfo.size());
                if (paramPack->activeTree->m_traceInitInfo.empty()) {
                    printf("no more boundary curve.\n");
                    MAKEPROCESSSTATUS(resSta, false, className_, "no more boundary curve.");
                    return resSta;
                }
                //pop trace init
                PopTraceInitInfoCurrentBlock();
                if (nextInit_.first.cwiseAbs().sum() < 1.0) return false;
                //prepare data*/
                auto &it = nextInit_;
                paramPack->activeTree->m_curInitPt =it.first;//update curInitPt
                Vec3d &tmpBoundaryPt = it.first;
                tmpBoundaryPtDirection = it.second;
                //stop when out of volume
                if (tmpBoundaryPt(0) < 0 || tmpBoundaryPt(0) >= largeVolumeBoundary(0)
                    || tmpBoundaryPt(1) < 0 || tmpBoundaryPt(1) >= largeVolumeBoundary(1)
                    || tmpBoundaryPt(2) < 0 || tmpBoundaryPt(2) >= largeVolumeBoundary(2)) {
                    printf("out of data boundary.\n");
                    MAKEPROCESSSTATUS(resSta, false, className_, "out of data boundary.");
                    return resSta;
                }
                //soma
                NG_SMART_POINTER_NEW_DEFAULT(Soma, tmpSoma);
                tmpSoma->push_back(Cell(0, tmpBoundaryPt(0) ,tmpBoundaryPt(1), tmpBoundaryPt(2), 1,1,1,1,1,1));
                printf("seed: %lf %lf %lf\n", tmpBoundaryPt(0) ,tmpBoundaryPt(1), tmpBoundaryPt(2));
                //save old paramPack. If failed, restore status
                //if (paramPack->isLocalTrace_) {
                //}
                //Calculate New image Range
                paramPackBakup_ = std::make_shared<ParamPack>();
                *paramPackBakup_ = *paramPack;
                if (!CalcImageRange()) {
                    //*paramPack = *paramPackBakup_;
                    MAKEPROCESSSTATUS(resSta, false, className_, "Image range error.");
                    return resSta;
                    //return false;
                }
                if (paramPack->dataMode_ == READMODE::MOSTD) {
                    /*read data*/
                    mostdReader->SetParam(paramPack);
                    //mostdReader->UpdateROI();
                    auto res = mostdReader->Update();
                    if (!res->success()){
                        *paramPack = *paramPackBakup_;
                        return res;
                    }
                    paramPack->OrigImage = mostdReader->ReleaseData();
                    paramPackBakup_->OrigImage = paramPack->OrigImage;
					paramPackBakup_->xMin_ = paramPack->xMin_;
					paramPackBakup_->xMax_ = paramPack->xMax_;
					paramPackBakup_->yMin_ = paramPack->yMin_;
					paramPackBakup_->yMax_ = paramPack->yMax_;
					paramPackBakup_->zMin_ = paramPack->zMin_;
					paramPackBakup_->zMax_ = paramPack->zMax_;
                }
                else{//crop image
                    NG_CREATE_DYNAMIC_CAST(SVolume, m_OrigImg, paramPack->OrigImage);
                    NG_SMART_POINTER_DEFINE_NEW(SVolume, localImage,0);
                    ExtractArea(*m_OrigImg, paramPack->xMin_, paramPack->xMax_, paramPack->yMin_, paramPack->yMax_, paramPack->zMin_, paramPack->zMax_, *localImage);
                    paramPack->OrigImage = localImage;
                }
                
                //map to local scaled range*/
                NG_CREATE_DYNAMIC_CAST(SVolume, m_OrigImg, paramPack->OrigImage);
                //map soma to local*/
                tmpSoma->GetCell(0).x -= paramPack->xMin_ ;
                tmpSoma->GetCell(0).y -= paramPack->yMin_ ;
                tmpSoma->GetCell(0).z -= paramPack->zMin_ ;
                //
                
                if (paramPack->xScale_ != 1 || paramPack->yScale_ != 1) {
                    IDataPointer scaleImage;
                    {
                        NGImageScaleFiter scaler = ImageScaleFiter::New();
                        scaler->SetInput(paramPack->OrigImage);
                        scaler->SetParam(paramPack);
                        auto scres = scaler->Update();
                        if (!scres->success()){
                            *paramPack = *paramPackBakup_;
                            return scres;
                        }
                        scaleImage = scaler->ReleaseData();
                        tmpSoma->GetCell(0).x /= paramPack->xScale_;
                        tmpSoma->GetCell(0).y /= paramPack->yScale_;
                        tmpSoma->GetCell(0).z /= paramPack->zScale_;
                    }
                    paramPack->OrigImage = scaleImage;
                }
                //else scaleImage = paramPack->OrigImage;
                //binary*/
                if (!Binary()){
                    *paramPack = *paramPackBakup_;
                    MAKEPROCESSSTATUS(resSta, false, className_, "Binary failed.");
                    return resSta;
                }
                //test
                //NGImageWriter writer = ImageWriter::New();
                //writer->SetOutputFileName("F:/qiao/nimei.tif");
                //writer->SetInput(paramPack->OrigImage);//OrigImage
                //writer->Update();
                //trace*/
                if(!sparseAxonFilter) sparseAxonFilter = SparseTraceAxonFilter::New();
                sparseAxonFilter->SetInput(paramPack->OrigImage);
                sparseAxonFilter->SetInputBack(paramPack->BackImage);
                sparseAxonFilter->SetInputBin(paramPack->BinImage);
                sparseAxonFilter->SetSoma(tmpSoma);
                sparseAxonFilter->SetParam(paramPack);
                sparseAxonFilter->SetInitialDirection(tmpBoundaryPtDirection);
                //clock_t debugbeg = clock();
                auto sres = sparseAxonFilter->Update();
                //clock_t debugend = clock();
                //printf("%d ms eclipsed in Total trace. \n", int(debugend - debugbeg));
                if (!sres->success()) {
                    *paramPack = *paramPackBakup_;
                    return sres;
                }
                tmpRelease = sparseAxonFilter->ReleaseData();
                NG_DYNAMIC_CAST(TreeCurve, tmpTree,  tmpRelease);
                Vec3d offSet(paramPack->xMin_, paramPack->yMin_, paramPack->zMin_);
                std::vector<VectorVec5d>& tmpCurve5d = tmpTree->GetCurve();
                if (tmpCurve5d.empty()){
                    *paramPack = *paramPackBakup_;
                    MAKEPROCESSSTATUS(resSta, false, className_, "trace line is empty.");
                    return resSta;
                }
                if (paramPack->xScale_ != 1 || paramPack->yScale_ != 1) {
                    for (auto &it : tmpCurve5d) {
                        for (auto &it2 : it) {
                            it2(0) *= paramPack->xScale_;
                            it2(1) *= paramPack->yScale_;
                        }
                    }
                }
                //map to global position
                size_t addID = paramPack->activeTree->m_curveList.size();
                tree<LineID> &topo = paramPack->activeTree->m_Topo;
                std::vector<VectorVec5d> &globalCurves = paramPack->activeTree->m_curveList;
                tree<LineID> traceTopo; traceTopo.set_head(LineID(std::numeric_limits<size_t>::max()));
                std::vector<VectorVec5d> traceLines;
                Vec3d root = tmpCurve5d.front().front().block(0, 0, 3, 1) + offSet;
                for(size_t k = 0; k < tmpCurve5d.size(); ++k){
                    if (tmpCurve5d[k].size() < 2) continue;
                    VectorVec5d tmpCurve;
                    Vec5d tmp;
                    for(size_t l = 0; l < tmpCurve5d[k].size(); ++l){
                        tmp = tmpCurve5d[k][l];
                        tmp.block(0,0,3,1) = tmp.block(0,0,3,1) + offSet;
                        tmpCurve.push_back(tmp);
                    }
                    traceLines.push_back(VectorVec5d());
                    traceLines.back().swap(tmpCurve);
                    traceTopo.append_child(traceTopo.begin(), LineID(k));
                }
                //sort tree
                if (traceLines.empty()) {
                    printf("There is no more signal to trace.\n");
                    *paramPack = *paramPackBakup_;
                    MAKEPROCESSSTATUS(resSta, true, className_, "There is no more signal to trace.");
                    return resSta;
                }
                KPTREE::SetTreeRoot(traceLines, traceTopo, root);
                for (auto treeIt = traceTopo.begin() + 1; treeIt != traceTopo.end(); ++treeIt) {
                    treeIt->id += addID;
                }
                //KPTREE::print_tree_bracketed(traceTopo);
                CollectAllowableBoundCheckPt(traceLines);
                //merge tree
                paramPack->curTraceTopoID_ = std::numeric_limits<size_t>::max();//parent id
                for (auto initIter = globalCurves.rbegin(); initIter != globalCurves.rend(); ++ initIter) {//find parent curve
                    if ( (initIter->back().block(0,0,3,1) - paramPack->activeTree->m_curInitPt).cwiseAbs().sum() < 1.0 ) {
                        paramPack->curTraceTopoID_ = std::distance(globalCurves.begin(), (++initIter).base());
                        break;
                    }
                }
                if (paramPack->curTraceTopoID_ == std::numeric_limits<size_t>::max()) {
                    NG_ERROR_MESSAGE("Why I cannot find parent curve?");
                    *paramPack = *paramPackBakup_;
                    MAKEPROCESSSTATUS(resSta, false, className_, "Why I cannot find parent curve?");
                    return resSta;
                    //throw MyException();
                }
                else{
                    auto topIter = std::find(topo.begin() + 1, topo.end(), LineID(paramPack->curTraceTopoID_));
                    auto neurontype = topIter->type;
                    for (auto iter = traceTopo.begin() + 1; iter != traceTopo.end(); ++iter)
                        iter->type = neurontype;
                }
                globalCurves.reserve(globalCurves.size() + traceLines.size());
                for (size_t k = 0; k < traceLines.size(); ++k) {
                    globalCurves.push_back(Line5d());
                    traceLines[k].swap(globalCurves.back());
                }
                //KPTREE::print_tree_bracketed(paramPack->activeTree->m_Topo);
                auto pIter = std::find(topo.begin() + 1, topo.end(), LineID(paramPack->curTraceTopoID_));//iter
                traceTopo.reparent(pIter, traceTopo.begin());
                //merge lines, new traced first line and trace initial lines.
                size_t traceFirstID = addID;//new traced line id
                if (!KPTREE::MergeTwoLine(*(paramPack->activeTree), traceFirstID, paramPack->curTraceTopoID_)){
                    NG_ERROR_MESSAGE("MergeTwoLine error.");
                    MAKEPROCESSSTATUS(resSta, false, className_, "MergeTwoLine error.");
                    system("pause");
                    return resSta;
                }
                //KPTREE::print_tree_bracketed(paramPack->activeTree->m_Topo);
                errorcheck();

                if (paramPack->OrigImage) paramPack->OrigImage.reset();//delete tmp image
                if (paramPack->BinImage) paramPack->BinImage.reset();
                if (paramPack->BackImage) paramPack->BackImage.reset();
                if (paramPack->isLocalTrace_) {
                    *paramPack = *paramPackBakup_;//restore completely
                    paramPackBakup_->allowableBoundCheckPtList_.swap(paramPack->allowableBoundCheckPtList_);
                }
				else{
					paramPack->OrigImage = paramPackBakup_->OrigImage;//just update
				}
                
            }//switch
    }
    printf("current block traced.\n");
    MAKEPROCESSSTATUS(resSta, true, className_, "");
    return resSta;
}

void LargeSparseTraceFilter::errorcheck()
{
    size_t maxId = 0;
    for (auto it = paramPack->activeTree->m_Topo.begin() + 1lu; it != paramPack->activeTree->m_Topo.end(); ++it) {
        auto curLineID = it->id;
        maxId = std::max(curLineID, maxId);
        const Line5d &line = paramPack->activeTree->m_curveList[curLineID];
        if (paramPack->activeTree->m_Topo.parent(it) == paramPack->activeTree->m_Topo.begin()){
        }
        else{
            size_t parentLineID = paramPack->activeTree->m_Topo.parent(it)->id;
            const Vec3d &headNode = line.front().block(0, 0, 3, 1);
            size_t nearestID = KPTREE::FindNearestID(headNode, paramPack->activeTree->m_curveList[parentLineID]);
            if (nearestID > 1000000){
                char line[256]; sprintf_s(line, "call zhouhang. parent id : %lu cur id : %lu\n", parentLineID, curLineID);
                NG_ERROR_MESSAGE(line);
                //printf("call zhouhang. parent id : %lu cur id : %lu\n", parentLineID, curLineID);
                system("pause");
            }
        }
    }
    if (maxId!= (paramPack->activeTree->m_curveList.size() - 1)) {
        printf("topo and curve are not right\n");
        system("pause");
    }
}

IDataPointer LargeSparseTraceFilter::ReleaseCurrentImage()
{
    return paramPack->OrigImage;
}

void LargeSparseTraceFilter::Stop()
{
    paramPack->activeTree->m_traceInitInfo.clear();
    m_Source.reset();
    paramPack->manualLabelCurve_.clear();
}

//
void LargeSparseTraceFilter::Step()
{
    Run();
}

void LargeSparseTraceFilter::CalcPtCurveDirection(const VectorVec3d &ptCurve, Vec3d &fir)
{
    int nx = int(ptCurve.size());
    Vec3d tmpFir;
    tmpFir.setZero();
    double tmp(0.0);
    Vec3d tmpVec;
    for (int i = 0; i < nx - 1; ++i){
        tmpVec = ptCurve[i + 1] - ptCurve[i];
        tmp = tmpVec.norm();
        tmpFir += tmp * tmpVec;
    }
    fir = tmpFir.normalized();
}

void LargeSparseTraceFilter::SetInputOldTree(IDataPointer op)
{
    m_Source = op;
}

bool LargeSparseTraceFilter::IsInCurrentBlock(const Vec3d& arg)
{
    return arg(0) >= double(paramPack->xMin_) && arg(0) < double(paramPack->xMax_) && arg(1) >= double(paramPack->yMin_) && arg(1) < double(paramPack->yMax_) && arg(2) >= double(paramPack->zMin_) && arg(2) < double(paramPack->zMax_);
}

bool LargeSparseTraceFilter::IsInCurrentBlock(const Vec5d& arg)
{
    return arg(0) >= double(paramPack->xMin_) && arg(0) < double(paramPack->xMax_) && arg(1) >= double(paramPack->yMin_) && arg(1) < double(paramPack->yMax_) && arg(2) >= double(paramPack->zMin_) && arg(2) < double(paramPack->zMax_);
}

bool LargeSparseTraceFilter::IsInOrNearCurrentBlock(const Vec3d& arg)
{
    return arg(0) >= double(paramPack->xMin_ - paramPack->boundaryDistanceThreshold_)
        && arg(0) < double(paramPack->xMax_ + paramPack->boundaryDistanceThreshold_)
        && arg(1) >= double(paramPack->yMin_ - paramPack->boundaryDistanceThreshold_)
        && arg(1) < double(paramPack->yMax_ + paramPack->boundaryDistanceThreshold_)
        && arg(2) >= double(paramPack->zMin_ - paramPack->boundaryDistanceThreshold_)
        && arg(2) < double(paramPack->zMax_ + paramPack->boundaryDistanceThreshold_);
}

void LargeSparseTraceFilter::CollectAllowableBoundCheckPt(const std::vector<VectorVec5d>&arg)
{
    auto &list = paramPack->allowableBoundCheckPtList_;
    list.clear();
    for (auto &curve : arg) {
        list.push_back(curve.back().block(0, 0, 3, 1));
    }
}

void LargeSparseTraceFilter::Train()
{
    if (sparseAxonFilter && paramPack->activeTree) {
        const auto &allLines = paramPack->activeTree->m_curveList;
        std::vector<VectorVec5d> tmp;
        Vec5d tmp5d;
        //search node in current image area
        bool flag = false;
        for (size_t k = 0; k < allLines.size(); ++k) {
            flag = false;
            for (size_t kk = 0; kk < allLines[k].size(); ++kk) {
                if (IsInCurrentBlock(allLines[k][kk])) {
                    if (!flag) {
                        tmp.push_back(VectorVec5d());
                        flag = true;
                    }
                    tmp5d = allLines[k][kk];
                    tmp5d(0) -= paramPack->xMin_;
                    tmp5d(1) -= paramPack->yMin_;
                    tmp5d(2) -= paramPack->zMin_;
                    for (size_t kkk = 0; kkk < 2llu; ++kkk) tmp.back().push_back(tmp5d);//to get more train samples
                }
                else{
                    if (flag) flag = false;
                }
            }
        }
        sparseAxonFilter->Train(tmp);
    }
}

void  LargeSparseTraceFilter::PopTraceInitInfoCurrentBlock()
{
    auto &traceInfo = paramPack->activeTree->m_traceInitInfo;
    auto revIt = traceInfo.rbegin();
    for (; revIt != traceInfo.rend(); ++revIt){
        if (IsInOrNearCurrentBlock((*revIt).first)) break;
    }
    if (revIt != traceInfo.rend()) {
        auto iter = (++revIt).base();
        std::swap(nextInit_, *iter);
        traceInfo.erase(iter);//
    }
    else{
        if (!traceInfo.empty()) {
            std::swap(nextInit_, traceInfo.back());
            traceInfo.pop_back();
        }
        else {
            nextInit_.first.setZero();
            nextInit_.second.setZero();
        }

    }
}

bool LargeSparseTraceFilter::Binary()
{
    NGBinaryFilter filter = BinaryFilter::New();
    filter->SetInput(paramPack->OrigImage);
    filter->SetParam(paramPack);
    filter->SetThreshold(paramPack->axonBinaryThreshold_);//it is independent. Maybe I should set as flag.
    if (!filter->Update()->success())
        return false;
    double iter = paramPack->axonBinaryThreshold_;
    while (filter->GetBinPtSet().size() > 100000){//2016-6-7
        iter += 1.0;
        filter = BinaryFilter::New();
        filter->SetInput(paramPack->OrigImage);
        filter->SetParam(paramPack);
        filter->SetThreshold(iter);
        if(!filter->Update()->success()) return false;
    }
    paramPack->BinImage = filter->ReleaseData();
    paramPack->BackImage = filter->ReleaseBackNoiseImage();
    if (!paramPack->BinImage || !paramPack->BackImage){
        printf("cannot binary \n");
        return false;
    }
    return true;
}

bool LargeSparseTraceFilter::CalcImageRange()
{
    int xOffset = paramPack->xExtract_;
    int yOffset = paramPack->yExtract_;
    int zOffset = paramPack->zExtract_;
    paramPack->xMin_ = tmpSoma->GetCell(0).x - xOffset + int(tmpBoundaryPtDirection(0) * double(paramPack->xExtract_) * 0.95);
    paramPack->xMax_ = tmpSoma->GetCell(0).x + xOffset + int(tmpBoundaryPtDirection(0) * double(paramPack->xExtract_) * 0.95);
    paramPack->yMin_ = tmpSoma->GetCell(0).y - yOffset + int(tmpBoundaryPtDirection(1) * double(paramPack->yExtract_) * 0.95);
    paramPack->yMax_ = tmpSoma->GetCell(0).y + yOffset + int(tmpBoundaryPtDirection(1) * double(paramPack->yExtract_) * 0.95);
    paramPack->zMin_ = tmpSoma->GetCell(0).z - zOffset + int(tmpBoundaryPtDirection(2) * double(paramPack->zExtract_) * 0.95);
    paramPack->zMax_ = tmpSoma->GetCell(0).z + zOffset + int(tmpBoundaryPtDirection(2) * double(paramPack->zExtract_) * 0.95);
    if (!paramPack->isLocalTrace_) {
        LimitRange(0, largeVolumeBoundary(0) - 1, 0, largeVolumeBoundary(1) - 1, 0, largeVolumeBoundary(2) - 1,
            paramPack->xMin_, paramPack->xMax_, paramPack->yMin_, paramPack->yMax_, paramPack->zMin_, paramPack->zMax_);
    }
    else{//limit in old image range
        LimitRange(paramPackBakup_->xMin_, paramPackBakup_->xMax_-1, paramPackBakup_->yMin_,
            paramPackBakup_->yMax_-1, paramPackBakup_->zMin_, paramPackBakup_->zMax_-1,
            paramPack->xMin_, paramPack->xMax_, paramPack->yMin_, paramPack->yMax_, paramPack->zMin_, paramPack->zMax_);
        if (paramPack->xMax_ - paramPack->xMin_ < 10 || paramPack->yMax_ - paramPack->yMin_ < 10
            || paramPack->zMax_ - paramPack->zMin_ < 10) {
            printf("Local Run point near boundary\n");
            return false;
        }
    }
    
    Vec6i tmpVec6i;
    tmpVec6i << paramPack->xMin_, paramPack->xMax_, paramPack->yMin_, paramPack->yMax_, paramPack->zMin_, paramPack->zMax_;
    printf(" range: %d %d %d %d %d %d \n",
        paramPack->xMin_, paramPack->xMax_, paramPack->yMin_, paramPack->yMax_, paramPack->zMin_, paramPack->zMax_);
    return true;
}


