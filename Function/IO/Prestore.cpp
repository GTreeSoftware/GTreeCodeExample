#include <QMessageBox>
#include <iostream>
#include "../../ngtypes/volume.h"
#include "Prestore.h"
#include "../../Thread/StoreThread.h"
#include "Function/IO/MOSTDReader.h"

using namespace std;
Prestore::Prestore()
{
    className_ = std::string("Prestore");
    paramPackCp = std::make_shared<ParamPack>();
    worker = Q_NULLPTR;
}


Prestore::~Prestore()
{
    if (worker) {
        if (worker->isRunning()){
            worker->quit();
            worker->wait();
        }
        delete worker;
    }
}

bool Prestore::IsImageDequeValid()
{
    return !imageDeque_.empty();
}


void Prestore::StartThread()
{
    if (worker!= Q_NULLPTR && !worker->isFinished()){
        return;
    }
    FowardTrace();
    paramPackCp->xMin_ = FrameCoo[0] - paramPack->xExtract_;
    paramPackCp->yMin_ = FrameCoo[2] - paramPack->yExtract_;
    paramPackCp->zMin_ = FrameCoo[4] - paramPack->zExtract_;
    paramPackCp->xMax_ = FrameCoo[1] + paramPack->xExtract_;
    paramPackCp->yMax_ = FrameCoo[3] + paramPack->yExtract_;
    paramPackCp->zMax_ = FrameCoo[5] + paramPack->zExtract_;

    printf("INPre %d %d %d %d %d %d\n", paramPackCp->xMin_, paramPackCp->xMax_, paramPackCp->yMin_, paramPackCp->yMax_, paramPackCp->zMin_, paramPackCp->zMax_);
    mostdReader->SetParam(paramPackCp);
    if (worker == Q_NULLPTR) {
        worker = new StoreThread(this);
        connect(worker, SIGNAL(finished()), this, SLOT(EndThread_Slot()));
        worker->setParam(paramPackCp, mostdReader);
    }

    worker->start();
}

void Prestore::EndThread_Slot()
{
    if (!worker->RetVal()) {
        printf("please restart \n");
        return;
    }
    if (forceThreadStop_) {
        mostdReader->ReleaseData();
        forceThreadStop_ = false;
        return;
    }
    ImageData = mostdReader->ReleaseData();
    Vec6i tmp; tmp << paramPackCp->xMin_, paramPackCp->yMin_, paramPackCp->zMin_, paramPackCp->xMax_, paramPackCp->yMax_, paramPackCp->zMax_;
    imageDeque_.push_front(std::tie(dynamic_pointer_cast<SVolume>(ImageData), tmp, Vec3i(traverseIter->curLineID, traverseIter->curBeg, traverseIter->curEnd)));
    //++readingBox;
    cout << imageDeque_.size() << endl;
    if (imageDeque_.size() > storageImageMaxNum_)
    {
        //imageDeque_.pop_back();
        printf("hello\n");
        return;
    }
    if (imageDeque_.size() < storageImageMaxNum_){
        cout << imageDeque_.size() << "nima" << endl;
        StartThread();
    }
    emit CacheComplete_Signal();
}

void  Prestore::FowardTrace()
{
    if (directionFlag_) {
        traverseIter->Next(numPoints_);
    }
    else{
        traverseIter->Back(numPoints_);
    }
    auto &beg = traverseIter->curBeg;
    auto &fin = traverseIter->curEnd;
    auto &curLineID = traverseIter->curLineID;
    auto &allLines = paramPack->activeTree->m_curveList;
    //auto &flagList = paramPack->activeTree->traverseFlag_;
    auto &curLine = allLines[curLineID];
    
    FrameCoo << 10000000, 0, 10000000, 0, 10000000, 0;
    for (int j = beg; j < fin; ++j){
        if (curLine[j][0] < FrameCoo[0]){
            FrameCoo[0] = curLine[j][0];
        }
        if (curLine[j][0] > FrameCoo[1]){
            FrameCoo[1] = curLine[j][0];
        }
        if (curLine[j][1] < FrameCoo[2]){
            FrameCoo[2] = curLine[j][1];
        }
        if (curLine[j][1] > FrameCoo[3]){
            FrameCoo[3] = curLine[j][1];
        }
        if (curLine[j][2] < FrameCoo[4]){
            FrameCoo[4] = curLine[j][2];
        }
        if (curLine[j][2] > FrameCoo[5]){
            FrameCoo[5] = curLine[j][2];
        }
    }
    //AddCurLineID();
    printf("curId:%d  curPoint%d\n", curLineID, beg);
}


bool Prestore::GetNextImage(IDataPointer &ptr)//, size_t arg
{
    if (directionFlag_ == false) {
        if (worker->isRunning()){
            printf("Please waiting for back end thread finished.\n");
            forceThreadStop_ = true;
            return false;
        }
        imageDeque_.clear();
        traverseIter->curLineID = traverseIterCp->curLineID;
        traverseIter->curBeg = traverseIterCp->curBeg;
        traverseIter->curEnd = traverseIterCp->curEnd;
    }
    directionFlag_ = true;
    if (IsImageDequeValid()==false){
        StartThread();
        return false;
    }
    else
    {
        auto &curImageBlock = imageDeque_.back();
        ptr = std::get<0>(curImageBlock);
        paramPack->xMin_ = std::get<1>(curImageBlock)(0);
        paramPack->xMax_ = std::get<1>(curImageBlock)(3);
        paramPack->yMin_ = std::get<1>(curImageBlock)(1);
        paramPack->yMax_ = std::get<1>(curImageBlock)(4);
        paramPack->zMin_ = std::get<1>(curImageBlock)(2);
        paramPack->zMax_ = std::get<1>(curImageBlock)(5);
        traverseIterCp->curLineID = std::get<2>(curImageBlock)(0);
        traverseIterCp->curBeg = std::get<2>(curImageBlock)(1);
        traverseIterCp->curEnd = std::get<2>(curImageBlock)(2);
        //stamp
        auto &curLine = paramPack->activeTree->m_curveList[traverseIterCp->curLineID];
        for (int k = traverseIterCp->curBeg; k <= traverseIterCp->curEnd; ++k)
            curLine[k](4) = 1.0;
        //auto &flagList = paramPack->activeTree->traverseFlag_;
        //flagList[traverseIterCp->curLineID].first = std::min(traverseIterCp->curBeg, flagList[traverseIterCp->curLineID].first);
        //flagList[traverseIterCp->curLineID].second = std::max(traverseIterCp->curEnd, flagList[traverseIterCp->curLineID].second);
        imageDeque_.pop_back();

        StartThread();

        return true;
    }
}


void Prestore::SetParam(NGNeuronBigReader arg2, NGParamPack arg3)
{
    mostdReader = arg2;
    paramPack = arg3;
    *paramPackCp = *paramPack;
}

bool Prestore::GetPreviousImage(IDataPointer &ptr)
{
    if (directionFlag_ == true) {
        if (worker->isRunning()){
            printf("Please waiting for back end thread finished.\n");
            forceThreadStop_ = true;
            return false;
        }
        imageDeque_.clear();//the sort is error
        traverseIter->curLineID = traverseIterCp->curLineID;
        traverseIter->curBeg = traverseIterCp->curBeg;
        traverseIter->curEnd = traverseIterCp->curEnd;
    }
    directionFlag_ = false;
    if (IsImageDequeValid() == false){
        StartThread();
        return false;
    }
    else
    {
        auto &curImageBlock = imageDeque_.back();
        ptr = std::get<0>(curImageBlock);
        paramPack->xMin_ = std::get<1>(curImageBlock)(0);
        paramPack->xMax_ = std::get<1>(curImageBlock)(3);
        paramPack->yMin_ = std::get<1>(curImageBlock)(1);
        paramPack->yMax_ = std::get<1>(curImageBlock)(4);
        paramPack->zMin_ = std::get<1>(curImageBlock)(2);
        paramPack->zMax_ = std::get<1>(curImageBlock)(5);
        traverseIterCp->curLineID = std::get<2>(curImageBlock)(0);
        traverseIterCp->curBeg = std::get<2>(curImageBlock)(1);
        traverseIterCp->curEnd = std::get<2>(curImageBlock)(2);
        //stamp
        auto &curLine = paramPack->activeTree->m_curveList[traverseIterCp->curLineID];
        for (int k = traverseIterCp->curBeg; k <= traverseIterCp->curEnd; ++k)
            curLine[k](4) = 1.0;
        //auto &flagList = paramPack->activeTree->traverseFlag_;
        //flagList[traverseIterCp->curLineID].first = std::min(traverseIterCp->curBeg, flagList[traverseIterCp->curLineID].first);
        //flagList[traverseIterCp->curLineID].second = std::max(traverseIterCp->curEnd, flagList[traverseIterCp->curLineID].second);
        imageDeque_.pop_back();
        StartThread();

        return true;
    }
}

bool Prestore::Initial()
{
    if (!paramPack->activeTree) return false;
    auto &at = paramPack->activeTree;
    auto &mc = at->m_curveList;
    traverseIter = new RangeIterator(&(at->m_Topo), &mc);
    traverseIterCp = new RangeIterator(&(at->m_Topo), &mc);
    traverseIter->curBeg = preBeg;
    traverseIter->curEnd = preEnd;
    traverseIter->curLineID = preLineID;
    if (traverseIter->curBeg != std::numeric_limits<int>::max()) {
        auto &topo = at->m_Topo;
        auto curIter = std::find(topo.begin() + 1, topo.end(), traverseIter->curLineID);
        traverseIter->treeIter = curIter;
    }
    //auto &flag_ = paramPack->activeTree->traverseFlag_;
    //if (flag_.size() != mc.size()) {
    //    flag_.resize(mc.size());
    //    std::for_each(flag_.begin(), flag_.end(), [](std::pair<int, int> &arg){arg.first = 100000; arg.second = 0; });
    //}
    //else{//check error and update. If tree is modified, then traverse flag should be resize
    //    for (size_t k = 0; k < mc.size(); ++k) {
    //        if (flag_[k].first >= int(mc[k].size())) {
    //            flag_[k].first = 100000;
    //            flag_[k].second = 0;
    //        }
    //        if (flag_[k].second >= int(mc[k].size()) && flag_[k].first < int(mc[k].size())) {
    //            flag_[k].second = int(mc[k].size()) - 1;
    //        }
    //    }
    //}
    forceThreadStop_ = false;
    return true;
}

bool Prestore::Finish()
{
    forceThreadStop_ = true;
    if (worker && worker->isRunning()) {
        worker->wait();
    }
    if (!paramPack->activeTree){
        NG_ERROR_MESSAGE("nima");
        return false;
    }
    preBeg = traverseIterCp->curBeg;
    preEnd = traverseIterCp->curEnd;
    preLineID = traverseIterCp->curLineID;
    delete traverseIter;
    traverseIter = NULL;
    delete traverseIterCp;
    traverseIterCp = NULL;
    //mostdReader->ClearCache();
    if (mostdReader && paramPack->dataMode_ == READMODE::MOSTD) {
        NG_CREATE_DYNAMIC_CAST(MOSTDReader, tmpMost, mostdReader);
        if (tmpMost) {
            tmpMost->ClearCache();
        }
    }
    //paramPack->activeTree->traverseFlag_.clear();
    imageDeque_.clear();
    return true;
}

void Prestore::SetTraversePosition(int lineID, int beg)
{
    if (traverseIter) {
        if (worker && worker->isRunning()){
            printf("Please waiting for back end thread finished.\n");
            forceThreadStop_ = true;
            return;
        }
        directionFlag_ = true;
        imageDeque_.clear();
        traverseIter->curLineID = lineID;
        auto &topo = paramPack->activeTree->m_Topo;
        auto curIter = std::find(topo.begin() + 1, topo.end(), traverseIter->curLineID);
        if (curIter == topo.end() || curIter.node == NULL) {
            traverseIter->curLineID = 0;
            curIter = std::find(topo.begin() + 1, topo.end(), traverseIter->curLineID);
        }
        traverseIter->treeIter = curIter;
        int sz = int(paramPack->activeTree->m_curveList[lineID].size());
        traverseIter->curBeg = std::max(0, std::min(sz - 41, beg));
        traverseIter->curEnd = std::min(sz-1, beg + numPoints_);
        auto &beg = traverseIter->curBeg;
        auto &fin = traverseIter->curEnd;
        auto &curLineID = traverseIter->curLineID;
        //auto &topo = paramPack->activeTree->m_Topo;
        auto &allLines = paramPack->activeTree->m_curveList;
        //auto &flagList = paramPack->activeTree->traverseFlag_;
        auto &curLine = allLines[curLineID];
        //stamp
        for (int k = beg; k <= fin; ++k)
            curLine[k](4) = 1.0;
        //flagList[curLineID].first = std::min(beg, flagList[curLineID].first);
        //flagList[curLineID].second = std::min(fin, flagList[curLineID].second);
        FrameCoo << 10000000, 0, 10000000, 0, 10000000, 0;
        for (int j = beg; j < fin; ++j){
            if (curLine[j][0] < FrameCoo[0]){
                FrameCoo[0] = curLine[j][0];
            }
            if (curLine[j][0] > FrameCoo[1]){
                FrameCoo[1] = curLine[j][0];
            }
            if (curLine[j][1] < FrameCoo[2]){
                FrameCoo[2] = curLine[j][1];
            }
            if (curLine[j][1] > FrameCoo[3]){
                FrameCoo[3] = curLine[j][1];
            }
            if (curLine[j][2] < FrameCoo[4]){
                FrameCoo[4] = curLine[j][2];
            }
            if (curLine[j][2] > FrameCoo[5]){
                FrameCoo[5] = curLine[j][2];
            }
        }
        //AddCurLineID();
        printf("curId:%d  curPoint%d\n", curLineID, beg);
        paramPackCp->xMin_ = FrameCoo[0] - 20;
        paramPackCp->yMin_ = FrameCoo[2] - 20;
        paramPackCp->zMin_ = FrameCoo[4] - 10;
        paramPackCp->xMax_ = FrameCoo[1] + 20;
        paramPackCp->yMax_ = FrameCoo[3] + 20;
        paramPackCp->zMax_ = FrameCoo[5] + 10;

        printf("INPre %d %d %d %d %d %d\n", paramPackCp->xMin_, paramPackCp->xMax_, paramPackCp->yMin_, paramPackCp->yMax_, paramPackCp->zMin_, paramPackCp->zMax_);
        mostdReader->SetParam(paramPackCp);
        if (worker == Q_NULLPTR) {
            worker = new StoreThread(this);
            connect(worker, SIGNAL(finished()), this, SLOT(EndThread_Slot()));
            worker->setParam(paramPackCp, mostdReader);
        }

        worker->start();
    }
}

void Prestore::Reset()
{
    preBeg = preEnd = preLineID = 0;
    if (paramPack && paramPack->activeTree) {
        auto &curveList = paramPack->activeTree->m_curveList;
        for (auto &it1 : curveList) {
            for (auto &it2 : it1) {
                it2(4) = 0.0;
            }
        }
    }
}


void Prestore::RangeIterator::Next(int arg)
{
    if (treePointer == NULL || curvePointer == NULL) return;
    if (treeIter == treePointer->end()) return;
    if (std::numeric_limits<int>::max() == curBeg) {
        treeIter = treePointer->begin() + 1;
        curLineID = int(treeIter->id);
        int sz = int((*curvePointer)[curLineID].size()) - 1;
        curBeg = 0;
        curEnd = std::min(curBeg + arg, sz);
    }
    else{
        int sz = int((*curvePointer)[curLineID].size()) - 1;
        curBeg += arg;
        if (curBeg > sz) {
            ++treeIter;
            if (treeIter == treePointer->end()) {//loop
                treeIter = treePointer->begin() + 1;
            }
            curLineID = int(treeIter->id);
            sz = int((*curvePointer)[curLineID].size()) - 1;
            curBeg = 0;
            curEnd = std::min(curBeg + arg, sz);
        }
        else{
            curEnd = std::min(curBeg + arg, sz);
        }
    }
}

void Prestore::RangeIterator::Back(int arg)
{
    if (treePointer == NULL || curvePointer == NULL) return;
    if (std::numeric_limits<int>::max() == curBeg) {
        return;
    }
    else{
        //int sz = int(curvePointer[curLineID].size()) - 1;
        curBeg -= arg;
        curEnd -= arg;
        if (curBeg <0) {
            --treeIter;
            if (treeIter == treePointer->begin()) {
                treeIter = --treePointer->end();
            }
            curLineID = int(treeIter->id);
            int sz = int((*curvePointer)[curLineID].size()) - 1;
            curEnd = sz;
            curBeg = std::max(0, curEnd - arg);            
        }
    }
}
