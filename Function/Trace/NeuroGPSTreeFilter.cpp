#include "NeuroGPSTreeFilter.h"
#include "../binaryfilter.h"
#include "tracefilter.h"
#include "bridgebreaker.h"
#include "neurotreecluster.h"
#include "../../ngtypes/soma.h"
#include "../ImageScaleFiter.h"


NeuroGPSTreeFilter::NeuroGPSTreeFilter()
{
}


NeuroGPSTreeFilter::~NeuroGPSTreeFilter()
{
}

ProcStatPointer NeuroGPSTreeFilter::Update()
{
    if (!param_->OrigImage || !soma_) {
        printf("trace data not prepared.\n");
        MAKEPROCESSSTATUS(resSta, false, className_, "trace data not prepared.");
        return resSta;
    }
    NG_CREATE_DYNAMIC_CAST(const Soma, tmpSoma, soma_);
    if (tmpSoma->IsEmpty()) {
        printf("trace data: soma not prepared\n");
        //NG_ERROR_MESSAGE("trace data: soma not prepared.");
        //return false;
    }
    NG_CREATE_DYNAMIC_CAST(SVolume, tmpOrigImg, param_->OrigImage);
    //scale image
    NGImageScaleFiter scaler = ImageScaleFiter::New();
    scaler->SetInput(param_->OrigImage);
    scaler->SetParam(param_);
    scaler->Update();
    scaleImg_ = scaler->ReleaseData();
    //int xScale = SetScale(param_->xRes_);
    //int yScale = SetScale(param_->yRes_);
    //int zScale = SetScale(param->zRes_);
    //if (xScale * yScale ==1) {
    //    scaleImg_ = param_->OrigImage;
    //} else{
    //    scaleImg_ = std::make_shared<SVolume>();
    //    NG_CREATE_DYNAMIC_CAST(SVolume, tmpScaleImg, scaleImg_);
    //    tmpScaleImg->SetSize(tmpOrigImg->x() / xScale, tmpOrigImg->y() / yScale, tmpOrigImg->z() );
    //    tmpScaleImg->SetResolution(tmpOrigImg->XResolution() * xScale, tmpOrigImg->YResolution() * yScale, tmpOrigImg->ZResolution());
    //    unsigned short temp;
    //    double ratio = 1.0 / double(xScale * yScale);
    //    for (int z = 0; z < tmpScaleImg->z(); ++z) {
    //        for (int i = 0; i < tmpScaleImg->x(); ++i)//
    //            for (int j = 0; j < tmpScaleImg->y(); ++j){
    //                temp = 0;
    //                for (int ii = 0; ii < xScale; ++ii){
    //                    for (int jj = 0; jj < yScale; ++jj){
    //                        temp += tmpOrigImg->GetPixel(i*xScale + ii, j * yScale + jj, z);
    //                    }
    //                }
    //                tmpScaleImg->operator ()(i, j, z) = (unsigned short)(double(temp) * ratio );//2015-6-8//int(double(temp) * ratio);
    //                if (param_->dataType_ == IMAGE8) {
    //                    tmpScaleImg->operator ()(i, j, z) *= 4;
    //                }
    //            }
    //    }
    //}

    //scale soma
    IDataPointer scaleSoma;
    if (param_->xScale_ == 1 && param_->yScale_ == 1) {
        scaleSoma = soma_;
    }
    else{
        scaleSoma = std::make_shared<Soma>();
        NG_CREATE_DYNAMIC_CAST(Soma, tmpScaleSoma, scaleSoma);
        for (size_t i = 0; i < tmpSoma->size(); ++i) {
            tmpScaleSoma->push_back(Soma::MakeCell(tmpSoma->GetCell(i).x / double(param_->xScale_), 
                tmpSoma->GetCell(i).y / double(param_->yScale_), tmpSoma->GetCell(i).z));
        }
    }
    
    //binary
    binaryFilter = BinaryFilter::New();
    binaryFilter->SetInput(scaleImg_);//param->OrigImage
    binaryFilter->SetParam(param_);
    binaryFilter->SetThreshold(param_->binThreshold_);
    //binaryFilter->SetThreValue(255 * 4);//2015-8-13
    auto bres = binaryFilter->Update();
    if (!bres->success())
        return bres;
    param_->BinImage = binaryFilter->ReleaseData();
    param_->BackImage = binaryFilter->ReleaseBackNoiseImage();
    printf("Image Binarization complete.\n");
    //
    traceFilter = TraceFilter::New();
    traceFilter->SetInput(scaleImg_);//param->OrigImage
    traceFilter->SetInputBack(param_->BackImage);
    traceFilter->SetInputBin(param_->BinImage);
    traceFilter->SetSoma(scaleSoma);//soma
    traceFilter->SetThreValue(param_->traceValue_);
    traceFilter->SetParam(param_);
    auto tres = traceFilter->Update();
    if (!tres->success())
        return tres;
    printf("Trace First Step complete.\n");
    //
    tree_ = traceFilter->ReleaseData();
    treeConInfo_ = traceFilter->GetConnect();
    {
       /* NG_CREATE_DYNAMIC_CAST(TreeCurve, nima, tree_);
        FILE * fp = fopen("F:/video/hehe.swc", "w");
        int id = 0;

        for (auto &curve : nima->GetCurve()) {
            fprintf(fp, "%d 1 %lf %lf %lf 1 -1\n", ++id, curve[0](0), curve[0](1), curve[0](2));
            for (size_t k = 0; k < curve.size(); ++k) {
                fprintf(fp, "%d 1 %lf %lf %lf 1 %d\n", ++id, curve[k](0), curve[k](1), curve[k](2), id - 1);
            }
        }
        fclose(fp);
        system("pause");*/
    }

    breaker = BridgeBreaker::New();
    breaker->SetInput(traceFilter->GetConnect());
    breaker->SetInputSoma(scaleSoma);//soma
    breaker->SetInputOrigImage(scaleImg_);//param->OrigImage
    breaker->SetInputTree(tree_);
    auto brres = breaker->Update();
    if (!brres->success())
        return brres;
    printf("Trace Second Step complete.\n");
    //
    NG_SMART_POINTER_DEFINE(INeuronDataObject, breakCon);
    breakCon = breaker->ReleaseData();
    treeCluster = NeuroTreeCluster::New();
    treeCluster->SetInput(breakCon);
    treeCluster->SetInputDeleteID(breaker->getResultIndex());
    treeCluster->SetInputCurve(tree_);
    treeCluster->SetInputSoma(scaleSoma);//soma
    treeCluster->SetInputOrigImage(scaleImg_);//param->OrigImage
    auto trres = treeCluster->Update();
    if (!trres->success())
        return trres;
    printf("Trace complete.\n");
    m_Source = treeCluster->ReleaseData();
    //add offset
    NG_CREATE_DYNAMIC_CAST(NeuronPopulation, tmpTreeList, m_Source);
    std::vector<NeuronTreePointer> &treeList = tmpTreeList->m_pop;
    //make soma
    NG_CREATE_DYNAMIC_CAST(Soma, tmpScaleSoma, scaleSoma);
    VectorVec3d copySoma;
    if (false) {//tmpScaleSoma->size() > 0
        for (size_t i = 0; i < tmpScaleSoma->m_Source.size(); ++i) {
            copySoma.push_back(Vec3d(tmpScaleSoma->GetCell(i).x, 
                tmpScaleSoma->GetCell(i).y, tmpScaleSoma->GetCell(i).z));
        }
    }
    else{
        for (size_t i = 0; i < treeList.size(); ++i) {
            copySoma.push_back(Vec3d(treeList[i]->m_curveList[0][0](0),
                treeList[i]->m_curveList[0][0](1), treeList[i]->m_curveList[0][0](2) ));
        }
    }
    if (tmpTreeList->m_pop.size() != copySoma.size()) {
        printf("whether you are debugging or the soma tree is wrong ? \n");
        MAKEPROCESSSTATUS(resSta, false, className_, "whether you are debugging or the soma tree is wrong ?");
        return resSta;
    }
    {
        /*char line[256]; int fileId = 0;
        for (auto &it : tmpTreeList->m_pop) {
        sprintf_s(line, "F:/nima_%d.swc", fileId++);
        auto &allLine = it->m_curveList;
        int xuhao = 0;
        FILE* fp = fopen(line, "w");
        for (auto &line : allLine) {
        for (size_t k = 0; k < line.size(); ++k) {
        if (k == 0) fprintf(fp, "%d 1 %lf %lf %lf 1 -1\n", ++xuhao, line[k](0)*param_->xRes_, line[k](1)*param_->yRes_, line[k](2));
        else fprintf(fp, "%d 1 %lf %lf %lf 1 %d\n", ++xuhao, line[k](0)*param_->xRes_, line[k](1)*param_->yRes_, line[k](2), xuhao - 1);
        }
        }
        fclose(fp);
        }*/
        /*std::cout << "tmpScaleSoma" << std::endl;
        for (auto i = 0; i < tmpScaleSoma->size(); ++i) {
        std::cout << tmpScaleSoma->GetCell(i).x << " " << tmpScaleSoma->GetCell(i).y << " " << tmpScaleSoma->GetCell(i).z << " " << std::endl;
        }

        std::cout << "copySoma£º" << std::endl;
        for (auto &it : copySoma) {
        std::cout << it(0) << " " << it(1) << " " << it(2) << std::endl;
        }
        system("pause");*/
        /*std::cout << "tree size" << std::endl;
        for (size_t i = 0; i < treeList.size(); ++i) {
        std::cout << treeList[i]->m_curveList.size() << std::endl;
        }*/
    }

    //make tree
    std::vector<size_t > treeIDList(tmpTreeList->m_pop.size());
    {
        size_t i = 0;
        std::for_each(treeIDList.begin(), treeIDList.end(), [&](size_t& arg){arg = i; ++i; });
    }
    while (copySoma.size() >0lu && treeIDList.size() >0lu) {
        if (copySoma.size() != treeIDList.size()) {
            printf("error in NeuroGPSTreeFilter Make dendrites.\n");
            return false;
        }
        if (!GenerateDendritePopulation(tmpTreeList->m_pop, treeIDList, copySoma))
        {
            NG_ERROR_MESSAGE("error in GenerateDendritePopulation");
            MAKEPROCESSSTATUS(resSta, false, className_, "error in GenerateDendritePopulation.");
            return resSta;
        }
    }
    //map  to global
    for (size_t i = 0; i < treeList.size(); ++i) {
        for (size_t j = 0; j < (*treeList[i]).m_curveList.size(); ++j) {
            for (size_t ij = 0; ij < (*treeList[i]).m_curveList[j].size(); ++ij) {
                (*treeList[i]).m_curveList[j][ij](0) *= param_->xScale_;
                (*treeList[i]).m_curveList[j][ij](1) *= param_->yScale_;
                (*treeList[i]).m_curveList[j][ij](0) += param_->xMin_;
                (*treeList[i]).m_curveList[j][ij](1) += param_->yMin_;
                (*treeList[i]).m_curveList[j][ij](2) += param_->zMin_;
            }

        }
    }
    MAKEPROCESSSTATUS(resSta, true, className_, "");
    return resSta;
}

INEURONPROCESSOBJECT_GETOUTPUT_IMPLE(NeuroGPSTreeFilter, NeuronPopulation)
//ConstIDataPointer NeuroGPSTreeFilter::GetOutput()
//{
//    if (!m_Source)
//        NG_SMART_POINTER_NEW(NeuronPopulation, m_Source, this);
//    return m_Source;
//}
INEURONPROCESSOBJECT_RELEASEDATA_IMPLE(NeuroGPSTreeFilter)
//IDataPointer NeuroGPSTreeFilter::ReleaseData()
//{
//    m_Source->ReleaseProcessObject();
//    IDataPointer tData(m_Source);
//    m_Source.reset();
//    return tData;
//}

void NeuroGPSTreeFilter::SetSoma(ConstIDataPointer arg)
{
    if (!param_ || !arg) {
        printf("cannont input Soma because param is emtpy.\n");
        return;
    }
    //set offset
    NG_SMART_POINTER_NEW_DEFAULT(Soma, soma_);
    NG_CREATE_DYNAMIC_CAST(Soma, tmpArg, arg);
    NG_CREATE_DYNAMIC_CAST(Soma, tmpSoma, soma_);
    tmpSoma->GetAllCell() = tmpArg->GetAllCell();
    for (size_t i = 0; i < tmpSoma->size(); ++i) {
        tmpSoma->GetCell(i).x -= param_->xMin_;
        tmpSoma->GetCell(i).y -= param_->yMin_;
        tmpSoma->GetCell(i).z -= param_->zMin_;
    }
}

int NeuroGPSTreeFilter::SetScale(double res)
{
    //2015-8-13
    int scale;
    if (res > 0.25 && res <= 0.5) {
        scale = 2;
    }
    else if (res <= 0.25) {
        scale = 3;
    }
    else{
        scale = 1;
    }
    return scale;
}

bool NeuroGPSTreeFilter::GenerateDendritePopulation(std::vector<NeuronTreePointer> &neuronTree, std::vector<size_t> &treeIDList, VectorVec3d& somaList)
{
    bool flag = false;
    for (size_t i = 0; i < somaList.size(); ++i) {
        //printf("%d\n", i);
        for (size_t j = 0; j < treeIDList.size(); ++j) {
            if (KPTREE::SetTreeRoot(*neuronTree[treeIDList[j]], somaList[i])) {
                neuronTree[treeIDList[j]]->m_curInitPt =somaList[i];
                flag = true;
                std::cout << "G complete " << somaList[i](0) << " " << somaList[i](1) << " " << somaList[i](2) << std::endl;
                somaList.erase(somaList.begin() + i);
                treeIDList.erase(treeIDList.begin() + j);
                break;
            }
        }
        if (flag) break;
        else{
            bool allempty = true;
            for (size_t j = 0; j < treeIDList.size(); ++j) {
                if (!neuronTree[treeIDList[j]]->m_curveList.empty()) {
                    allempty = false;
                    break;
                }
            }
            if (!allempty) {
                std::cout << "nima size:" << neuronTree[treeIDList[0]]->size() << std::endl;
                std::cout << "nima :" << somaList[i](0) << " " << somaList[i](1) << " " << somaList[i](2) << std::endl;
                Vec3d root(somaList[i]);
                for (auto &it : neuronTree[treeIDList[0]]->m_curveList) {
                    for (int k = 0; k < it.size(); ++k){
                        if ((it[k].block(0, 0, 3, 1) - root).norm() < 3.0 ) {
                            std::cout << "wocao " << k << " :" << it[k](0) << " " << it[k](1) << " " << it[k](2) << std::endl;
                            break;
                        }
                    }
                }
            }
            else{
                std::cout << "rest trees are empty. " << std::endl;
            }
        }
    }
    return flag;
}

