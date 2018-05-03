/*
 * Copyright (c)2013-2015  Zhou Hang, Shaoqun Zeng, Tingwei Quan
 * Britton Chance Center for Biomedical Photonics, Huazhong University of Science and Technology
 * All rights reserved.
 */
#ifndef TREEREADER_H
#define TREEREADER_H
#include <memory>
#include <string>
#include <vector>
#include "../ineuronio.h"
#include "ngtypes/ParamPack.h"
#pragma warning(disable:4996)
class INeuronDataObject;
class NeuronPopulation;
class TreeReader;
typedef std::shared_ptr<TreeReader> NGTreeReader;
class TreeReader :
    public INeuronIO
{
public:
    typedef std::shared_ptr<NeuronPopulation> NeuronPopulationPointer;
    static NGTreeReader New(){return NGTreeReader(new TreeReader());}
    TreeReader(void);
    ~TreeReader(void);
    void SetInputFileNameList(const std::vector<std::string>&);
    void SetTraverseFile(const std::vector<std::string>&);
    INEURONPROCESSOBJECT_DEFINE
        /*bool Update();
        IDataPointer ReleaseData();
        ConstIDataPointer GetOutput();*/

    
    bool ReadTraverseFile(const std::string &file, std::vector<Line5d>&);
    void SetParam(NGParamPack arg){ param = arg; }
private:
    void SetInput(ConstIDataPointer){}
    bool SetInputFileName(const std::string &){ return false; }
    bool SetOutputFileName(const std::string&){ return false; }
    bool isFileValid;
    std::vector<std::string> fileNameList;
    std::vector<std::string> traverseNameList;
    NGParamPack param;
    //double xRes_, yRes_, zRes_;
    //2015-8-13
    //double xScale_, yScale_, zScale_;
    //NGParamPack paramPack;//TODO
};
#endif
