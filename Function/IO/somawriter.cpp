/* SparseTracer
*	function: trace single neuron
*	Author : zhouhang
*	2015-10-28
*/
#include "somawriter.h"
#include <cstdlib>
#include "../../ngtypes/soma.h"
#pragma warning(disable:4996)

SomaWriter::SomaWriter()
{
    className_ = std::string("SomaWriter");
}

SomaWriter::~SomaWriter()
{

}

bool SomaWriter::SetOutputFileName(const std::string &path)
{
    fileName = path;
    return true;
}

ProcStatPointer SomaWriter::Update()
{
    if(!m_Input || m_Input->GetIdentifyName() != std::string("Soma")){
        printf("error occurred in %s\n", className_.c_str());
        MAKEPROCESSSTATUS(resSta, false, className_, "no input data.");
        return resSta;
    }
    if(m_Input->GetProcessObject()){
        ProcStatPointer res = m_Input->GetProcessObject()->Update();
        if (!res->success()){
            printf("error occurred in %s\n", className_.c_str());
            //LOG(ERROR) << "error occurred in "<<className_<<" when check m_Input.";
            return res;
        }
    }
    FILE* fp = fopen(fileName.c_str(), "w");
    if(!fp) {
        printf("error occurred in %s\n", className_.c_str());
        return false;
    }
    std::tr1::shared_ptr<const Soma> tmpsoma = std::dynamic_pointer_cast<const Soma>(m_Input);

    size_t num = tmpsoma->size();
    int index(1);
    const Cell& tmpcell = tmpsoma->GetCell(0);
    fprintf(fp,"%d 1 %lf %lf %lf %lf -1", index++, tmpcell.x, tmpcell.y, tmpcell.z, tmpcell.r);
    for (size_t i = 1; i < num; ++i){
        fprintf(fp, "\n");
        const Cell& tmpCell = tmpsoma->GetCell(i);
        fprintf(fp,"%d %lf %lf %lf %lf 1 -1", index++, tmpCell.x, tmpCell.y, tmpCell.z, tmpCell.r);
    }
    fclose(fp);
    MAKEPROCESSSTATUS(resSta, true, className_, "");
    return resSta;
}

