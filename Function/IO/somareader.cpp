/* SparseTracer
*	function: trace single neuron
*	Author : zhouhang
*	2015-10-28
*/
#include "somareader.h"
#include <cstdlib>
#include <cmath>
#include "../../ngtypes/soma.h"
#ifdef _WIN32
#pragma warning(disable : 4996)
#endif

SomaReader::SomaReader()
{
    className_ = std::string("SomaReader");
    m_Source = std::tr1::shared_ptr<Soma>(new Soma(this));
    isFileValid = false;
}

SomaReader::~SomaReader()
{

}

bool SomaReader::SetInputFileName(const std::string &path)
{
    if(!path.empty()){
        filename = path;
        isFileValid = true;
        return true;
    }
    return false;
//    FILE* fp = fopen(path.c_str(), "r");
//    if(NULL == fp)
//        printf("can not open file : %s.\n error from %s", path.c_str(), identifyName.c_str());
//    else {
//        filename = path;
//        isFileValid = true;
//    }
//    fclose(fp);
}

ProcStatPointer SomaReader::Update()
{
    Cell tempCell;
    char line[256];
	SomaPointer tempSoma = std::dynamic_pointer_cast<Soma >( m_Source );
    
    if(isFileValid){
        FILE* fp = fopen(filename.c_str(), "r");
        if(!fp){
            printf("error occurred in %s\n", className_.c_str());
            MAKEPROCESSSTATUS(resSta, false, className_, "Cannot create file.");
            return resSta;
        }

		int except;
        while (!feof(fp)) {
            if (fgets(line, 256, fp) == 0) continue;
            if (!isdigit( line[0] )) continue;
            sscanf(line, "%d %d %lf %lf %lf %lf -1", &tempCell.ID, &except, &tempCell.x, &tempCell.y, &tempCell.z, &tempCell.r);
            //2015-2-5
            //--tempCell.x;
            //--tempCell.y; 
            //--tempCell.z;
			tempCell.x /= paramPack->xRes_;
			tempCell.y /= paramPack->yRes_;
			tempCell.z /= paramPack->zRes_;
            tempSoma->push_back(tempCell);
        }
		fclose(fp);
    }
    else {
        printf("error occurred in %s\n", className_.c_str());
        MAKEPROCESSSTATUS(resSta, false, className_, "no input file path.");
        return resSta;
    }
    MAKEPROCESSSTATUS(resSta, true, className_, "");
    return resSta;
}

IDataPointer SomaReader::ReleaseData()
{
    m_Source->ReleaseProcessObject();
    IDataPointer tData = m_Source;
    m_Source.reset();
    return tData;
}

ConstIDataPointer SomaReader::GetOutput()
{
    if(!m_Source)
#ifdef _WIN32
        m_Source = std::tr1::shared_ptr<Soma>(new Soma(this));
#else
        m_Source = std::shared_ptr<Soma>(new Soma(this));
#endif
    return m_Source;
}

//void SomaReader::SetImageResolution( double x, double y, double z )
//{
//    xres_ = x;
//    yres_ = y;
//    zres_ = z;
//}

