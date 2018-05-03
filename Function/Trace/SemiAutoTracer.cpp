#include "SemiAutoTracer.h"
#include "Function/Trace/traceutil.h"
#include "Function/IO/imagewriter.h"
#include "Function/NGUtility.h"
#include "Function\Trace\SemiAutoConnect.h"

#include <iostream>
#include <algorithm>


SemiAutoTracer::SemiAutoTracer()
{
    className_ = std::string("SemiAutoTracer");
    m_Source = std::shared_ptr<TreeCurve>(new TreeCurve(this));//no use
}


SemiAutoTracer::~SemiAutoTracer()
{
}

INEURONPROCESSOBJECT_RELEASEDATA_IMPLE(SemiAutoTracer)
//IDataPointer SemiAutoTracer::ReleaseData()
//{
//    m_Source->ReleaseProcessObject();
//    IDataPointer tData(m_Source);
//    m_Source.reset();
//    return tData;
//}
INEURONPROCESSOBJECT_GETOUTPUT_IMPLE(SemiAutoTracer, NeuronPopulation)
//ConstIDataPointer SemiAutoTracer::GetOutput()
//{
//    if (!m_Source)
//        m_Source = std::shared_ptr<NeuronPopulation>(new NeuronPopulation(this));
//    return m_Source;
//}

ProcStatPointer SemiAutoTracer::Update()
{
    if (!m_Input){
        printf("error occurred in %s\n", className_.c_str());
        MAKEPROCESSSTATUS(resSta, false, className_, "no input data.");
        return resSta;
    }
    if (m_Input->GetProcessObject()){//|| !m_Soma
        ProcStatPointer res = m_Input->GetProcessObject()->Update();
        if (!res->success()){
            printf("error occurred in %s\n", className_.c_str());
            //LOG(ERROR) << "error occurred in "<<className_<<" when check m_Input.";
            return res;
        }
    }
    origImgPointer =  std::dynamic_pointer_cast<const Volume<unsigned short>>(m_Input);//2015-6-8
    if (!origImgPointer ) {
        printf("error occurred in %s\n", className_.c_str());
        MAKEPROCESSSTATUS(resSta, false, className_, "Input image is wrong.");
        return resSta;
    }
	
    beginPt_ -= Vec3d(xOffSet_, yOffSet_, zOffSet_);
    endPt_ -= Vec3d(xOffSet_, yOffSet_, zOffSet_);
    Vec3d iBegPt; iBegPt << NGUtility::Round(beginPt_(0)), NGUtility::Round(beginPt_(1)), NGUtility::Round(beginPt_(2));
    Vec3d iEndPt; iEndPt << NGUtility::Round(endPt_(0)), NGUtility::Round(endPt_(1)), NGUtility::Round(endPt_(2));

	SemiAutoConnect connectTrace(iBegPt, iEndPt, origImgPointer,paramPack->axonSemiautoConnectResampleThreshold);
	connectTrace.normThreshold = this->paramPack->axonSemiautoTraceThreshold / 100;
	connectTrace.search();//2018-1-18

	//MyHistogram _his(*origImgPointer);
	//_his.exec();//for test 2018-2-8
	
	if (!connectTrace.errorFlag)
	{
		traceLine_.clear();
		for (int i = 0; i <= connectTrace.resampleLine.size() - 1; ++i)
		{
			Vec5d tem = NGUtility::MakeVec5d(connectTrace.resampleLine[i](0), connectTrace.resampleLine[i](1), connectTrace.resampleLine[i](2));
			traceLine_.push_back(tem);
		}

		Vec5d offSet = NGUtility::MakeVec5d(xOffSet_, yOffSet_, zOffSet_);
		for (auto & it : traceLine_) {
			it += offSet;
		}
		MAKEPROCESSSTATUS(resSta, true, className_, "");
		return resSta;
	}
	else
	{
		traceLine_.clear();
		MAKEPROCESSSTATUS(resSta, false, className_, "Threshold value error! Please reset!");
		return resSta;
	}
}

