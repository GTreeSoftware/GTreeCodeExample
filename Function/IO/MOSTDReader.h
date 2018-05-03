/* SparseTracer
*	function: trace single neuron
*	Author : zhouhang, liyuxin
*	2015-10-28
*/
#ifndef MOSTDREADER_H
#define MOSTDREADER_H
#include <memory>
#include <string>
#include "../ineuronio.h"
#include "ngtypes/basetypes.h"
#include "ngtypes/ParamPack.h"

template<typename T> class Volume;
typedef Volume<unsigned char> CVolume;
typedef Volume<unsigned short> SVolume;
class INeuronDateObject;
class MostImage;
class MostImage16;
class MOSTDReader;
typedef std::shared_ptr<MOSTDReader> NGMOSTDReader;

//the output is unsigned char volume
class MOSTDReader : public INeuronBigReader
{
public:
#ifdef _WIN32
	typedef std::tr1::shared_ptr<SVolume> VolumePointer;
#else
	typedef std::shared_ptr<SVolume> VolumePointer;//2015-6-8
#endif
	static NGMOSTDReader New(){return NGMOSTDReader(new MOSTDReader());}
	MOSTDReader();
	~MOSTDReader();

	virtual bool SetInputFileName(const std::string&);
	virtual ProcStatPointer Update();
    virtual IDataPointer ReleaseData();
    virtual ConstIDataPointer GetOutput();
	void SetParam(NGParamPack arg){paramPack = arg;}
	
	virtual void GetMaxRange( int& x, int& y, int& z );/*the size of mostd dataq*/
	//void SetDataType(IMAGETYPE arg){dataType_ = arg;}/*1 is read 8bit image, 2 is read 16bit image*/
	//DATATYPE GetDataType()const {return paramPack->dataType_;}
	void SetMaxValue(int arg){paramPack->maxValue_ = arg;}/*the max value of mostd data. maybe 4096*/
	//int GetLevel()const{return paramPack->mostdLevel_;}
    void ClearCache();
    virtual DATATYPE GetImageType()const{ return dataType; }

protected:
	//--- forbidden use--//
	void SetInput(ConstIDataPointer&){}
	virtual bool SetOutputFileName(const std::string&){ return false; }
	void SetResolution(double, double, double);/**/
    bool UpdateROI();
	bool TransferMOSTD2Volume();/*transfer 1d data into volume structure*/
	//void SetRealResolution( double res, int origSz, int& scale, int& sz );/*set the resolution of result image */
	//int SetScale(double res);
    
private:
    DATATYPE dataType = DATATYPE::IMAGE16;
	MostImage* reader;/*mostd reader written by liyuxin*/
	MostImage16* reader16;/*16bit mostd reader written by zhouhang*/
	NGParamPack paramPack;
};
#endif
