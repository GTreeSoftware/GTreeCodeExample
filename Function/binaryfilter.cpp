#include "binaryfilter.h"
#include "../ngtypes/volume.h"
#ifdef _WIN32
#include <ctime>
#define GLOG_NO_ABBREVIATED_SEVERITIES
#include <glog/logging.h>
#else
#include <sys/time.h>
#endif

BinaryFilter::BinaryFilter()
{
    className_ = std::string("BinaryFilter");
    NG_SMART_POINTER_NEW(CVolume, m_Source, this);
    NG_SMART_POINTER_NEW(SVolume, m_Back, this);
    m_BinPtSet = VectorVec3i();//BinPtSetPointer(new VectorVec3i);
    binThreshold_ = 6.0;
}

BinaryFilter::~BinaryFilter()
{

}

ProcStatPointer BinaryFilter::Update()
{
    if(!m_Input){
        printf("error occurred in %s\n", className_.c_str());
        //LOG(ERROR) << "error occured in " << className_ << "when check m_Input.";
        MAKEPROCESSSTATUS(resSta, false, className_, "no input data.");
        return resSta;
    }
    if(m_Input->GetProcessObject()){
        ProcStatPointer res = m_Input->GetProcessObject()->Update();
        if(! res->success()){
            printf("error occurred in %s\n", className_.c_str());
            //LOG(ERROR) << "error occurred in "<<className_<<" when check m_Input.";
            return res;
        }
    }
    if(!Binary()) {
        printf("error occurred in %s\n", className_.c_str());
        //LOG(ERROR) << "error occured in " << className_ << "when check Binary.";
        MAKEPROCESSSTATUS(resSta, false, className_, "Binary operation failed.");
        return resSta;
    }
    MAKEPROCESSSTATUS(resSta, true, className_, "");
    return resSta;
}

ConstIDataPointer BinaryFilter::GetOutput()
{
    if(!m_Source) NG_SMART_POINTER_NEW(CVolume, m_Source, this);
    if(!m_Back) NG_SMART_POINTER_NEW(SVolume, m_Back, this);
    //if(!m_BinPtSet) m_BinPtSet = BinPtSetPointer(new VectorVec3i);
    return m_Source;
}

ConstIDataPointer BinaryFilter::GetBackNoiseImage()
{
    return m_Back;
}

VectorVec3i &BinaryFilter::GetBinPtSet()
{
    return m_BinPtSet;
}

INEURONPROCESSOBJECT_RELEASEDATA_IMPLE(BinaryFilter)

IDataPointer BinaryFilter::ReleaseBackNoiseImage()
{
    m_Back->ReleaseProcessObject();
    IDataPointer tData(m_Back);
    m_Back.reset();
    return tData;
}

//BinaryFilter::BinPtSetPointer BinaryFilter::ReleaseBinPtSet()
//{
//    return m_BinPtSet;
//}

//void BinaryFilter::SetThreadNum(int t)
//{
//    threadNum = t;
//}

void BinaryFilter::SetThreshold(double t)
{
    binThreshold_ = t;
}

bool BinaryFilter::Binary()
{
    if ( !m_Input || m_Input->GetIdentifyName() != std::string("Volume")) return false;
    NG_CREATE_DYNAMIC_CAST(const SVolume, tmpImg, m_Input);
    NG_CREATE_DYNAMIC_CAST(CVolume, tmpBin, m_Source);
    NG_CREATE_DYNAMIC_CAST(SVolume, tmpBack, m_Back);
    /*range*/
    int i, j, ij, num,k,l;
    int w = tmpImg->x();
    int h = tmpImg->y();
    int f = tmpImg->z();
    int wh =w * h;
    //int filterNum = 5;//2015-6-8
#ifdef _WIN32
    clock_t beg = clock();
#else
    timeval start1, end1;
    gettimeofday(&start1, 0);
#endif
    if (tmpBin->x() < 10  && tmpBack->x() < 10) {//if not computed
        /*initialization*/
        tmpBin->SetSize(w,h,f);
        tmpBin->SetResolution(tmpImg->XResolution(),tmpImg->YResolution(),tmpImg->ZResolution());
        tmpBack->SetSize(w, h, f);
        tmpBack->SetResolution(tmpImg->XResolution(),tmpImg->YResolution(),tmpImg->ZResolution());
        /*parameter*/
        int radius = 4;
        double tmp_value = 0.0f;
        double accu = double( (2*radius+1) * (2*radius+1) );/*model sum num*/
        double mean; CalcRandomMean(*tmpImg, mean);
        double minThrev = m_Input->DataType() == DATATYPE::IMAGE8 ? 400.0 : 600.0;

        ///---------------------------------20 filters--------------------------///
        omp_set_num_threads(paramPack->threadNum_);
#pragma omp parallel
#pragma omp for  private( i, j, num,k,l,tmp_value)
        for (ij = 0 ; ij < f ; ++ij){			//all frame
            double *m1 = new double[wh];
            double *sum = new double[ w - 2*radius ];
            memset(sum, 0, sizeof(double) * (w-2*radius));//model
            ///----------------------------------------------------
            for (j = 0; j < h; ++j)
                for (i = 0 ; i < w ; ++i)
                    m1[j * w + i] = std::min((double)tmpImg->GetPixel(i, j, ij), minThrev);//2016-4-14

            double *m2 = new double[wh];			//memset(m2, 0, sizeof(double) * H * W);

            for(num = 0 ; num < paramPack->filterNum_; ++num)
            {
                ///------------------first mean filter---------------------------------

                //inner part
                for(i = radius ; i < w - radius ; ++i){   //initialization array template
                    tmp_value = 0;
                    for ( k = i - radius ; k <= i + radius ; ++k){     // linear scan
                        for ( l = 0 ; l <= 2 * radius ; ++l)
                            tmp_value += m1[ l * w + k ];
                    }
                    m2[ radius * w + i ] = tmp_value / accu;
                    sum[ i - radius ] = tmp_value;							//finish linear template
                }

                for (j = radius+1 ; j < h - radius ; ++j){
                    for(i = radius; i < w - radius; ++i){
                        for( k = i - radius ; k <= i + radius ; ++k)
                            sum[i-radius] += m1[ k+ (j+radius)* w] - m1[k+(j-radius-1)*w ];
                        m2[j * w + i] = sum[ i-radius ] / accu;//get mean value
                    }
                }

                // boundray
                for (i = 0 ; i < w ; ++i){
                    for (l = 0 ; l < radius ;++ l){
                        //upper
                        m2[i+l*w] = (m1[i + l*w ] + m1[i  + (l+1)*w ]) / 2;
                        //lower
                        m2[i + (h - 1 -l ) * w]=(m1[i + (h - 1 -l) * w ]
                        + m1[i + (h - 2 -l) * w ]) / 2;
                    }
                }

                //left and right boundary
                for (j = 0 ; j < h  ; ++j){
                    for (l = 0 ; l < radius ;++ l ){
                        //left
                        m2[l + j * w] =(m1[l + j * w ]
                        + m1[l+1 + j * w ])/ 2 ;
                        //right
                        m2[w - 1-l + j * w] =(m1[w - 1-l + j * w ] +
                            m1[w - 2-l + j * w])/ 2 ;
                    }
                }

                memcpy(m1, m2, sizeof(double) * wh);		//there are some questions.
            }
            ///-------------------------------------finish filtering and Binary----------------------------
            // YY=MMt-YYs
            for ( j = 0; j < h; ++j ){//YY0=(YY>(1+threv*sqrt(YYs)));%(key points
                for (i = 0 ; i < w ; ++i){
                    if ( ( double(tmpImg->GetPixel(i,j,ij)) - m1[ j * w + i ]) > (1.0f + binThreshold_ * sqrt( m1[ j * w + i ] ) ) 
                        )//|| double(tmpImg->GetPixel(i,j,ij)) - m1[ j * w + i ] > threValue|| double(tmpImg->GetPixel(i,j,ij)) > ratio * m1[ j * w + i ]
                        tmpBin->GetPixel(i,j,ij) = 255;
                    //tmpBack->GetPixel(i,j,ij) = std::max(8,(int)m1[ j * w + i ]);//2015-6-8/*save background noise*/
                    tmpBack->GetPixel(i,j,ij) = std::max(1,(int)m1[ j * w + i ]);//2015-6-15 TDI079
                }
            }

            /*clear temp data*/
            delete[] m1;
            delete[] m2;
            delete[] sum;
        }//completed
    }else{
        //backImg has been computed
        omp_set_num_threads(paramPack->threadNum_);
#pragma omp parallel
#pragma omp for  private( i, j)
        for (ij = 0 ; ij < f ; ++ij){
            for ( j = 0; j < h; ++j ){//YY0=(YY>(1+threv*sqrt(YYs)));%(key points
                for (i = 0 ; i < w ; ++i){
                    if (  double(tmpImg->GetPixel(i,j,ij)) - tmpBack->GetPixel(i,j,ij) > (1.0f + binThreshold_ * sqrt( double(tmpBack->GetPixel(i,j,ij)) ) ) 
                        )//|| double(tmpImg->GetPixel(i,j,ij)) - double(tmpBack->GetPixel(i,j,ij)) > threValue|| double(tmpImg->GetPixel(i,j,ij)) > ratio * double(tmpBack->GetPixel(i,j,ij))
                        tmpBin->GetPixel(i,j,ij) = 255;
                }
            }
        }
    }

    ///----------------dont use multi-thread----------------///
    for(int ij = 0; ij < f; ++ij)
        for(int j = 0; j < h; ++j)
            for(int i = 0; i < w; ++i){
                if(0 < tmpBin->GetPixel(i,j,ij))
                    m_BinPtSet.push_back(Vec3i(i,j,ij));
            }

    ///--------------------get eclipsed time--------------------------
#ifdef _WIN32
    clock_t end = clock();
    printf("%d ms eclipsed in binarization. \n", int(end - beg));
    //LOG(INFO) << int(end - beg) <<"ms eclipsed in binarization.";
    
#else
    gettimeofday(&end1, 0);
    double timeuse=1000000*(end1.tv_sec-start1.tv_sec)+end1.tv_usec-start1.tv_usec;
    timeuse/=1000000;
    printf("%lf s eclipsed in binarization.\n", timeuse);
#endif
    /*update status*/
    printf("Binarying Finished!\n");
    LOG(INFO) <<"Binarying Finished!";
    printf("There are %d dots before Erosion.\n", (int)m_BinPtSet.size());
    //LOG(INFO) << "There are "<< (int)m_BinPtSet.size() << " dots before Erosion.";
    return true;
}

void BinaryFilter::CalcRandomMean(const SVolume& vol, double& mean)
{
    srand(0);
    int x, y, z;
    VecXd tmpData(1000);
    for (int i = 0; i < 1000; ++i) {
        x = int(double(rand()) / double(RAND_MAX) * double(vol.x()));
        y = int(double(rand()) / double(RAND_MAX) * double(vol.y()));
        z = int(double(rand()) / double(RAND_MAX) * double(vol.z()));
        x = std::max(std::min(x, vol.x() - 1), 0);
        y = std::max(std::min(y, vol.y() - 1), 0);
        z = std::max(std::min(z, vol.z() - 1), 0);
        tmpData(i) = double(vol(x, y, z));
    }
    mean = 1.5 * tmpData.mean();
    mean = mean < 600.0 ? 600.0 : 3096;
}


void BinaryFilter::CalcBackImage(const SVolume&src, SVolume&backimage, int filterNum){
	/*range*/
	int i, j, ij, num, k, l;
	int w = src.x();
	int h = src.y();
	int f = src.z();
	int wh = w * h;
	//int filterNum = 5;//2015-6-8
	/*initialization*/
	backimage.SetSize(w, h, f);
	backimage.SetResolution(src.XResolution(), src.YResolution(), src.ZResolution());
	backimage.SetZero();
	/*parameter*/
	int radius = 4;
	double tmp_value = 0.0f;
	double accu = double((2 * radius + 1) * (2 * radius + 1));/*model sum num*/
	//double mean; CalcRandomMean(src, mean);

	///---------------------------------20 filters--------------------------///
	omp_set_num_threads(omp_get_max_threads());
#pragma omp parallel
#pragma omp for  private( i, j, num,k,l,tmp_value)
	for (ij = 0; ij < f; ++ij){			//all frame
		double *m1 = new double[wh];
		double *sum = new double[w - 2 * radius];
		memset(sum, 0, sizeof(double) * (w - 2 * radius));//model
		///----------------------------------------------------
		for (j = 0; j < h; ++j)
			for (i = 0; i < w; ++i)
				m1[j * w + i] = std::min((double)src.GetPixel(i, j, ij), 1000.0);//2016-4-14

		double *m2 = new double[wh];			//memset(m2, 0, sizeof(double) * H * W);

		for (num = 0; num < filterNum; ++num)
		{
			///------------------first mean filter---------------------------------

			//inner part
			for (i = radius; i < w - radius; ++i){   //initialization array template
				tmp_value = 0;
				for (k = i - radius; k <= i + radius; ++k){     // linear scan
					for (l = 0; l <= 2 * radius; ++l)
						tmp_value += m1[l * w + k];
				}
				m2[radius * w + i] = tmp_value / accu;
				sum[i - radius] = tmp_value;							//finish linear template
			}

			for (j = radius + 1; j < h - radius; ++j){
				for (i = radius; i < w - radius; ++i){
					for (k = i - radius; k <= i + radius; ++k)
						sum[i - radius] += m1[k + (j + radius)* w] - m1[k + (j - radius - 1)*w];
					m2[j * w + i] = sum[i - radius] / accu;//get mean value
				}
			}

			// boundray
			for (i = 0; i < w; i++){
				for (l = 0; l < radius; l++){
					//upper
					m2[i + l*w] = m1[i + l*w];
					//lower
					m2[i + (h - 1 - l) * w] = m1[i + (h - 1 - l) * w];
				}
			}

			//left and right boundray
			for (j = 0; j < h; j++){
				for (l = 0; l < radius; l++){
					//left
					m2[l + j * w] = m1[l + j * w];
					//right
					m2[w - 1 - l + j * w] = m1[w - 1 - l + j * w];
				}
			}
			// boundray
			for (i = 0; i < w; i++){
				for (l = 0; l < radius; l++){
					//upper
					m2[i + l*w] = (m2[i + l*w] + m2[i + (l + 1)*w]) / 2;
					//lower
					m2[i + (h - 1 - l) * w] = (m2[i + (h - 1 - l) * w]
						+ m2[i + (h - 2 - l) * w]) / 2;
				}
			}

			//left and right boundray
			for (j = 0; j < h; j++){
				for (l = 0; l < radius; l++){
					//left
					m2[l + j * w] = (m2[l + j * w]
						+ m2[l + 1 + j * w]) / 2;
					//right
					m2[w - 1 - l + j * w] = (m2[w - 1 - l + j * w] +
						m2[w - 2 - l + j * w]) / 2;
				}
			}

			memcpy(m1, m2, sizeof(double) * wh);		//there are some questions.
		}
		///-------------------------------------finish filtering and Binary----------------------------
		// YY=MMt-YYs
		for (j = 0; j < h; ++j){//YY0=(YY>(1+threv*sqrt(YYs)));%(key points
			for (i = 0; i < w; ++i){
				backimage.GetPixel(i, j, ij) = std::round(std::max(0.0, m1[j * w + i]));//2015-6-15 TDI079
			}
		}

		/*clear temp data*/
		delete[] m1;
		delete[] m2;
		delete[] sum;
	}//completed	

}