#include <omp.h>
#include "most_data16.h"
#include <vector>
#include "Function/IO/imagewriter.h"
#include "ngtypes/volume.h"


bool MostImage16::select_IOR_file( int xbegin,int xend,int ybegin,int yend,int zbegin,int zend,int level )
{

    if (xend-xbegin+1>this->sz0||yend-ybegin+1>this->sz1||zend-zbegin+1>sz2)
    {
        //qDebug("out of size,error");
        std::cout << "out of size,error" << std::endl;
        return false;
    }
    if (xend<xbegin||yend<ybegin||zend<zbegin)
    {
        std::cout << "out of size,error" << std::endl;
        return false;
    }

    
    //const char *path;
    //QByteArray text;
    std::string text;
    //QStringList list_of_index;
    std::vector<std::string> list_of_index;
    int x_index_b=xbegin;
    int x_index_e=xend;
    int y_index_b=ybegin;
    int y_index_e=yend;
    int z_index_b=zbegin;
    int z_index_e=zend;
    /*int x_temp_b;
    int x_temp_e;
    int y_temp_b;
    int y_temp_e;
    int z_temp_b;
    int z_temp_e;*/
    int z_shift;
    int y_shift;
    int x_shift;

    for (int i=1;i<level;i++)
    {
        x_index_b /=2;
        x_index_e /=2;
        y_index_b /=2;
        y_index_e /=2;
        z_index_b /=2;
        z_index_e /=2;
    }

    int width_IOR = x_index_e - x_index_b+1;
    int height_IOR = y_index_e - y_index_b +1;
    int depth_IOR = z_index_e - z_index_b +1;
    this->IOR_sz0=width_IOR;
    this->IOR_sz1=height_IOR;
    this->IOR_sz2=depth_IOR;
    this->origin_IOR_x = xbegin;
    this->origin_IOR_y = ybegin;
    this->origin_IOR_z = zbegin;
    this->rez_IOR_x = std::pow(2.0,level-1)*rez_x;
    this->rez_IOR_y = std::pow(2.0, level-1)*rez_y;
    this->rez_IOR_z = std::pow(2.0, level-1)*rez_z;
    
    uint16 *buffer = new uint16[width_IOR*height_IOR*depth_IOR];  //申请总的内存开销
    //uint32 *slice = new uint32[int(BLOCK_SIZE*BLOCK_SIZE)];   //申请一张的内存开销，即每个slice的内存
    //uint16 *data = new uint16[int(BLOCK_SIZE*BLOCK_SIZE)];
    //uint16 *block_temp = new uint16[int(BLOCK_SIZE*BLOCK_SIZE*DEPTH_SIZE)];
    int compress(0);

    int x_block_b = x_index_b/int(BLOCK_SIZE);
    int x_block_e = x_index_e/int(BLOCK_SIZE);
    int y_block_b = y_index_b/int(BLOCK_SIZE);
    int y_block_e = y_index_e/int(BLOCK_SIZE);
    int z_block_b = z_index_b/int(DEPTH_SIZE);
    int z_block_e = z_index_e/int(DEPTH_SIZE);


    x_shift = x_index_b%int(BLOCK_SIZE);
    y_shift = y_index_b%int(BLOCK_SIZE);
    z_shift = z_index_b%int(DEPTH_SIZE);

    //std::cout<<x_block_b<<" "<<x_block_e<<" "<<y_block_b<<" "<<y_block_e<<" "<<z_block_b<<" "<<z_block_e<<std::endl;

    //int index_of_xy;
    int i, j, k;
//#pragma omp parallel for num_threads(2)  private( i, j, k)
    for (k = x_block_b; k <= x_block_e; k++)
    {
        uint16 *data = new uint16[int(BLOCK_SIZE*BLOCK_SIZE)];
        std::string str;
        int x_temp_b;
        int x_temp_e;
        int y_temp_b;
        int y_temp_e;
        int z_temp_b;
        int z_temp_e;
        for ( j = y_block_b; j <=y_block_e; j++)
        {
            for (i = z_block_b; i <= z_block_e; i++)
            {
                str=this->file_path;
                str.append("/level");
                std::stringstream a8;
                a8<<int(level);
                str.append(a8.str().data());
                str.append("_data/z");
                std::stringstream a9;
                a9<<int(i);
                str.append(a9.str().data());
                str.append("/y");
                std::stringstream a11;
                a11<<int(j);
                str.append(a11.str().data());
                str.append("/");
                std::stringstream a12;
                a12<<int(k);
                str.append(a12.str().data());
                str.append("_");
                std::stringstream a19;
                a19<<int(j);
                str.append(a19.str().data());
                str.append("_");
                std::stringstream a20;
                a20<<int(i);
                str.append(a20.str().data());
                str.append(".tif");
#pragma omp critical
                {
                //list_of_index.push_back(str);   //list
                }
                if (k==x_block_b&&k!=x_block_e)
                {
                    x_temp_b=x_index_b%int(BLOCK_SIZE);
                    x_temp_e=int(BLOCK_SIZE)-1;
                }
                if (k!=x_block_b&&k!=x_block_e)
                {
                    x_temp_b=0;
                    x_temp_e=int(BLOCK_SIZE)-1;
                }
                if (k!=x_block_b&&k==x_block_e)
                {
                    x_temp_b=0;
                    x_temp_e=x_index_e%int(BLOCK_SIZE);
                }
                if (k==x_block_b&&k==x_block_e)
                {
                    x_temp_b=x_index_b%int(BLOCK_SIZE);
                    x_temp_e=x_index_e%int(BLOCK_SIZE);
                }
                if (j==y_block_b&&j!=y_block_e)
                {
                    y_temp_b=y_index_b%int(BLOCK_SIZE);
                    y_temp_e=int(BLOCK_SIZE)-1;
                }
                if (j!=y_block_b&&j!=y_block_e)
                {
                    y_temp_b=0;
                    y_temp_e=int(BLOCK_SIZE)-1;
                }
                if (j!=y_block_b&&j==y_block_e)
                {
                    y_temp_b=0;
                    y_temp_e=y_index_e%int(BLOCK_SIZE);
                }
                if (j==y_block_b&&j==y_block_e)
                {
                    y_temp_b=y_index_b%int(BLOCK_SIZE);
                    y_temp_e=y_index_e%int(BLOCK_SIZE);
                }
                if (i==z_block_b&&i!=z_block_e)
                {
                    z_temp_b=z_index_b%int(DEPTH_SIZE);
                    z_temp_e=int(DEPTH_SIZE)-1;
                }
                if (i!=z_block_b&&i!=z_block_e)
                {
                    z_temp_b=0;
                    z_temp_e=int(DEPTH_SIZE)-1;
                }
                if (i!=z_block_b&&i==z_block_e)
                {
                    z_temp_b=0;
                    z_temp_e=z_index_e%int(DEPTH_SIZE);
                }
                if (i==z_block_b&&i==z_block_e)
                {
                    z_temp_b=z_index_b%int(DEPTH_SIZE);
                    z_temp_e=z_index_e%int(DEPTH_SIZE);
                }
                uint16 *block_temp = NULL;
                bool isCacheExisted;
#pragma omp critical  
                {
                    isCacheExisted = imageCaches_.isCacheExisting(str);
                    SmartCache tmpCache = imageCaches_.GetSmartCache(str);
                    tmpCache->GetImagePtr(block_temp);
                }
                
                if(isCacheExisted){
#pragma omp critical  
                    {
                        cout << str << " is read! thread:" <<omp_get_thread_num() <<" loop id:" << k << endl;
                    }
                }else{
#pragma omp critical  
                    {
                        cout << str << " is reading! thread:" << omp_get_thread_num() << " loop id:" << k << endl;
                    }
                    TIFFSetWarningHandler(0);
                    TIFF *tif = TIFFOpen(str.c_str(), "r");
                    if (tif == 0) {
                        printf(" open TIFF file %s failed\n", str.c_str());
                        //return false;
                        continue;
                    }
                    uint16 bitspersample=16;
                    TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bitspersample);
                    if(bitspersample != 16){
                        printf(" not 16bit tiff\n", str.c_str());
                        //return false;
                        continue;
                    }
                    int width(0), height(0), f(0);
                    //int wh(0);
                    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);  
                    TIFFGetField(tif, TIFFTAG_IMAGELENGTH,&height);
                    TIFFGetField(tif, TIFFTAG_COMPRESSION, &compress);
                    f = TIFFNumberOfDirectories(tif);
                    //wh = width * height;

                    //unsigned short *data = new unsigned short[width * height * f];
                    memset(data, 0, sizeof(short) * width * height);
                    TIFFSetDirectory(tif,0);
                    int StripSize = TIFFStripSize(tif);
                    //int Vstripsize = TIFFVStripSize(tif, height);
                    //int rowsperpix = TIFFComputeStrip(tif, height,0);
                    int NumberofStrips = TIFFNumberOfStrips(tif);


                    for (int i_block_temp = 0; i_block_temp<int(DEPTH_SIZE); ++i_block_temp) {
                        for (int i = 0; i < NumberofStrips; ++i)
                        {
                            TIFFReadEncodedStrip(tif, i, &(data[i * StripSize/2]), StripSize);//
                        }

                        //transform RGBA image to 8Bit
                        for (int j_block_voxel = 0; j_block_voxel<int(BLOCK_SIZE*BLOCK_SIZE); j_block_voxel++)
                        {
                            block_temp[i_block_temp*int(BLOCK_SIZE*BLOCK_SIZE) + j_block_voxel] = data[j_block_voxel];
                        }
                        TIFFReadDirectory(tif);
                    }
                    TIFFClose(tif);
                    //2017-2-25
#pragma omp critical 
                    {
                        /*if (compress != COMPRESSION_NONE) {
                            printf("NGP will transfer compressed signal image %s to uncompressed image.\n", (const char*)str.c_str());
                            std::string newName = str.substr(0, str.find_first_of('.', 0)) + std::string("_bak.tif");
                            rename(str.c_str(), newName.c_str());
                            NGImageWriter writer = ImageWriter::New();
                            writer->SetOutputFileName(str);
                            std::shared_ptr<SVolume> saveImg = std::make_shared<SVolume>();
                            saveImg->SetResolution(rez_IOR_x, rez_IOR_y, rez_IOR_z);
                            saveImg->SetSize(width, height, int(DEPTH_SIZE));
                            memcpy(saveImg->GetPointer(), block_temp, sizeof(uint16) * int(BLOCK_SIZE*BLOCK_SIZE * DEPTH_SIZE));
                            writer->SetInput(saveImg);
                            writer->Update();
                            printf("bak image:%s\n", newName.c_str());
                            }*/
                    }
                }

                //crop
                for (uint16 i_read_voxel_z=z_temp_b; i_read_voxel_z<=z_temp_e; ++i_read_voxel_z) {
                    for (uint32 j_read_voxel_y=y_temp_b;j_read_voxel_y<=y_temp_e; j_read_voxel_y++)
                    {
                        for (uint32 k_read_voxel_x=x_temp_b;k_read_voxel_x<=x_temp_e; k_read_voxel_x++)
                        {
                            buffer[((i-z_block_b)*int(DEPTH_SIZE)+(i_read_voxel_z-z_shift))*width_IOR*height_IOR+((j-y_block_b)*int(BLOCK_SIZE)+(j_read_voxel_y-y_shift))*width_IOR+(k-x_block_b)*int(BLOCK_SIZE)+k_read_voxel_x-x_shift]=block_temp[i_read_voxel_z*int(BLOCK_SIZE)*int(BLOCK_SIZE)+j_read_voxel_y*int(BLOCK_SIZE)+k_read_voxel_x];
                        }
                    }
                }
                cout << str << " read complete. thread:" << omp_get_thread_num() << endl;
                str.clear();
            }
        }
        delete[] data;
    }
        
    this->IORdata=buffer;
    //delete[] slice;
    //delete[] block_temp;
    return true;

}



bool MostImage16::loadDataInfo(const char *path)
{
    //int id_t;
    //int img_type_t;           //图像类型，区分血管数据，胞体数据，荧光数据
    //long sz0_t;            //三个方向的大小，像素
    //long sz1_t;
    //long sz2_t;
    //int level_size_t;
    //int block_num_t[30][3];
    //double rez_x_t, rez_y_t, rez_z_t; //the resolution of a image pixel along the 3 axes
    char  file_path_t[150];

    FILE *fp = fopen(path, "r");
    if (!feof(fp))
    {

        fscanf(fp, "%d\n",&id);
        fscanf(fp, "%d\n",&img_type);
        fscanf(fp, "%d\n",&sz0);
        fscanf(fp, "%d\n",&sz1);
        fscanf(fp, "%d\n",&sz2);
        fscanf(fp, "%d\n",&level_size);

        for (int i =1;i<=level_size;i++)
        {
            fscanf(fp, "%d\n",&block_num[i][0]);
            fscanf(fp, "%d\n",&block_num[i][1]);
            fscanf(fp, "%d\n",&block_num[i][2]);

        }
        fscanf(fp, "%lf\n",&rez_x);
        fscanf(fp, "%lf\n",&rez_y);
        fscanf(fp, "%lf\n",&rez_z);
        fscanf(fp, "%s\n",file_path_t);
        this->file_path = file_path_t;

    }

    fclose(fp);

    return true;

}

