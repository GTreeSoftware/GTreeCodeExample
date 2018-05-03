#include "most_data.h"


bool MostImage::select_IOR_file( int xbegin,int xend,int ybegin,int yend,int zbegin,int zend,int level )
{

    if (xend-xbegin+1>this->sz0||yend-ybegin+1>this->sz1||zend-zbegin+1>sz2)
    {
        qDebug("out of size,error");
        return false;
    }
    if (xend<xbegin||yend<ybegin||zend<zbegin)
    {
        qDebug("out of size,error");
        return false;
    }

    
    const char *path;
    QByteArray text;
    QStringList list_of_index;
    int x_index_b=xbegin;
    int x_index_e=xend;
    int y_index_b=ybegin;
    int y_index_e=yend;
    int z_index_b=zbegin;
    int z_index_e=zend;
    int x_temp_b;
    int x_temp_e;
    int y_temp_b;
    int y_temp_e;
    int z_temp_b;
    int z_temp_e;
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
    this->rez_IOR_x = std::pow(2.0, level)*rez_x;
    this->rez_IOR_y = std::pow(2.0, level)*rez_y;
    this->rez_IOR_z = std::pow(2.0, level)*rez_z;
    
    //uint8 *buffer = new uint8[width_IOR*height_IOR*depth_IOR];  //申请总的内存开销
    uint16 *buffer = new uint16[width_IOR*height_IOR*depth_IOR];  //申请总的内存开销
    uint32 *slice = new uint32[int(BLOCK_SIZE*BLOCK_SIZE)];   //申请一张的内存开销，即每个slice的内存
    //uint8 *block_temp = new uint8[int(BLOCK_SIZE*BLOCK_SIZE*DEPTH_SIZE)];

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

    QString str;
    //char temppp[6];
    str.clear();

    //int index_of_xy;
    for (int i = z_block_b; i <= z_block_e; i++)
    {
        for (int j = y_block_b; j <=y_block_e; j++)
        {
            for (int k = x_block_b; k <= x_block_e; k++)
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
                list_of_index.push_back(str);   //将每个要读取的block块保存进list
                //qDebug(str.toLatin1());
            
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
                //std::cout<<x_temp_b<<" "<<x_temp_e<<" "<<y_temp_b<<" "<<y_temp_e<<" "<<z_temp_b<<" "<<z_temp_e<<" "<<std::endl;
                QByteArray text = str.toLocal8Bit();
                path = text.data();
                //uint8 *block_temp = NULL;
                uint16 *block_temp = NULL;
                bool isCacheExisted = imageCaches_.isCacheExisting(str.toStdString());
                SmartCache tmpCache = imageCaches_.GetSmartCache(str.toStdString());
                tmpCache->GetImagePtr(block_temp);
                if(isCacheExisted){
                    cout<<str.toStdString()<<" is readed!"<<endl;
                }else{
                    cout<<str.toStdString()<<" is reading!"<<endl;
                    TIFFSetWarningHandler(0);
                    TIFF *tif = TIFFOpen(path, "r");
                    if (tif == 0) {
                        printf(" open TIFF file %s failed\n", path);
                        return false;
                    }

                    for (int i_block_temp=0; i_block_temp<int(DEPTH_SIZE); ++i_block_temp) {
                        TIFFReadRGBAImageOriented(tif, int(BLOCK_SIZE), int(BLOCK_SIZE), slice, ORIENTATION_TOPLEFT);
                        for (int  j_block_voxel=0; j_block_voxel<int(BLOCK_SIZE*BLOCK_SIZE); j_block_voxel++)
                        {
                            //block_temp[i_block_temp*int(BLOCK_SIZE*BLOCK_SIZE)+j_block_voxel] = TIFFGetR(slice[j_block_voxel]);
                            block_temp[i_block_temp*int(BLOCK_SIZE*BLOCK_SIZE)+j_block_voxel] = uint16(TIFFGetR(slice[j_block_voxel]));
                        }
                        TIFFReadDirectory(tif);
                    }

                    TIFFClose(tif);
                }

                for (uint16 i_read_voxel_z=z_temp_b; i_read_voxel_z<=z_temp_e; ++i_read_voxel_z) {
                    for (uint32 j_read_voxel_y=y_temp_b;j_read_voxel_y<=y_temp_e; j_read_voxel_y++)
                    {
                        for (uint32 k_read_voxel_x=x_temp_b;k_read_voxel_x<=x_temp_e; k_read_voxel_x++)
                        {
                            buffer[((i-z_block_b)*int(DEPTH_SIZE)+(i_read_voxel_z-z_shift))*width_IOR*height_IOR+((j-y_block_b)*int(BLOCK_SIZE)+(j_read_voxel_y-y_shift))*width_IOR+(k-x_block_b)*int(BLOCK_SIZE)+k_read_voxel_x-x_shift]=block_temp[i_read_voxel_z*int(BLOCK_SIZE)*int(BLOCK_SIZE)+j_read_voxel_y*int(BLOCK_SIZE)+k_read_voxel_x];
                        }
                    }
                }
                str.clear();
            }
        }
    }
        
    this->IORdata=buffer;
    delete[] slice;
    //delete[] block_temp;
    return true;

}



bool MostImage::loadDataInfo(const char *path)
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

