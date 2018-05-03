#ifndef VOLUMECREATOR_H
#define VOLUMECREATOR_H
#include <string.h>
template<class T>
void Create2DArray(T **ptr, const int z, const int xy){
    if (xy < 1 || z < 1)
        return;
    if (*ptr){
        printf("error!!\n");
    }
    else{
        *ptr = new T[xy * z];
        memset(*ptr, 0, sizeof(T) * xy*z);
    }
}

template<class T>
void Set2DArray(T *ptr, const int z, const int xy, const int val){
    if (ptr){
        memset(ptr, val, sizeof(T) * xy*z);
    }
    else{
        printf("NULL Volume!!\n");
    }
}

template<class T>
void Destroy2DArray(T **ptr, const int z, const int xy){
    if (xy < 1 || z < 1) return;
    if (*ptr != NULL){
        delete *ptr;
    }
    *ptr = NULL;/*warning!!!pointer is dangerous*/
}

//template<class T>
//void Create2DArray(T ***ptr, const int z, const int xy){
//    if (xy < 1 || z < 1)
//        return;
//    if (*ptr){
//        printf("error!!\n");
//    }
//    else{
//        *ptr = new T*[z];
//        for (int i = 0; i < z; ++i){
//            (*ptr)[i] = new T[xy];
//            memset((*ptr)[i], 0, sizeof(T) * xy);
//        }
//    }
//}
//
//template<class T>
//void Set2DArray(T **ptr, const int z, const int xy, const int val){
//    if (ptr){
//        for (int i = 0; i < z; ++i){
//            memset(ptr[i], val, sizeof(T) * xy);
//        }
//    }
//    else{
//        printf("NULL Volume!!\n");
//    }
//}


//template<class T>
//void Destroy2DArray(T ***ptr, const int z, const int xy){
//    if (xy < 1 || z < 1) return;
//    if (*ptr != NULL){
//        if(z > 1){
//            for (int i = 0; i < z; ++i){
//                if(xy > 1)
//                    delete[] (*ptr)[i];
//                else
//                    delete (*ptr)[i];
//            }
//            delete[] *ptr;
//        }
//        else{
//            if(xy > 1)
//                delete[] (*ptr)[0];
//            else
//                delete (*ptr)[0];
//            delete *ptr;
//        }
//    }
//    *ptr = NULL;/*warning!!!pointer is dangerous*/
//}

#endif // VOLUMECREATOR_H
