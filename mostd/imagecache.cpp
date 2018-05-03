#include "imagecache.h"


ImageCache::ImageCache(size_t bitsNumArg):
imageBits_(bitsNumArg)
{
    if(imageBits_ == 8){
        imageData8_ = new unsigned char[512 * 512 * 512];
        imageData16_ = NULL;
    }else if(imageBits_ == 16){
        imageData16_ = new unsigned short[512 * 512 * 512];
        imageData8_ = NULL;
    }else{
        imageData16_ = NULL;
        imageData8_ = NULL;
    }
    lastUsedTime_ = time(0);
    //long long int cTime;
    //QueryPerformanceCounter(&cTime);
    //buildTime_ = cTime.QuadPart;
    //lastUsedTime_ = cTime.QuadPart;
}

ImageCache::~ImageCache()
{
    if(imageData8_){
        delete []imageData8_;
    }
    if(imageData16_){
        delete []imageData16_;
    }
}

void ImageCache::Modified()
{
    lastUsedTime_ = time(0);
    //printf("%ld", lastUsedTime_);
    //LARGE_INTEGER cTime;
    //QueryPerformanceCounter(&cTime);
    //lastUsedTime_ = cTime.QuadPart;
}

ImageCaches::ImageCaches( size_t numArg, size_t bitsArg ):
cacheMaxNum_(numArg),
imageBits_(bitsArg)
{

}

ImageCaches::~ImageCaches()
{

}

bool ImageCaches::CachSortFun( const SmartCache &arg1, const SmartCache &arg2 )
{
    return arg1->GetLastUsedTime() > arg2->GetLastUsedTime();
}

bool ImageCaches::isCacheExisting( const string &nameArg )
{
    if(cachePtrs_.size() != 0){
        for(deque<SmartCache>::size_type i = 0; i < cachePtrs_.size(); ++i){
            SmartCache cachePtrI = cachePtrs_.at(i);
            if(nameArg.compare(cachePtrI->GetImageName()) == 0){
                return true;
            }
        }
    }
    return false;
}

void ImageCaches::SortCaches()
{
    sort(cachePtrs_.begin(), cachePtrs_.end(), CachSortFun);
    /*printf("sortcaches: ");
    for (auto it = cachePtrs_.begin(); it != cachePtrs_.end(); ++it) {
    printf("%ld ", (*it)->GetLastUsedTime());
    }
    printf("\n");*/
}

void ImageCaches::DeleteUnusedCache()
{
    if (!cachePtrs_.empty()) {
        cachePtrs_.pop_back();
    }
}

void ImageCaches::AddCache( const string &nameArg )
{
    SmartCache tmpCache(new ImageCache(imageBits_));
    tmpCache->SetImageName(nameArg);
    if(cachePtrs_.size() >= cacheMaxNum_ && cacheMaxNum_ >= 2){
        SortCaches();
        DeleteUnusedCache();
    }else if(cacheMaxNum_ == 1){
        DeleteUnusedCache();
    }
    cachePtrs_.push_front(tmpCache);
}

SmartCache  ImageCaches::GetSmartCache( const string &nameArg )
{
    if(cachePtrs_.size() != 0){
        for(deque<SmartCache>::size_type i = 0; i < cachePtrs_.size(); ++i){
            SmartCache cachePtrI = cachePtrs_.at(i);
            if(nameArg.compare(cachePtrI->GetImageName()) == 0){
                cachePtrI->Modified();
                return cachePtrI;
            }
        }
        AddCache(nameArg);
    }else{
        AddCache(nameArg);
    }
    return *cachePtrs_.begin();
}

void ImageCaches::TravelCaches()
{
    cout<<"-------------------------------"<<endl;
    for(deque<ImageCache>::size_type i = 0; i < cachePtrs_.size(); ++i){
        SmartCache tmpCache = cachePtrs_.at(i);
        cout<<i<<endl;
        cout<<tmpCache->GetImageName()<<endl;
        cout<<tmpCache->GetLastUsedTime()<<endl;
    }
    cout<<"##########################"<<endl;
}

void ImageCaches::ClearCache()
{
    cachePtrs_.clear();
}