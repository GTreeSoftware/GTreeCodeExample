#ifndef IMAGECACHE_H
#define IMAGECACHE_H
#include <sstream>
#include <iostream>
#include <math.h>
#include <memory>
#include <algorithm>
#include <list>
#include <deque>
#include <vector>
#include <time.h>
using namespace std;
class ImageCache;
#ifdef _WIN32
typedef std::tr1::shared_ptr<ImageCache> SmartCache;
#else
typedef std::shared_ptr<ImageCache> SmartCache;
#endif

class ImageCache{
public:
    ImageCache(size_t bitsNumArg);
    ~ImageCache();
    void Modified();
    void GetImagePtr(unsigned char* &arg){arg = imageData8_;}
    void GetImagePtr(unsigned short* &arg){arg = imageData16_;}
    void SetImageName(const string &arg){imageName_ = arg;}
    size_t GetImageBits(){return imageBits_;}
    time_t GetLastUsedTime(){return lastUsedTime_;}
    string GetImageName(){return imageName_;}
    //LONGLONG GetBuildTime(){return buildTime_;}

private:
    size_t imageBits_;
    unsigned char *imageData8_;
    unsigned short *imageData16_;
    time_t lastUsedTime_;
    //LONGLONG buildTime_;
    string imageName_;
};


class ImageCaches{
public:
    ImageCaches(size_t numArg, size_t bitsArg);
    ~ImageCaches();
    SmartCache GetSmartCache(const string &nameArg);
    void TravelCaches();
    bool isCacheExisting(const string &nameArg);  
    void ClearCache();

private:
    static bool CachSortFun(const SmartCache &arg1, const SmartCache &arg2);
    void SortCaches();
    void DeleteUnusedCache();
    void AddCache(const string &nameArg);
    
    std::deque<SmartCache> cachePtrs_;
    size_t cacheMaxNum_;
    size_t imageBits_;
};

#endif