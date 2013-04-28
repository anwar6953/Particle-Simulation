//
//  ColorAndVector.h

#ifndef _performance_h
#define _performance_h
//{
#include <ctime>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <time.h>
#include <math.h>
#include <list>
#include <stdlib.h>
#include <cstdlib>

#ifdef _WIN32
#include "performance.h"
#include <windows.h>
#else
#include <sys/time.h>
#endif

#ifdef OSX
#include <GLUT/glut.h>
#include <OpenGL/glu.h>
#else
#include <GL/glut.h>
#include <GL/glu.h>
#endif

#endif
using namespace std;
//}
const __int64 DELTA_EPOCH_IN_MICROSECS= 11644473600000000;


struct timezone2 
{
  __int32  tz_minuteswest; /* minutes W of Greenwich */
  bool  tz_dsttime;     /* type of dst correction */
};

struct timeval2 {
__int32    tv_sec;         /* seconds */
__int32    tv_usec;        /* microseconds */
};



int gettimeofday(struct timeval2 *tv/*in*/, struct timezone2 *tz/*in*/)
{
  FILETIME ft;
  __int64 tmpres = 0;
  TIME_ZONE_INFORMATION tz_winapi;
  int rez=0;

   ZeroMemory(&ft,sizeof(ft));
   ZeroMemory(&tz_winapi,sizeof(tz_winapi));

    GetSystemTimeAsFileTime(&ft);

    tmpres = ft.dwHighDateTime;
    tmpres <<= 32;
    tmpres |= ft.dwLowDateTime;

    /*converting file time to unix epoch*/
    tmpres /= 10;  /*convert into microseconds*/
    tmpres -= DELTA_EPOCH_IN_MICROSECS; 
    tv->tv_sec = (__int32)(tmpres*0.000001);
    tv->tv_usec =(tmpres%1000000);


    //_tzset(),don't work properly, so we use GetTimeZoneInformation
    rez=GetTimeZoneInformation(&tz_winapi);
    tz->tz_dsttime=(rez==2)?true:false;
    tz->tz_minuteswest = tz_winapi.Bias + ((rez==2)?tz_winapi.DaylightBias:0);

  return 0;
}
