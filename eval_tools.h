/* 
* evaluation tools, mainly papi and timer
*/

#ifndef EVAL_TOOLS_H
#define EVAL_TOOLS_H

#include <sys/time.h>

double get_cur_time(){
    struct timeval tv;
    struct timezone tz;
    double cur_time;
    gettimeofday(&tv, &tz);
    cur_time= tv.tv_sec + tv.tv_usec/1000000.0;
    return cur_time;
}


#endif