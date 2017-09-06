#ifndef FILTER_H
#define FILTER_H

void unifilter(const float * data , float * result, const int fsz , const int xsz, const int ysz, const int tsz); 
void maxfilter(const float * data , float * result, const int fsz , const int xsz, const int ysz, const int tsz);

#endif