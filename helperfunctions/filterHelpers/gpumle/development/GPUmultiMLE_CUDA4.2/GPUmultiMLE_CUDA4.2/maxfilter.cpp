/*! \file maxfilter.cu
 *  \author Fang Huang
 *  \date October 10, 2010
 *  \brief This code provides a set of functions to do a maximum filter.
 */

#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include "definitions.h"
#include "filter.h"

void maxfilter(const float * data , float * result, const int fsz , const int xsz, const int ysz, const int tsz)
{
/*! \brief Perform a maximum filter on the original dataset.
 *  \param data the original dataset
 *  \param result the maximum filtered output
 *  \param fsz filter size fsz x fsz
 *  \param xsz data's x dimension size
 *  \param ysz data's y dimension size
 *  \param tsz data's time dimension size
 */

	int ii,jj,tt,ll;
	int x1,x2;
	float *temp,maxfp;
	int step=(int)floor((float)(fsz-1)/2);
	temp=(float*) malloc (tsz*xsz*ysz*sizeof(float));
	memset(temp,0,tsz*xsz*ysz*sizeof(float));
	for(tt=0;tt<tsz;tt++)
	{
		for(jj=0;jj<ysz;jj++)
		{
			for(ii=0;ii<xsz;ii++)
			{
				x1=max(0,ii-step);
				x2=min(xsz-1,ii+step);
				maxfp=0;
				for(ll=x1;ll<=x2;ll++)
					maxfp=max(maxfp,data[tt*xsz*ysz+jj*xsz+ll]);

				temp[tt*xsz*ysz+jj*xsz+ii]=maxfp;
			}
		}

		for(ii=0;ii<xsz;ii++)
		{
			for(jj=0;jj<ysz;jj++)
			{

				x1=max(0,jj-step);
				x2=min(ysz-1,jj+step);
				maxfp=0;
				for(ll=x1;ll<=x2;ll++)
					maxfp=max(maxfp,temp[tt*xsz*ysz+ll*xsz+ii]);

				result[tt*xsz*ysz+jj*xsz+ii]=maxfp;
			}
		}
	}

	free (temp);
	return;
}