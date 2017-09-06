/*! \file unifilter.cu
 *  \author Fang Huang
 *  \date October 10, 2010
 *  \brief This code provides a set of functions that perform a uniform filter.
 */

#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include "definitions.h"
#include "filter.h"

void unifilter(const float * data , float * result, const int fsz , const int xsz, const int ysz, const int tsz) 
{
/*! \brief Perform a uniform filter on the original dataset.
 *  \param data the original dataset
 *  \param result the maximum filtered output
 *  \param fsz ???
 *  \param xsz data's x dimension size
 *  \param ysz data's y dimension size
 *  \param tsz data's time dimension size
 */
	int ii,jj,tt,ll;
	int x1, x2;
	float *temp,unifsum;
	int step=(int)floor((float)(fsz-1)/2+0.5f);
	const size_t temp_size = tsz*xsz*ysz*sizeof(float);
	temp=(float*) malloc (temp_size);
	memset(temp,0,temp_size);
	for(tt=0;tt<tsz;tt++)
	{
		for(jj=0;jj<ysz;jj++)
		{
			for(ii=0;ii<xsz;ii++)
			{
				x1=max(0,ii-step);
				x2=min(xsz-1,ii+step);
				unifsum=0;
				for(ll=x1;ll<=x2;ll++)
					unifsum+=data[tt*xsz*ysz+jj*xsz+ll];

				temp[tt*xsz*ysz+jj*xsz+ii]=unifsum/(x2-x1+1);
			}
		}

		for(ii=0;ii<xsz;ii++)
		{
			for(jj=0;jj<ysz;jj++)
			{

				x1=max(0,jj-step);
				x2=min(ysz-1,jj+step);
				unifsum=0;
				for(ll=x1;ll<=x2;ll++)
					unifsum+=temp[tt*xsz*ysz+ll*xsz+ii];

				result[tt*xsz*ysz+jj*xsz+ii]=unifsum/(x2-x1+1);
			}
		}
	}

	free (temp);

	return;
}