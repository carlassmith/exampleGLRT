/*! \file cMakeNdSubregions.c
 * \brief make subregions for N-dimensional data set
 *
 * [subregions subSampling]=cMakeNdSubregions(position,boxSize,data,sampling,verbose)
 *
 * \param position center positions (zero based) of subregions to make from data. This is a
 * 2 dimensional array with dimensions # of subregions by
 * N-dimensions. DOUBLE INPUT.
 * \param boxSize vector giving size of each dimension for subregions. Size of
 * vector is 1 by N-dimensions. DOUBLE INPUT
 * \param data N-dimensional array of data. SINGLE INPUT.
 * \param sampling cell array with dimensions N by 1. Each cell contains a
 * vector with length of dimension size(data,ii) by 1 shere the vector
 * describes subregions sampling in demension ii. DOUBLE INPUT.
 *
 * \return dataOut (N+1)-dimensional array of subregions. The added dimension
 * contains each subregion. Example for a 3-D array
 * dataOut(:,:,:,ii) selects the ith subregion. SINGLE OUTPUT.
 * \return subSampling cell array with dimensions N by 1. Each cell contains a
 * 2D array with length of dimension size(dataOut,ii) by # of
 * subregions describing sampling in demension ii. DOUBLE OUTPUT.
 *
 *
 * NOTE: convention used for all inputs is standard matlab notation for
 * indexing. Dimension order is rows, columns, etc.
 *
 * See 'cMakeNdSubregionsTesting.m' for sample usage and testing
 *
 * \Author Pat Cutler
 * \Date September 2010 (UNM)
 *
 * Compile in matlab with 'mex -O cMakeNdSubregions.c'
 *
 */

#include <stdlib.h>
#include <string.h>
#include "mex.h"
#include <math.h>

#define max(a,b) ( (a) >= (b) ? (a) : (b) )  


// Thread block size

void getbox(const float*,int,const double*,int,const int*,const double*,const int*,float*,float*,const double**,double**);

void mexFunction(int nlhs, mxArray *plhs[],	int	nrhs, const	mxArray	*prhs[])
{
    const int *positionsSize=0;
    const int *boxSizeSize=0;
    const int *dataSize=0;
    const int *samplingSize=0;
    const int *tempSize=0;
    int ii, *outSizeShift, *outSizeSubRegions, positionsNd, boxSizeNd, dataNd, samplingNd;
    const double *boxSize=0;
    const double *positions=0;
    const double *boxsize=0;
    const float *data=0;
    float *dataOut, *shift;
    const double **sampling=0;
    double **subSampling=0;
    const mxArray *samplingCell=0;
    mxArray *subSamplingCell=0;
    char text[1024];
    
    //input checks
    if (nrhs!=4)
        mexErrMsgTxt("cMakeNdSubregions:  4 inputs required - 'positions', 'boxSize', 'data', & sampling\n");
    
    if (mxIsCell(prhs[3]) == 0)
        mexErrMsgTxt("cMakeNdSubregions: 'sampling' input must be a cell\n");
    
    //get dimension and size information for each input
    positionsSize=(int*) mxGetDimensions(prhs[0]);
    boxSizeSize=(int*) mxGetDimensions(prhs[1]);
    dataSize=(int*) mxGetDimensions(prhs[2]);
    samplingSize=(int*) mxGetDimensions(prhs[3]);
    positionsNd=(int) mxGetNumberOfDimensions(prhs[0]);
    boxSizeNd=(int) mxGetNumberOfDimensions(prhs[1]);
    dataNd=(int) mxGetNumberOfDimensions(prhs[2]);
    samplingNd=(int) mxGetNumberOfDimensions(prhs[3]);
    
    //mexPrintf("dataNd %i\n", dataNd);
    
    //input checks
    if (positionsNd != 2 || dataNd != positionsSize[1]) {
        sprintf(text, "cMakeNdSubregions: input 'positions' must be # subregions by %i", dataNd);
        mexErrMsgTxt(text);
    }
    if (boxSizeNd != 2 || dataNd != boxSizeSize[1]) {
        sprintf(text,"cMakeNdSubregions:  input 'boxSize' must be 1 by %i", dataNd);
        mexErrMsgTxt(text);
    }
    if (samplingNd != 2 || dataNd != samplingSize[1]) {
        sprintf(text,"cMakeNdSubregions:  input 'boxSize' must be 1 by %i", dataNd);
        mexErrMsgTxt(text);
    }

    //mexPrintf("cMakeNdSubregions: # of subregions - %d\n",positionsSize[0]);
    
    //get variables
    
    positions=(double *) mxGetPr(prhs[0]);
    boxSize=(double*) mxGetPr(prhs[1]); //matlab-dip_image convention
    data=(float *) mxGetPr(prhs[2]);
    samplingCell = prhs[3];
    sampling =(const double**) mxMalloc(sizeof(double*)*dataNd); //allocate pointers to point to memory in cell array
    subSampling =(double**) mxMalloc(sizeof(double*)*dataNd); //allocate pointers to point to memory in cell array
    
    //mexPrintf("line 105\n");
    plhs[1] = mxCreateCellArray(samplingNd, samplingSize);
    subSamplingCell = plhs[1];
    //mexPrintf("line 109\n");
    for (ii = 0;ii < dataNd;ii++) {
        tempSize = (int*)mxGetDimensions(mxGetCell(prhs[3],ii));
        //mexPrintf("tempSize %i %i\n", tempSize[0], tempSize[1]);
    }
    for (ii = 0;ii < dataNd;ii++) {
        //mexPrintf("Class Name: %s\n", mxGetClassName(mxGetCell(samplingCell,ii)));
        if (mxGetClassID(mxGetCell(samplingCell,ii))!=mxDOUBLE_CLASS)
            mexErrMsgTxt("cMakeNdSubregions: 'sampling' cell must be contain of double floats\n");
        tempSize = (int*)mxGetDimensions(mxGetCell(samplingCell,ii));
        //mexPrintf("tempSize %i %i\n",tempSize[0],tempSize[1]);
        if (tempSize[0] != dataSize[ii] || tempSize[1] != 1) {
            sprintf(text, "cMakeNdSubregions: input 'sampling' cell %i must be have dimensions %i by 1 instead of %i by %i", ii, dataSize[ii],tempSize[0],tempSize[1]);
            mexErrMsgTxt(text);
        }
        //mexPrintf("line 113\n");
        sampling[ii] =(double*) mxGetPr(mxGetCell(samplingCell,ii)); //get pointer to sampling
        //sampling[ii] =(double*) mxGetPr(mxGetData(samplingCell[ii]));
        mxSetCell(subSamplingCell, ii, mxCreateNumericMatrix(boxSize[ii], positionsSize[0], mxDOUBLE_CLASS, mxREAL));//allocate memmory for subSampling cells
        subSampling[ii] =(double*) mxGetPr(mxGetCell(subSamplingCell,ii)); //get pointer to subSampling
        //subSampling[ii] =(double*) mxGetPr(mxGetData(subSamplingCell[ii]));
    }
    //mexPrintf("line 122\n");

    if (mxGetClassID(prhs[2])!=mxSINGLE_CLASS)
        mexErrMsgTxt("cMakeNdSubregions: data must be comprised of single floats\n");
    
    outSizeSubRegions = mxMalloc(sizeof(int)*(dataNd+1)); //allocate memory for size of subregion output
    
    outSizeShift = mxMalloc(sizeof(int*)*2);
    outSizeShift[0] = positionsSize[0]; //define size for shift output
    outSizeShift[1] = dataNd; //define size for shift output
    outSizeSubRegions[dataNd] = positionsSize[0]; //define size for subregions output
    for (ii = 0;ii<dataNd;ii++) {
        if (boxSize[ii]>dataSize[ii]) {
            sprintf(text,"cMakeNdSubregions: size of data dimension %i is %i and must be bigger than boxSize[%i] which is %i\n",ii,dataSize[ii],ii,boxSize[ii]);
            mexErrMsgTxt(text);
        }
        outSizeSubRegions[ii]=boxSize[ii]; //define size for subregions output
        //mexPrintf("outSizeSubRegions[%i] %i\n",ii,outSizeSubRegions[ii]);
    }
    //mexPrintf("outSizeSubRegions[%i] %i\n", dataNd, outSizeSubRegions[dataNd]);
    //mexPrintf("outSizeShift[0] %i\n", outSizeShift[0]);
    //mexPrintf("outSizeShift[1] %i\n", outSizeShift[1]);
    
    plhs[0]=mxCreateNumericArray((mwSize) dataNd+1,outSizeSubRegions,mxSINGLE_CLASS,mxREAL); //allocate memory for subregions output
    dataOut=(float *)mxGetData(plhs[0]);
    //plhs[1]=mxCreateNumericArray(2,outSizeShift,mxSINGLE_CLASS,mxREAL); //allocate memory for shift ouput
    //shift=(float *)mxGetData(plhs[1]);
    shift=(float *)mxGetData(mxCreateNumericArray(2, outSizeShift, mxSINGLE_CLASS, mxREAL)); //allocate memory for shift ouput
    
    //mexPrintf("line159");
    
    for (ii=0;ii<positionsSize[0];ii++) //main loop to get subregions
        getbox(data,ii,boxSize,dataNd,dataSize,positions,positionsSize,dataOut,shift,sampling,subSampling);
  
   return;
};


void determineNDlocation(int element, int numberOfElements, int numberOfDimensions,const double* dimensions, int* NdLocation) {
    /*! \breif determine location of element in relation to dimensions
     * \param element index of element
     * \param numberOfElements total number of elements in data set
     * \param numberOfDimensions number of Dimensions in data set
     * \param dimensions pointer to size of each dimension
     * \param NdLocation pointer to location in each dimension that is determined
     */
    int ii;
    for (ii = numberOfDimensions-1;ii>-1;ii--) {
        numberOfElements = numberOfElements/(int)dimensions[ii];
        NdLocation[ii] = floor((float) (element/numberOfElements));
        element = element-NdLocation[ii]*numberOfElements;
        //mexPrintf("dimension %i NdLocation %i\n",ii,NdLocation[ii]);
    }
}



void getbox(const float* data,int ii,const double* boxSize,int dataNd,const int* dataSize,const double* positions,const int* positionsSize,float* dataOut,float* shift,const double** sampling,double** subSampling)
{
    /*! \breif get specific subregion
     * \param data pointer to data
     * \param ii subregion number
     * \param boxSize size of each dimension for subregions.
     * \param dataNd number of dimensions in data
     * \param dataSize pointer to dimension sizes of data
     * \param positions pointer to positions input
     * \param positionsSize pointer dimension size positions
     * \param dataOut pointer to output data
     * \param shift pointer to shift in each dimension for subregion
     * \param sampling pointer pointer to sampling input
     * \param subSampling pointer pointer to subSampling output
     * \
     */
    int jj,kk,count, center,*numberOfSubRegionElements;
    int *NdLocation,*NdLocation1, *numberOfNdElements,idxIn1,idxOut1,idxIn,idxOut;
    double *tempSize;
    
    tempSize = (double*) mxMalloc(sizeof(double)*(dataNd)); //allocate memory for total number of subregion elements
    numberOfNdElements = (int*) mxMalloc(sizeof(int)*(dataNd)); //allocate memory for total number elements in an N dimensional slice
    numberOfSubRegionElements = (int*) mxMalloc(sizeof(int*)*(dataNd)); //allocate memory for total number of subregion elements
    NdLocation = (int*) mxMalloc(sizeof(int*)*(dataNd)); //allocate memory for total number of subregion elements
    NdLocation1 = (int*) mxMalloc(sizeof(int*)*(dataNd)); //allocate memory for total number of subregion elements
    numberOfNdElements[0] = 1; //initialize variable for total number elements in an N dimensional slice
    idxIn1= 0; //initialize variable for index of first element in subregion for input data
    //mexPrintf("position %i: ", ii);
    for (jj = 0;jj < dataNd;jj++) {
        //mexPrintf("(positions idx %i) ", jj*positionsSize[0] + ii);
        center = (int) floor(positions[jj*positionsSize[0] + ii]+0.5)+1; //round and change to 1 based index??
        //mexPrintf("%i  ", center);
        //mexPrintf("(%i-floor(boxSize[%i]/2.0+0.5) %f) ",center,jj,center-floor(boxSize[jj]/2.0+0.5));
        shift[jj*positionsSize[0] + ii] = (float) max(center-floor(boxSize[jj]/2.0+0.5), 0);
        if (shift[jj*positionsSize[0] + ii]+boxSize[jj] > dataSize[jj]-1)
            shift[jj*positionsSize[0] + ii] = max(dataSize[jj]-boxSize[jj], 0);
        if (jj == 0) {
            idxIn1 += ((int)shift[jj*positionsSize[0] + ii])*numberOfNdElements[jj];
            numberOfNdElements[jj] = (int) dataSize[jj];
            numberOfSubRegionElements[jj] = (int) boxSize[jj];
        }
        else {
            idxIn1 += (int) (shift[jj*positionsSize[0] + ii])*numberOfNdElements[jj-1];
            numberOfNdElements[jj] = (int) numberOfNdElements[jj-1]*dataSize[jj];
            numberOfSubRegionElements[jj] = (int) numberOfSubRegionElements[jj-1]*boxSize[jj];
        }
        tempSize[jj] = (double) dataSize[jj];
        //mexPrintf("(shift[%i] %1.0f) ", jj*positionsSize[0] + ii, shift[jj*positionsSize[0] + ii]);
        //mexPrintf("(idxIn1 %i) ", idxIn1);
        //mexPrintf("(numberOfNdElements[%i] %i) ",jj, numberOfNdElements[jj]);
        //mexPrintf("(numberOfSubRegionElements[%i] %i) ",jj, numberOfSubRegionElements[jj]);
    }
    //mexPrintf("\n");
    
    //for (jj = 0;jj < dataNd;jj++)
    //    mexPrintf("numberOfNdElements[%i] %i\n",jj, numberOfNdElements[jj]);
    //for (jj = 0;jj < dataNd;jj++)
   //     mexPrintf("numberOfSubRegionElements[%i] %i\n",jj, numberOfSubRegionElements[jj]);
    
    idxOut1= ii*numberOfSubRegionElements[dataNd-1]; //initialize variable for index of first element in subregion for output data
    //mexPrintf("1st element idxIn %i\n", idxIn1);
    //mexPrintf("1st element idxOut %i\n", idxOut1);
    for (jj = 0;jj < numberOfSubRegionElements[dataNd-1];jj++) {
        idxOut = idxOut1+jj; //get index of current element for dataOut
        //mexPrintf("jj %i----------\n", jj);
        //mexPrintf("idxOut %i\n", idxOut);
        determineNDlocation(jj,numberOfSubRegionElements[dataNd-1],dataNd, boxSize, NdLocation);
        idxIn = idxIn1+NdLocation[0];
        for (kk = 1; kk < dataNd;kk++) {
            //subSampling[kk][ii*dataNd+NdLocation[kk]] = subSampling[kk][idxIn1+jj];
            idxIn += (int)NdLocation[kk]*numberOfNdElements[kk-1];
            //mexPrintf("numberOfNdElements[%i] %i\n",kk-1, numberOfNdElements[kk-1]);
            //mexPrintf("NdLocation[%i] %i\n",kk, NdLocation[kk]);
            //mexPrintf("idxIn %i\n", idxIn);
        }
        determineNDlocation(idxIn, numberOfNdElements[dataNd-1], dataNd, tempSize, NdLocation1);
        //mexPrintf("jj %i kk %i NdLocation %i NdLocation1 %i\n", jj, kk, NdLocation[1], NdLocation1[1]);
        for (kk = 0; kk < dataNd;kk++) {
            //mexPrintf("NdLocation1[%i] %i\n",kk, NdLocation1[kk]);
            //subSampling[kk][ii*dataNd+NdLocation[kk]] = subSampling[kk][NdLocation[kk]];
            subSampling[kk][(int)(ii*boxSize[kk]+NdLocation[kk])] = sampling[kk][NdLocation1[kk]];
        }
        //mexPrintf("idxIn %i\n", idxIn);
        //mexPrintf("idxOut %i\n", idxOut);
        //mexPrintf("data %f\n", data[idxIn]);
        //mexPrintf("dataOut %f\n", dataOut[idxOut]);
        dataOut[idxOut] = data[idxIn];
    }

    return;
}



