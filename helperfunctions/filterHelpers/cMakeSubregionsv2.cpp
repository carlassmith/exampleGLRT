/** \file cMakeSubregions.c
*
*Compile the mex file:
*
mex -O cMakeSubregions.c
*
*      [subregions leftcoord topcoord]=cMakeSubregions(X,Y,Z,sz,data)

 * $Rev: 54 $
 * $Date: 2011-09-06 13:18:16 -0600 (Tue, 06 Sep 2011) $
 * $Author: jbyars $
 * $HeadURL: https://abbe.phys.unm.edu/svn/cCode/cMakeSubregions.c $
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mex.h"
#include <math.h>

#ifndef max
//! not defined in the C standard used by visual studio
#define max(a,b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
//! not defined in the C standard used by visual studio
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif


//#define _DEBUG

// Thread block size

void getbox(const float *data,const int ii,const int sz,int szZ,const int szinX,const int szinY,
        const int szinZ,const double *X,const double *Y,const double *Z, float* dataout,float* leftcoord,float* topcoord,float* depthcoord);
void version();

void mexFunction(int nlhs, mxArray *plhs[],	int	nrhs, const	mxArray	*prhs[])
{
	/** \brief This function takes a list of x,y,z coordinates and crops out subregions of size 
	* sz x sz at the specified locations in the dataset.  It returns the array of subregions along with a list
	* of left and top coordinates for each. 
	*  \param nlhs number of left hand mxArrays to return
	*  \param plhs array of pointers to the output mxArrays
	*  \param nrhs number of input mxArrays
	*  \param prhs array of pointers to the input mxArrays.
	*/

	int sz,szZ,szinX,szinY,szinZ;
	//int N;
	int outsize[4];
	float *dataout;				//!< an array of 2D subregions
	float *leftcoord;			//!< an array of the left coordinates of the subregions
	float *topcoord;			//!< an array of the top coordinates of the subregions
    float *depthcoord;			//!< an array of the depth coordinates of the subregions
	int ii;						//!< loop counter
	const mwSize *Xdatasize;	//!< dimensions of x coordinates
	const mwSize *Ydatasize;	//!< dimensions of y coordinates
	const mwSize *Zdatasize;	//!< dimensions of z coordinates
	const mwSize *boxsize;		//!< dimensions of box size dimensions (just for error checking)
    const mwSize *boxsizeZ;		//!< dimensions of box size dimensions (just for error checking)
	const mwSize *datasize;		//!< dimensions of data set
	const double *X;			//!< array of x coordinates
	const double *Y;			//!< array of y coordinates
	const double *Z;			//!< array of z coordinates
	const float *data;			//!< dataset to crop regions out of
	int N, NdimX, NdimY, NdimZ, NdimData;

    //version();
	//input checks
	if (nrhs<5) 
		mexErrMsgIdAndTxt("MATLAB:WrongNumberOfInputs",
		 "Usage: cMakeSubregions(xpostitions,ypositions,zpositions, boxsize,data)\n$Rev: 54 $\n$Date: 2011-09-06 13:18:16 -0600 (Tue, 06 Sep 2011) $\n$Author: jbyars $\n$HeadURL: https://abbe.phys.unm.edu/svn/cCode/cMakeSubregions.c $");
    
	Xdatasize=mxGetDimensions(prhs[0]);
	Ydatasize=mxGetDimensions(prhs[1]);
	Zdatasize=mxGetDimensions(prhs[2]);
	boxsize=mxGetDimensions(prhs[3]);
    datasize=mxGetDimensions(prhs[4]);
    if (nrhs == 5)
        szZ=-1;
    else
        boxsizeZ=mxGetDimensions(prhs[5]);
    
	NdimX=mxGetNumberOfDimensions(prhs[0]);
	NdimY=mxGetNumberOfDimensions(prhs[1]);
	NdimZ=mxGetNumberOfDimensions(prhs[2]);
	NdimData=mxGetNumberOfDimensions(prhs[4]);

	//input checks
	if (mxGetClassID(prhs[4])!=mxSINGLE_CLASS)
		mexErrMsgIdAndTxt("cMakeSubregions:WrongPrecision","data must be comprised of single floats\n");
	if (Xdatasize[0] == 0)
		mexErrMsgIdAndTxt("cMakeSubregions:WrongDimensions","X coordinate array is empty.\n");
	if (Ydatasize[0] == 0)
		mexErrMsgIdAndTxt("cMakeSubregions:WrongDimensions","Y coordinate array is empty.\n");
	if (Zdatasize[0] == 0)
		mexErrMsgIdAndTxt("cMakeSubregions:WrongDimensions","Z coordinate array is empty.\n");
	if (NdimX < 1 || NdimX > 2 || Xdatasize[1] != 1) //min 2 dim
		mexErrMsgIdAndTxt("cMakeSubregions:WrongDimensions","X should be a 1 dimensional array.\n");
	if (NdimY < 1 || NdimY > 2 || Ydatasize[1] != 1)
		mexErrMsgIdAndTxt("cMakeSubregions:WrongDimensions","Y should be a 1 dimensional array.\n");
	if (NdimZ < 1 || NdimZ > 2 || Zdatasize[1] != 1)
		mexErrMsgIdAndTxt("cMakeSubregions:WrongDimensions","Z should be a 1 dimensional array.\n");
	if (Xdatasize[0]!=Ydatasize[0] || Xdatasize[0]!=Zdatasize[0])
		mexErrMsgIdAndTxt("cMakeSubregions:DimensionMismatch","Size of X and Y and Z coordinates must match.\n");
	if (boxsize[0]!= 1 || boxsize[1] != 1)
		mexErrMsgIdAndTxt("cMakeSubregions:WrongDimensions","boxsize should be a single value specifying the szxsz region size.\n");
    if (szZ != -1 && (boxsizeZ[0]!= 1 || boxsizeZ[1] != 1))
        mexErrMsgIdAndTxt("cMakeSubregions:WrongDimensions","boxsize should be a single value specifying the szZ region size.\n");   

 	//Matlab will return 2D for NxNx1 DipImage arrays
	if (NdimData != 2 && NdimData != 3) 
		mexErrMsgIdAndTxt("cMakeSubregions:WrongDimensions","data should be a 3 dimensional array not %d.\n", NdimData);


	N=Xdatasize[0];

	//get variables

	X=(double *) mxGetData(prhs[0]);
	Y=(double *) mxGetData(prhs[1]);
	Z=(double *) mxGetData(prhs[2]);
    sz=(int)mxGetScalar(prhs[3]); //matlab-dip_image convention
	data=(float *) mxGetData(prhs[4]);
    if(szZ==-1)
        szZ=1;
    else
        szZ=(int)mxGetScalar(prhs[5]);  
    
	//validation
	szinX=datasize[0];
	szinY=datasize[1];
	szinZ=datasize[2];

    if (sz<=0)
		mexErrMsgIdAndTxt("cMakeSubregions:InputOutofRange","boxsize must be a positive integer. (%d)\n",&sz);

	if (sz>szinX)
		mexErrMsgIdAndTxt("cMakeSubregions:InputOutofRange","data size must be bigger than boxsize (%d)\n",szinX);

	if (sz>szinY)
		mexErrMsgIdAndTxt("cMakeSubregions:InputOutofRange","data size must be bigger than boxsize\n");
    if (szZ>szinZ)
		mexErrMsgIdAndTxt("cMakeSubregions:InputOutofRange","data size must be bigger than boxsize\n");
    
      if (szZ == 1){
         outsize[0]=sz;
         outsize[1]=sz;
         outsize[2]=N;
         plhs[0]=mxCreateNumericArray(3,outsize,mxSINGLE_CLASS,mxREAL);
         dataout=(float *)mxGetData(plhs[0]);   
     }else{        
        outsize[0]=sz;
        outsize[1]=sz;
        outsize[2]=szZ;
        outsize[3]=N;
        plhs[0]=mxCreateNumericArray(4,outsize,mxSINGLE_CLASS,mxREAL);
        dataout=(float *)mxGetData(plhs[0]);
     }
    
    mexPrintf("%dx%d",sz,szZ);
	outsize[0]=N;
	outsize[1]=1;
	outsize[2]=1;
	plhs[1]=mxCreateNumericArray(2,outsize,mxSINGLE_CLASS,mxREAL);
	plhs[2]=mxCreateNumericArray(2,outsize,mxSINGLE_CLASS,mxREAL);
    plhs[3]=mxCreateNumericArray(2,outsize,mxSINGLE_CLASS,mxREAL);
	leftcoord=(float *)mxGetData(plhs[1]);
	topcoord=(float *)mxGetData(plhs[2]);
    depthcoord=(float *)mxGetData(plhs[3]);
	
    for (ii=0;ii<N;ii++){
		getbox(data,ii,sz,szZ,szinX,szinY,szinZ,X,Y,Z,dataout,leftcoord,topcoord,depthcoord);
    }

	return;
};

void getbox(const float *data,const int ii,const int sz,int szZ,const int szinX,const int szinY,
        const int szinZ,const double *X,const double *Y,const double *Z, float* dataout,float* leftcoord,float* topcoord,float* depthcoord)
{
	/** \brief This function copies the specified subregion in data to data 
	*  \param data the original data set to crop from
	*  \param ii the point to copy
	*  \param sz the size of the subregion to copy
	*  \param szinX x dimension of data
	*  \param szinY y dimension of data
	*  \param szinZ z dimension of data
	*  \param X x coordinate of the center of the subregion
	*  \param Y y coordinate of the center of the subregion
	*  \param Z z coordinate of the center of the subregion
	*  \param dataout array of subregions to copy to
	*  \param leftcoord left coordinate of the subregion in the original image
	*  \param topcoord right coordinate of the subregion in the original image
	*/

	int yy,zz;
	int l,r,t,b,q,o;
	
	const int szl=(int) floor(sz/2.0+.5);
	//get coordinates
	const int x=(int) floor(X[ii]+.5)+1; 
	const int y=(int) floor(Y[ii]+.5)+1;
	const int z=(int) floor(Z[ii]+.5)+1;
    
 #ifdef _DEBUG
	// x and y get +1 so they can = szin
	if (x<1 || y<1 || z<0 || x>szinX || y>szinY || z>szinZ) {
		static char msgbuf[256];
		sprintf(msgbuf,"Point %d out of bounds position %d,%d,%d dataset size %d,%d,%d",ii,x,y,z,szinX, szinY, szinZ);
		mexErrMsgIdAndTxt("cMakeSubregions:InputOutofRange",msgbuf);
	}
    mexPrintf("ii: %d %d %d %d %d\n",ii,x,y,z,szl);
 #endif



	//if (z>szinZ-1)
	//{
	//	mexPrintf("%d %d %d %d\n",x,y,z,szl);
	//    mexErrMsgTxt(" Z index larger than datasize\n");
	//}

	l=max(x-szl,0);
	r=l+sz-1;
	if (r>(szinX-1)) {r=szinX-1;l=r-sz+1;}
	t=max(y-szl,0);
	b=t+sz-1;
	if (b>(szinY-1)) {b=szinY-1;t=b-sz+1;}

  	q=max(z-floor(szZ/2.0+.5),0);
	o=q+szZ-1;
	if (o>(szinZ-1)) {o=szinZ-1;q=o-szZ+1;}
    
    for (zz=0;zz<szZ;zz++) {
        for (yy=0;yy<sz;yy++){
            #ifdef _DEBUG
             mexPrintf("%d:%d\n",(l+0) + szinX * ((t+yy)+ szinY * (q+zz)),(int)data[(l+0) + szinX * ((t+yy)+ szinY * (q+zz))]);
            #endif                      
             
            memcpy(dataout+(sz*sz*szZ*ii+sz*(yy+sz*zz)+0),data+((l+0) + szinX * ((t+yy)+ szinY * (q+zz))),sz*sizeof(float));
            
            #ifdef _DEBUG
             mexPrintf("%d:%d\n",sz*sz*szZ*ii+sz*(yy+sz*zz)+0,(int)dataout[sz*sz*szZ*ii+sz*(yy+sz*zz)+0]);
            #endif           
        }
	}
	leftcoord[ii]=(float) l;
	topcoord[ii]=(float) t;
    depthcoord[ii]=(float)q;
    
	return;
}

/*
void version()
{
// 
//  \brief Add the svn version info of this mex to the Matlab global VERSION.
//  I.E. if the mexFunction is listGPUs, you would have VERSION.listGPUs.Rev, etc 
//  after running the mexfile.  
 
    mxArray *ver, *info;
    int field = 0;
    const char *field_names[] = {"Rev","Date","Author","HeadURL"};
    const int NUMBER_OF_FIELDS = 4;
	char svnrev[] = "$Rev: 54 $";
    char svndate[] = "$Date: 2011-09-06 13:18:16 -0600 (Tue, 06 Sep 2011) $";
    char svnauthor[] = "$Author: jbyars $";
    char svnheadurl[] = "$HeadURL: https://abbe.phys.unm.edu/svn/cCode/cMakeSubregions.c $";
	char path[] = __FILE__;
    char *rev,*date,*author,*headurl,*file,*eos;
    mwSize dims[2] = {1, 1};
	// REALLY UGLY string mangling
	rev = strchr(svnrev,':');
	if (rev == NULL) mexErrMsgTxt("VERSION: Rev is empty");
	rev+=2; //skip to revision
	rev[strlen(rev)-2] = '\0';
	date = strchr(svndate,':');
	if (date == NULL) mexErrMsgTxt("VERSION: Date is empty");
	date+=2;
	date[strlen(date)-2] = '\0';
	author = strchr(svnauthor,':');
	if (author == NULL) mexErrMsgTxt("VERSION: author is empty");
	author+=2;
	author[strlen(author)-2] = '\0';
	headurl = strchr(svnheadurl,':');
	if (headurl == NULL) mexErrMsgTxt("VERSION: headurl is empty");
	headurl+=2;
	headurl[strlen(headurl)-2] = '\0';
	file = strrchr(path,'\\');
	if (file == NULL) 
		file = strrchr(path,'/');
	if (file == NULL) 		
		mexErrMsgTxt("VERSION: filename is empty");
	file++; // move past the slash
	eos = strchr(file,'.');
	*eos = '\0';
	
    //create the struct of version info
    info = mxCreateStructArray(2, dims, NUMBER_OF_FIELDS, field_names); 
    if (info == NULL) mexErrMsgTxt("Could not create info struct for versioning");
    mxSetFieldByNumber(info,0,0,mxCreateString(rev));
    mxSetFieldByNumber(info,0,1,mxCreateString(date));
    mxSetFieldByNumber(info,0,2,mxCreateString(author));
    mxSetFieldByNumber(info,0,3,mxCreateString(headurl));
    
    // see if the global exists
    ver = mexGetVariable("global","VERSION");
    if (ver == NULL) {  //if not create the global VERSION
        ver = mxCreateStructArray(2, dims, 1, &file);
    }
    if (ver == NULL) mexErrMsgTxt("Could not access or create VERSION global");
    // if the field for this mex doesn't exist, add it
    field = mxGetFieldNumber(ver,file);
    if (field <0)
        field = mxAddField(ver,file);
    if (field <0) mexErrMsgTxt("Could not create a field in VERSION");
    mxSetField(ver,0,file,info);
    mexPutVariable("global","VERSION", ver);
}

*/



