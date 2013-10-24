/*
 * FDL_boxcountdim.h
 *
 *  Created on: Jan 31, 2011
 *      Author: abon
 */

#ifndef FDL_BOXCOUNTDIM_H_
#define FDL_BOXCOUNTDIM_H_

#include "libht.h"
#include "FDL_util.h"
#include "FDL_streamline.h"

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

using namespace cv;


class Comparator2
{
public:
	bool operator()( VECTOR3 v1, VECTOR3 v2 )
	{
		if( v1[2] < v2[2] ) return true;
		if( v1[2] > v2[2] ) return false;
		if( v1[1] < v2[1] ) return true;
		if( v1[1] > v2[1] ) return false;
		if( v1[0] < v2[0] ) return true;
		if( v1[0] >= v2[0] ) return false;
	}// end operator definition
};

class FDL_boxcountdim
{
public:

	VECTOR3 lowGrid;
	VECTOR3 highGrid;
	int nBox[2];
	VECTOR3 boxLen[2];
	VECTOR3 resolutions[2];

	bool toTransformGeometry;

public:

	FDL_boxcountdim( VECTOR3 *lowgrid, VECTOR3 *highgrid, bool transformflag );

	double computeBCD_SingleStreamline_SingleWindow( VECTOR3* trace, int traceLength );
	double computeBCD_SingleStreamline_Full( FDL_streamline* streamline );
	double computeBCD_SingleStreamline_InAxisAlignedBox( VECTOR3* trace, int traceLength, double *boxMin, double *boxMax );
	void computeBCD_SingleStreamline_windowed( FDL_streamline* streamline, int windowLength, double** bcdArray );
	void computeBCD_SingleStreamline_DisjointSegments( FDL_streamline* streamline, int segmentLength, double* bcdArray );

	double computeBCD_MultiStreamline_SingleWindow( FDL_streamline* streamlineList,
														  int nStreamline, int *windowIdList,
														  int* windowLengthList );
	void computeBCD_MultiStreamline_Full( FDL_streamline* streamlineList, int nStreamline, double* bcdArray );
	void computeBCD_MultiStreamline_DisjointSegments( FDL_streamline* streamlineList,
															int nStreamline, int segmentLength,
															double** segmentedFdArray );

	double computeBC_SingleResolution( VECTOR3* trace, int traceLength, int iRes );
	void computeBC_SingleResolution2( VECTOR3* trace, int traceLength, int iRes, set<VECTOR3,Comparator2> *boxidset );

	int computeNextWindowSize_VoxelBased( FDL_streamline* streamline, int startIndex, int voxelLimit );
	int computeNextWindowSize_ArclengthBased( FDL_streamline* streamline, int startIndex, double arclengthLimit );

	int positionToBox( VECTOR3 position, int iRes );
	VECTOR3 positionToBox3D( VECTOR3 position, int iRes );
	int boxId3DTo1D( VECTOR3 id, int iRes );

	double computeLengthAlongRay( VECTOR3 startpoint, VECTOR3 endpoint );
	int countBoxIntersectAlongRay_MD( VECTOR3 startpoint, VECTOR3 endpoint, int iRes );
	void countBoxIntersectAlongRay_DDA_Set( VECTOR3 startpoint, VECTOR3 endpoint, int iRes, set<VECTOR3, Comparator2> *boxidset );
	int countBoxIntersectAlongRay_DDA_Array( VECTOR3 startPoint, VECTOR3 endPoint, int iRes, bool *countCheckArray );

	void setGridResolutionPair(  double boxLenLow, double boxLenHigh );

	void transformGeometry( VECTOR3* trace, int traceLength );
	void transformHotelling( VECTOR3* trace, VECTOR3* newTrace, int traceLength );
	void transformHotelling2( VECTOR3* cHullTrace, VECTOR3* trace, VECTOR3* newTrace, int hullLength, int traceLength );
	void transformPCA( VECTOR3* trace, VECTOR3* newTrace, int traceLength );

	void transformToLocal_FirstPointBased( VECTOR3* trace, VECTOR3* newTrace, int traceLength );
	void transformToLocal_CenterBased( VECTOR3* trace, VECTOR3* newTrace, int traceLength );

	void transformTranslate( VECTOR3* trace, VECTOR3* newTrace, int traceLength, VECTOR3* delta );
	VECTOR3 transformMeanShift( int traceLength, VECTOR3* trace, VECTOR3* newTrace );
	void transformRotate( VECTOR3* trace, VECTOR3* newTrace, int traceLength, MATRIX3* rot );
};

#endif
/* FDL_BOXCOUNTDIM_H_ */
