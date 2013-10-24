#pragma once

#include "OSUFlow.h"

/*
void
htApply
(
	const list<VECTOR3*>& lv3Src,
	list<VECTOR3*>& lv3Dst
);
*/

#include <cxcore.h>

template<class T, int D>
void
htApply
(
	const list<T*>& lpv3Src,
	list<T*>& lpv3Dst
)
{
	// DEL-BY-LEETEN 08/27/2011:		assert(T::Dimension() < D);
	CvMat* cvSrc = cvCreateMat( 
		lpv3Src.size(), 
		D,			// 3D vector
		CV_32FC1	// 32 bit floating point in a single channel
		);
		
	CvMat* cvDst = cvCreateMat( 
		lpv3Src.size(), 
		D,		// 3D vector
		CV_32FC1	// 32 bit floating point in a single channel
		);
	
	// convert the data to an array of the CvMat type
	int p = 0;
	for(typename list<T*>::const_iterator
			ilpv3Iter = lpv3Src.begin(); 
		ilpv3Iter!= lpv3Src.end(); 
		ilpv3Iter++, p++) 
	{
		const T* pv3Coord = *ilpv3Iter;
		const T& v3Coord = *pv3Coord;
		for(int c = 0; c < D; c++)
		{
			CvScalar s;
			s.val[0] = v3Coord(c);
			cvSet2D( cvSrc, p, c, s);
		}
	}

	// allocate space 
	CvMat* cvMean = cvCreateMat( 
		1, D, // 3D vector
		CV_32FC1	// 32 bit floating point in a single channel
		);
	CvMat* cvEigValues = cvCreateMat( 
		1, D, // 3D vector
		CV_32FC1	// 32 bit floating point in a single channel
		);
	CvMat* cvEigVectors  = cvCreateMat( 
		D, D, // 3D vector
		CV_32FC1	// 32 bit floating point in a single channel
		);
		
	// apply PCA
	cvCalcPCA(cvSrc, cvMean, cvEigValues, cvEigVectors, CV_PCA_DATA_AS_ROW);
	cvProjectPCA( cvSrc, cvMean, cvEigVectors, cvDst );

	// find the lowest corner of the bounding box
	double pdMin[D];
	double pdMax[D];
	for(int c = 0; c < D; c++)
	{
		pdMin[c] = HUGE_VAL;
		pdMax[c] = HUGE_VAL;
	}
	p = 0;
	for(typename list<T*>::iterator
			ilpv3Iter = lpv3Dst.begin(); 
		ilpv3Iter!= lpv3Dst.end(); 
		ilpv3Iter++, p++) 
	{
		T* pv3Coord = *ilpv3Iter;
		T& v3Coord = *pv3Coord;
		for(int c = 0; c < D; c++)
		{
			CvScalar s;
			s = cvGet2D( cvDst, p, c);
			pdMin[c] = min(pdMin[c], s.val[0]);
			pdMax[c] = max(pdMax[c], s.val[0]);
		}
	}

	// convert the transformed points to the list<VECTOR4*> format
	p = 0;
	for(typename list<T*>::iterator
			ilpv3Iter = lpv3Dst.begin(); 
		ilpv3Iter!= lpv3Dst.end(); 
		ilpv3Iter++, p++) 
	{
		T* pv3Coord = *ilpv3Iter;
		T& v3Coord = *pv3Coord;
		for(int c = 0; c < D; c++)
		{
			CvScalar s;
			s = cvGet2D(cvDst, p, c);
			v3Coord[c] = s.val[0] - pdMin[c];
		}
	}

	// release the resource
	cvFree(&cvSrc);
	cvFree(&cvMean);
	cvFree(&cvEigValues);
	cvFree(&cvEigVectors );
	cvFree(&cvDst);
}

template<class T, int D>
void
htApply2
(
	const list<T*>& lpv3Src,
	const list<T*>& lpv3AllSrc,
	list<T*>& lpv3Dst,
	list<T*>& lpv3AllDst
)
{
	// DEL-BY-LEETEN 08/27/2011:		assert(T::Dimension() < D);
	CvMat* cvSrc = cvCreateMat(
		lpv3Src.size(),
		D,			// 3D vector
		CV_32FC1	// 32 bit floating point in a single channel
		);

	CvMat* cvDst = cvCreateMat(
		lpv3Src.size(),
		D,		// 3D vector
		CV_32FC1	// 32 bit floating point in a single channel
		);

	// convert the data to an array of the CvMat type
	int p = 0;
	for(typename list<T*>::const_iterator
			ilpv3Iter = lpv3Src.begin();
		ilpv3Iter!= lpv3Src.end();
		ilpv3Iter++, p++)
	{
		const T* pv3Coord = *ilpv3Iter;
		const T& v3Coord = *pv3Coord;
		for(int c = 0; c < D; c++)
		{
			CvScalar s;
			s.val[0] = v3Coord(c);
			cvSet2D( cvSrc, p, c, s);
		}
	}

	//////////////////////////////////////////////////
	CvMat* cvAllSrc = cvCreateMat(
		lpv3AllSrc.size(),
		D,			// 3D vector
		CV_32FC1	// 32 bit floating point in a single channel
		);

	CvMat* cvAllDst = cvCreateMat(
		lpv3AllSrc.size(),
		D,		// 3D vector
		CV_32FC1	// 32 bit floating point in a single channel
		);

	p = 0;
	for(typename list<T*>::const_iterator
			ilpv3Iter = lpv3AllSrc.begin();
		ilpv3Iter!= lpv3AllSrc.end();
		ilpv3Iter++, p++)
	{
		const T* pv3Coord = *ilpv3Iter;
		const T& v3Coord = *pv3Coord;
		for(int c = 0; c < D; c++)
		{
			CvScalar s;
			s.val[0] = v3Coord(c);
			cvSet2D( cvAllSrc, p, c, s);
		}
	}
	//////////////////////////////////////////////////

	// allocate space
	CvMat* cvMean = cvCreateMat(
		1, D, // 3D vector
		CV_32FC1	// 32 bit floating point in a single channel
		);
	CvMat* cvEigValues = cvCreateMat(
		1, D, // 3D vector
		CV_32FC1	// 32 bit floating point in a single channel
		);
	CvMat* cvEigVectors  = cvCreateMat(
		D, D, // 3D vector
		CV_32FC1	// 32 bit floating point in a single channel
		);

	// apply PCA
	cvCalcPCA( cvSrc, cvMean, cvEigValues, cvEigVectors, CV_PCA_DATA_AS_ROW );

	// Project all points using PCA
	cvProjectPCA( cvAllSrc, cvMean, cvEigVectors, cvAllDst );

	// find the lowest corner of the bounding box
	double pdMin[D];
	double pdMax[D];
	for(int c = 0; c < D; c++)
	{
		pdMin[c] = HUGE_VAL;
		pdMax[c] = HUGE_VAL;
	}
	p = 0;
	for(typename list<T*>::iterator
			ilpv3Iter = lpv3AllDst.begin();
		ilpv3Iter!= lpv3AllDst.end();
		ilpv3Iter++, p++)
	{
		T* pv3Coord = *ilpv3Iter;
		T& v3Coord = *pv3Coord;
		for(int c = 0; c < D; c++)
		{
			CvScalar s;
			s = cvGet2D( cvAllDst, p, c);
			pdMin[c] = min(pdMin[c], s.val[0]);
			pdMax[c] = max(pdMax[c], s.val[0]);
		}
	}

	// convert the transformed points to the list<VECTOR4*> format
	p = 0;
	for(typename list<T*>::iterator
			ilpv3Iter = lpv3AllDst.begin();
		ilpv3Iter!= lpv3AllDst.end();
		ilpv3Iter++, p++)
	{
		T* pv3Coord = *ilpv3Iter;
		T& v3Coord = *pv3Coord;
		for(int c = 0; c < D; c++)
		{
			CvScalar s;
			s = cvGet2D(cvAllDst, p, c);
			v3Coord[c] = s.val[0] - pdMin[c];
		}
	}

	// release the resource
	cvFree(&cvSrc);
	cvFree(&cvMean);
	cvFree(&cvEigValues);
	cvFree(&cvEigVectors );
	cvFree(&cvDst);
}

/*

$Log: libht.h,v $
Revision 1.1.1.1  2011-08-28 05:38:59  leeten

[08/27/2011]
1. [1ST] Add the library to do Hotelling transform for streamlines/pathlines.


*/

