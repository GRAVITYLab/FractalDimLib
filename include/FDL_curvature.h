/*
 * FDL_curvature.h
 *
 *  Created on: Aug 1, 2011
 *  Author: abon
 */

#ifndef FDL_CURVATURE_H_
#define FDL_CURVATURE_H_

#include "FDL_util.h"
#include "FDL_streamline.h"

class FDL_curvature
{
public:

	// Constructor
	FDL_curvature();

	// Destructor
	~FDL_curvature();

	double computeCurvature_SingleStreamline_Point( VECTOR3* trace );
	double computeCurvature_SingleStreamline_SingleWindow( VECTOR3* trace, int traceLength );
	double computeCurvature_SingleStreamline_Full( FDL_streamline* streamline,
														double* localCurvatureArray = NULL,
														double* minScore = NULL, double* maxScore = NULL );
	void computeCurvature_SingleStreamline_DisjointSegments( FDL_streamline* streamline,
																	int segmentLength, int nSegment,
																	double* curvatureArray );
	double computeCurvature_SingleStreamline_InAxisAlignedBox( VECTOR3* trace,
																	   int traceLength,
																	   double* boxMin, double* boxMax );
	double* computeCurvature_SingleStreamline_windowed( FDL_streamline* streamline, int windowLength );

	void computeCurvature_MultiStreamline_Full( FDL_streamline* streamlineList,
													  int nStreamline, double* curvatureArray,
													  double* globalMinScore = NULL, double* globalMaxScore = NULL );
	void computeCurvature_MultiStreamline_DisjointSegments( FDL_streamline* streamlineList,
																   int nStreamline, int segmentLength,
																   int *nSegmentList, double** segmentedCurvatureArray );
	double computeCurvature_MultiStreamline_SingleWindow( FDL_streamline* streamlineList,
																 int nStreamline,
																 int *windowIdList, int* windowLengthList );

};

#endif
