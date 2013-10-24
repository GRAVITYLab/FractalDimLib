/*
 * FDL_algorithm_curvature.h
 *
 *  Created on: Feb 2, 2011
 *  Author: abon
 */

#ifndef FDL_ALGORITHM_CURVATURE_H_
#define FDL_ALGORITHM_CURVATURE_H_

#include "FDL_header.h"
#include "FDL_streamline.h"
#include "FDL_boxcountdim.h"
#include "FDL_curvature.h"
#include "FDL_featuremap.h"

class FDL_algorithm_curvature
{
public:
	
	FDL_curvature *curvatureComputer;

	int nPartition[3];
	int nTotalBox;
	double *boxBoundary;
	double *scoreList;
	
	list<struct feature_in_box> partitionwiseFeatureList;
	double *sortedFeatureList;

public:

	// Constructor
	FDL_algorithm_curvature();

	void setCurvatureComputer( FDL_curvature *cvcomputer );

	// Feature map computation function 
	void compute_FeatureMap_SingleStreamline( FDL_streamline *streamline, int minLen );
	void compute_FeatureMap_SingleStreamline_Recursive( FDL_featuremap *featureMap,
															   VECTOR3 *trace,
															   int level, int startPoint,
															   int traceLength, int maxLevel );
	void compute_FeatureMap_MultiStreamline( FDL_streamline *streamlinelist, int nstreamline, int minlen );	

	// Segmentation based functions
	void segment_FeatureMapBased_Recursive_SingleStreamline(  FDL_streamline *streamline );
	void segment_FeatureMapBased_Iterative_SingleStreamline( FDL_streamline *streamline  );
	int getNextSegment( FDL_streamline *streamline, list<struct Segment> tempList, struct Segment segment, int totalLength );
	void segment_FeatureMapBased_MultiStreamline( FDL_streamline *streamlineList, int nStreamline, int minLen );
	void joinSegments_SingleStreamline( FDL_streamline *streamline );

	// File I/O
	void save_SegmentationResult( FDL_streamline *streamlineList, int nstreamline, const char *featureFileName  );
	
	
	void subdivide_FlowField( VECTOR3 *min, VECTOR3 *max, int npartx, int nparty, int npartz );

	void computeFeatureList_AllBox( FDL_streamline *streamlineList, int nstreamline );
	void computeFeatureList_SingleBox( FDL_streamline *streamlineList, int nstreamline, int boxid );

	void computeFeatureList_AllBox2( FDL_streamline *streamlineList, int nstreamline, double* scoreList );
	double computeFeatureList_SingleBox2( FDL_streamline *streamlineList, int nstreamline, int boxid );

	/*
	void computeFeatureVectors_AllBox_AllStreamline( FDL_streamline *streamlinelist, int nstreamline, float *featureVector );	
	void computeFeatureVectors_AllBox_SingleStreamline( FDL_streamline *streamline, float *featureVector  );
	void saveFeatureVectors( float *featureVector, int nstreamline, const char *featureFileName  );

	void projectFeatures_AllStreamlines( int nstreamline, int K, float *featureSpan, float *featureScore );
	list<struct feature_in_box> listFeatures_SingleStreamline( int streamlineID );
	list<struct feature_in_box> listTopKFeatures_SingleStreamline( list<struct feature_in_box> featureList_currentStreamline, int K );
	float compute_TopKFeaturesProjection_SingleStreamline(  list<struct feature_in_box>  featureList_topK_currentStreamline, int K );
	float compute_FeatureSpan_SingleStreamline(  list<struct feature_in_box> featureList_topK_currentStreamline, int K );
	
	void computeFeatureVectors_Multiresolution( float *featureVector, FDL_streamline *streamlineList, int nStreamline, int K, int nMaxLevel  );
	void saveFeatureVectors( float *featureVector, int nData, int  nDim, const char *featureFileName  );
	
	
	static bool comparator_feature( struct feature_in_box first, struct feature_in_box second );	
	static bool comparator_feature2( float first, float second );
	*/
		
};

#endif /* FDL_MULTIRESALGO_H_ */
