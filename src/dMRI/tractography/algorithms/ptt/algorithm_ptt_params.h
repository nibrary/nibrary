#pragma once

#include "dMRI/imageTypes/fod_image.h"

#define DEFAULT_PTT_STEPSIZE_IN_PIXELDIM                       0.025
#define DEFAULT_PTT_OUTPUTSTEPSIZE_IN_PIXELDIM                   0.5
#define DEFAULT_PTT_MINRADIUSOFCURVATURE_IN_PIXELDIM             0.5
#define DEFAULT_PTT_MINDATASUPPORT                              0.05
#define DEFAULT_PTT_DATASUPPORTEXPONENT                            1
#define DEFAULT_PTT_WEAKLINKTHRESH                                 0
#define DEFAULT_PTT_WRITEINTERVAL_IN_PIXELDIM                    0.5

#define DEFAULT_PTT_MAXESTINTERVAL                                 1
#define DEFAULT_PTT_INITMAXESTTRIALS                             100
#define DEFAULT_PTT_PROPMAXESTTRIALS                              20
#define DEFAULT_PTT_TRIESPERREJECTIONSAMPLING                   1000
#define DEFAULT_PTT_USEBESTATINIT                              false
#define DEFAULT_PTT_USELEGACYSAMPLING                          false
#define DEFAULT_PTT_SAMPLINGQUALITY                                2

#define DEFAULT_PTT_PROBELENGTH_IN_PIXELDIM                     0.25
#define DEFAULT_PTT_PROBERADIUS_IN_PIXELDIM                      0.5
#define DEFAULT_PTT_PROBECOUNT                                     1
#define DEFAULT_PTT_PROBECOUNT_WHEN_THEREIS_PROBERADIUS            4
#define DEFAULT_PTT_PROBEQUALITY                                   4

#define DEFAULT_PTT_MAXPOSTESTCOMPENS 	  		   			       2

namespace NIBR {

// Tracking parameters
class Params_PTT {

public:

	Params_PTT(){};
	~Params_PTT(){clear();};

	void needsUpdate() {pttIsReady = false;}
	void clear();
	bool update();
	void print();

	// FOD options in Trekker
	std::string           	 	img_FOD_path{""};
	std::string           	 	fod_sphere_path{""};
	bool 						fodIsSym{true};
	bool 						fodDiscretization{true};
	std::string           	    orderOfDirectionsTextInput{""};	
	
	// Tracking options
	float                		stepSize{NAN};
	float                		minRadiusOfCurvature{NAN};
	float                		minDataSupport{NAN};
	float                		dataSupportExponent{NAN};
	float                		weakLinkThresh{NAN};

	// Sampling options
	int             	 		maxEstInterval{-1};
	int			    	 		initMaxEstTrials{-1};
	int			    	 		propMaxEstTrials{-1};
	int 		         		triesPerRejectionSampling{-1};
	bool 		    	 		useBestAtInit{false}; 
	bool 		    	 		useLegacySampling{false}; 
	int			    	 		samplingQuality{-1};

	// Probe options
	float               		probeCount{NAN};
    float                		probeRadius{NAN};
    float               		probeLength{NAN};
	float                		probeQuality{NAN};

	// Output options
	bool 					 	saveFrame{false};
	float                       outputStep{NAN};

	// Derived parameters
	FOD_Image*           	 	img_FOD{NULL};
	float                		maxCurvature{NAN};
	bool                 		checkWeakLinks{false};
	float 				 		modMinDataSupport{NAN};

	int 								cdfVertCnt;
	int 								cdfFaceCnt;
	std::vector<std::array<int,3>>		cdfFace;
	std::vector<std::pair<float,float>> cdfk1k2;

private:

	void 						setDefaults();

	bool                        pttIsReady{false};

	std::vector<Point>			fodSphere;
	bool                        fodIsSpheresliced{false};
	float		    		    smallestPixDim{NAN};
	OrderOfDirections			orderOfDirections{XYZ};
	
};

}
