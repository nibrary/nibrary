#include "base/nibr.h"
#include "base/verbose.h"
#include "base/multithreader.h"
#include "sphericalHarmonics.h"

namespace NIBR {

    namespace SH {

        int 	sphericalHarmonicOrder             		= 0;
        bool    isEven                                  = true;
        int 	numberOfSphericalHarmonicCoefficients  	= 0;

        std::size_t  numberOfSamples_phi 			    = 0;
        std::size_t  numberOfSamples_theta 				= 0;

        float   scalingFactor_phi 						= 0;
        float   scalingFactor_theta 				    = 0;

        float  *precomputedPhiComponent 			    = NULL;
        float  *precomputedThetaComponent            	= NULL;

        OrderOfDirections orderOfDirs                   = ORDEROFDIRECTIONS_NOTSET;

        std::function<void(float *)> orderDirections    = NIBR::reorientFun(XYZ);

    }

}


using namespace NIBR;
using namespace NIBR::SH;
using namespace NIBR::MT;

int NIBR::SH::getNumberOfSHCoeffs() {return NIBR::SH::numberOfSphericalHarmonicCoefficients;}

int NIBR::getNumberOfSHCoeffs(int order) {
    
    if (order%2 == 0) {
        return (order+3)*order/2+1;
    } else {
        return (order+1)*(order+1);
    }

}

int NIBR::getSHOrderFromNumberOfCoeffs(int numberOfCoefficients) {

    int shOrder = sqrt(numberOfCoefficients);

    if ((shOrder*shOrder)!=numberOfCoefficients) {
        shOrder  = (sqrt(8*numberOfCoefficients+1)-3)/2;
    } else {
        shOrder  = shOrder-1;
    }

    disp(MSG_DETAIL,"Spherical harmonics order of input is %d", shOrder);

    return shOrder;

}

int sphPlmInd(int l,int m) {

	if (NIBR::SH::isEven)
		return (l*(l+1))/2+m;
	else
		return l*l+l+m;

}

// x is cos(theta)
void computeLegendrePolynomials(double *plm, double x, int order) {

	plm[0] = 1.0/std::sqrt(4.0*M_PI);

	for(double m=1; m <= order; m++) {
		plm[sphPlmInd(m,m)] = -std::sqrt((2*m+1)*(1-x*x)/(2*m))*plm[sphPlmInd(m-1,m-1)];
	}

	for(double m=0; m < order; m++) {
		plm[sphPlmInd(m+1,m)] = std::sqrt(2*m+3)*x*plm[sphPlmInd(m,m)];
	}

	for(double m=0; m <= order; m++)
		for(double l = m+2; l <= order; l++) {
			plm[sphPlmInd(l,m)] = std::sqrt(((2.0*l+1)*(2.0*l-1)) / ((l+m)*(l-m)))*x*plm[sphPlmInd(l-1,m)]-std::sqrt( (2.0*l+1)*(l-m-1.0)*(l+m-1.0) / ((2.0*l-3)*(l-m)*(l+m)))*plm[sphPlmInd(l-2,m)];
		}
            
}

void NIBR::SH::clean() {
	delete[] NIBR::SH::precomputedPhiComponent;
	delete[] NIBR::SH::precomputedThetaComponent;
    
    NIBR::SH::precomputedPhiComponent   = NULL;
    NIBR::SH::precomputedThetaComponent = NULL;
}


void NIBR::SH::precompute(int _sphericalHarmonicOrder, OrderOfDirections _orderOfDirs, std::size_t _num) {
    
    if (NIBR::SH::precomputedPhiComponent != NULL) {return;}

	NIBR::SH::numberOfSamples_phi 	    =   _num;
	NIBR::SH::numberOfSamples_theta     = 8*_num;

	NIBR::SH::sphericalHarmonicOrder    = _sphericalHarmonicOrder;
    NIBR::SH::orderOfDirs               = _orderOfDirs;

    NIBR::SH::orderDirections           = reorientFun(_orderOfDirs);

    
    if (NIBR::SH::sphericalHarmonicOrder%2 == 0) {
        NIBR::SH::isEven = true;
        NIBR::SH::numberOfSphericalHarmonicCoefficients 	= (NIBR::SH::sphericalHarmonicOrder+3)*NIBR::SH::sphericalHarmonicOrder/2+1;
    } else {
        NIBR::SH::isEven = false;
        NIBR::SH::numberOfSphericalHarmonicCoefficients 	= (NIBR::SH::sphericalHarmonicOrder+1)*(NIBR::SH::sphericalHarmonicOrder+1);
    }

	double delta_phi 	 			    = 2/(double)(NIBR::SH::numberOfSamples_phi   - 1);
	double delta_theta 			        = 2/(double)(NIBR::SH::numberOfSamples_theta - 1);
    
	NIBR::SH::scalingFactor_phi 	    = 1/delta_phi;
	NIBR::SH::scalingFactor_theta 	    = 1/delta_theta;

	NIBR::SH::precomputedPhiComponent   = new float[NIBR::SH::numberOfSphericalHarmonicCoefficients*NIBR::SH::numberOfSamples_phi*NIBR::SH::numberOfSamples_phi];
	NIBR::SH::precomputedThetaComponent = new float[NIBR::SH::numberOfSphericalHarmonicCoefficients*NIBR::SH::numberOfSamples_theta];    
     
    auto preComputePhi = [&](const NIBR::MT::TASK& task)->void {

        std::size_t c        = task.no*NIBR::SH::numberOfSamples_phi*NIBR::SH::numberOfSphericalHarmonicCoefficients;
		double x 		= (double)(task.no)*delta_phi-1;
        
		for (std::size_t j=0; j<NIBR::SH::numberOfSamples_phi; j++) {

			double y 		= (double)j*delta_phi-1;
			double phi 		= std::atan2(y,x);

			NIBR::SH::precomputedPhiComponent[c++] = 1;

			if (NIBR::SH::isEven) {
				for(double l = 2; l <= NIBR::SH::sphericalHarmonicOrder; l+=2) {
					for(double m = -l; m <= l; m++) {
						double ang = (fabs((double)m))*phi;
						if (m<0)  		NIBR::SH::precomputedPhiComponent[c++] = std::sin(ang);
						else if (m==0)  NIBR::SH::precomputedPhiComponent[c++] = 1;
						else 			NIBR::SH::precomputedPhiComponent[c++] = std::cos(ang);
					}
				}
			}
			else {
				for(double l = 1; l <= NIBR::SH::sphericalHarmonicOrder; l+=1) {
					for(double m = -l; m <= l; m++) {
						double ang = (fabs((double)m))*phi;
						if (m<0)  		NIBR::SH::precomputedPhiComponent[c++] = std::sin(ang);
						else if (m==0)  NIBR::SH::precomputedPhiComponent[c++] = 1;
						else 			NIBR::SH::precomputedPhiComponent[c++] = std::cos(ang);
					}
				}
			}

		}
	};
	NIBR::MT::MTRUN(NIBR::SH::numberOfSamples_phi,NIBR::MT::MAXNUMBEROFTHREADS(),preComputePhi);
    
	auto preComputeTheta = [&](const NIBR::MT::TASK& task)->void {

        std::size_t c        = task.no*NIBR::SH::numberOfSphericalHarmonicCoefficients;
        
        double *plm     = new double[NIBR::SH::numberOfSphericalHarmonicCoefficients];
		double theta    = (double)(task.no)*delta_theta-1;
		computeLegendrePolynomials(plm, theta, NIBR::SH::sphericalHarmonicOrder);

		NIBR::SH::precomputedThetaComponent[c++] = plm[sphPlmInd(0,0)];

		if (NIBR::SH::isEven) {
			for(float l = 2; l <= NIBR::SH::sphericalHarmonicOrder; l+=2) {
				for(float m = -l; m <= l; m++) {
					if (m<0)  		NIBR::SH::precomputedThetaComponent[c++] = M_SQRT2*plm[sphPlmInd(l,-m)];
					else if (m==0) 	NIBR::SH::precomputedThetaComponent[c++] =         plm[sphPlmInd(l, 0)];
					else  			NIBR::SH::precomputedThetaComponent[c++] = M_SQRT2*plm[sphPlmInd(l, m)];
				}
			}
		} else {
			for(float l = 1; l <= NIBR::SH::sphericalHarmonicOrder; l+=1) {
				for(float m = -l; m <= l; m++) {
					if (m<0)  		NIBR::SH::precomputedThetaComponent[c++] = M_SQRT2*plm[sphPlmInd(l,-m)];
					else if (m==0) 	NIBR::SH::precomputedThetaComponent[c++] =         plm[sphPlmInd(l, 0)];
					else  			NIBR::SH::precomputedThetaComponent[c++] = M_SQRT2*plm[sphPlmInd(l, m)];
				}
			}
		}
		
		delete[] plm;
        
	};
	NIBR::MT::MTRUN(NIBR::SH::numberOfSamples_theta,NIBR::MT::MAXNUMBEROFTHREADS(),preComputeTheta);

    disp(MSG_DEBUG,"Spherical harmonics precomputation completed");
    
}

// SH = Ylm  * SF
// SF = Ylm' * SH
// Note that in the above SF is the value of the spherical function at a SINGLE POINT, and not the value x the AREA for a point! This distinction is important!
// For example, consider the case, when we want to convert a spherical functions, i.e. all the values on the sphere, and their directions,
// in the discrete case, to compute the total energy using a spherical integral, we need to multiply the values by (4*pi)/(number of directions), and add them.
// This is needed to ensure that the energy in the SH domain is same as that one in the SF domain.
// This is the reason why sh2sf function in image_operators makes that scaling.

void NIBR::SH_basis(std::vector<std::vector<float>>& Ylm, std::vector<std::vector<float>>& inpCoords, int order) {
    
    int coeffCount = (order%2 == 0) ? (order+3)*order/2+1 : (order+1)*(order+1);
    double *plm    = new double[coeffCount];
       
    for (std::size_t i=0; i<inpCoords.size(); i++) {

        std::vector<float> basis;
        
        double phi   = std::atan2(inpCoords[i][1],inpCoords[i][0]);
        double theta = inpCoords[i][2]; //for plm we need acos theta
        
        computeLegendrePolynomials(plm, theta, order);
        
        basis.push_back(plm[sphPlmInd(0,0)]);
        
        if (NIBR::SH::isEven) {
			for(float l = 2; l <= order; l+=2) {
                for (int m=-l; m<=l; m++) {
                    double ang = (fabs((double)m))*phi;
                    if (m<0)  		basis.push_back(SQRT2*std::sin(ang)*plm[sphPlmInd(l,-m)]);
                    else if (m==0) 	basis.push_back(plm[sphPlmInd(l, 0)]);
                    else  			basis.push_back(SQRT2*std::cos(ang)*plm[sphPlmInd(l, m)]);
                }
			}
			
		} else {
			for(float l = 1; l <= order; l+=1) {
                for (int m=-l; m<=l; m++) {
                    double ang = (fabs((double)m))*phi;
                    if (m<0)  		basis.push_back(SQRT2*std::sin(ang)*plm[sphPlmInd(l,-m)]);
                    else if (m==0) 	basis.push_back(plm[sphPlmInd(l, 0)]);
                    else  			basis.push_back(SQRT2*std::cos(ang)*plm[sphPlmInd(l, m)]);
                }
			}
		}

        Ylm.push_back(basis);
        
    }
    
    delete[] plm;

}


// Note that the below SF is the value of the spherical function at a SINGLE POINT, and not the value x the AREA for a point. 
// (Read above for details, and see the two sh2sf functions in image_operators, where one of them uses the above basis function and the other uses the precomputed values as below.)
float NIBR::SH::toSF(float *sh, float *dir) {
    
    float* unit_dir = new float[3];
    unit_dir[0]     = dir[0];
    unit_dir[1]     = dir[1];
    unit_dir[2]     = dir[2];

    auto getPhiIndex = [&]()->std::size_t {
        return  numberOfSphericalHarmonicCoefficients*((std::size_t)((unit_dir[0]+1)*scalingFactor_phi)*numberOfSamples_phi + (std::size_t)((unit_dir[1]+1)*scalingFactor_phi));
    };

    auto getThetaIndex = [&]()->std::size_t {
        return  (int)((unit_dir[2]+1)*scalingFactor_theta)*numberOfSphericalHarmonicCoefficients;
    };
    
    NIBR::SH::orderDirections(unit_dir);
    verifyUnitRange(unit_dir);
    float *phiComp 		= NIBR::SH::precomputedPhiComponent   +   getPhiIndex();
    float *thetaComp  	= NIBR::SH::precomputedThetaComponent + getThetaIndex();
    delete[] unit_dir;
    
    float amp = 0;
    for (int i=0; i<NIBR::SH::numberOfSphericalHarmonicCoefficients; i++)
        amp += sh[i]*phiComp[i]*thetaComp[i];

    if (amp>0) 	return amp;
    else 		return 0;

}
