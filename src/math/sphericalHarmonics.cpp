#include "base/verbose.h"
#include "sphericalHarmonics.h"
#include "sphericalHarmonics_aux.h"

using namespace NIBR;

SH::SH(int _order, bool _noOddCoeffs, OrderOfDirections _orderOfDirs, size_t _numberOfSamples) {
    order                 = _order;
    noOddCoeffs           = _noOddCoeffs;
    orderOfDirs           = _orderOfDirs;
    numberOfSamples       = _numberOfSamples;
    numberOfSamples_phi   = numberOfSamples;
    numberOfSamples_theta = 8*numberOfSamples;
    orderDirections       = reorientFun(orderOfDirs);
    precompute();
}

void SH::clean() {
    if (precomputedPhiComponent   != NULL) {delete[] precomputedPhiComponent;}
    if (precomputedThetaComponent != NULL) {delete[] precomputedThetaComponent;}
    
    precomputedPhiComponent   = NULL;
    precomputedThetaComponent = NULL;
}


void SH::precompute() {

	if (order > 33) {
		disp(MSG_FATAL,"Order %d is too large. Max allowed order is 33.", order);
		order = 33;
		return;
	}
    
    clean();

    coeffCount          = getNumberOfSHCoeffs(order, noOddCoeffs);

	double delta_phi    = 2.0/(double)(numberOfSamples_phi   - 1);
	double delta_theta  = 2.0/(double)(numberOfSamples_theta - 1);
    
	scalingFactor_phi   = 1.0/delta_phi;
	scalingFactor_theta = 1.0/delta_theta;

	precomputedPhiComponent   = new float[coeffCount*numberOfSamples_phi*numberOfSamples_phi];
	precomputedThetaComponent = new float[coeffCount*numberOfSamples_theta];    
    
	// Precompute phi components
	for (std::size_t i=0; i<numberOfSamples_phi; i++) {

		std::size_t c 	= i*numberOfSamples_phi*coeffCount;
		double x 		= (double)(i)*delta_phi-1;
        
		for (std::size_t j=0; j<numberOfSamples_phi; j++) {

			double y 		= (double)j*delta_phi-1;
			double phi 		= std::atan2(y,x);

			precomputedPhiComponent[c++] = 1;

			if (noOddCoeffs) {    
				for(int l = 2; l <= order; l+=2) {
					for(int m = -l; m <= l; m++) {
						double ang = (std::fabs((double)m))*phi;
						if (m<0)  		precomputedPhiComponent[c++] = std::sin(ang);
						else if (m==0)  precomputedPhiComponent[c++] = 1;
						else 			precomputedPhiComponent[c++] = std::cos(ang);
					}
				}
			}
			else {
				for(int l = 1; l <= order; l+=1) {
					for(int m = -l; m <= l; m++) {
						double ang = (std::fabs((double)m))*phi;
						if (m<0)  		precomputedPhiComponent[c++] = std::sin(ang);
						else if (m==0)  precomputedPhiComponent[c++] = 1;
						else 			precomputedPhiComponent[c++] = std::cos(ang);
					}
				}
			}

		}
	}

    // Precompute theta components
	for (std::size_t i=0; i<numberOfSamples_theta; i++) {
		
		std::size_t c   = i*coeffCount;
        
        double *plm     = new double[coeffCount];
		double theta    = (double)(i)*delta_theta-1;

		computeLegendrePolynomials(plm, theta, order, noOddCoeffs);

		precomputedThetaComponent[c++] = plm[0];

		if (noOddCoeffs) {
			for(int l = 2; l <= order; l+=2) {
				for(int m = -l; m <= l; m++) {
					if (m<0)  		precomputedThetaComponent[c++] = M_SQRT2*plm[(l*(l+1))/2-m];
					else if (m==0) 	precomputedThetaComponent[c++] =         plm[(l*(l+1))/2];
					else  			precomputedThetaComponent[c++] = M_SQRT2*plm[(l*(l+1))/2+m];
				}
			}
		} else {
			for(int l = 1; l <= order; l+=1) {
				for(int m = -l; m <= l; m++) {
					if (m<0)  		precomputedThetaComponent[c++] = M_SQRT2*plm[l*l+l-m];
					else if (m==0) 	precomputedThetaComponent[c++] =         plm[l*l+l];
					else  			precomputedThetaComponent[c++] = M_SQRT2*plm[l*l+l+m];
				}
			}
		}
		
		delete[] plm;

	}

    disp(MSG_DEBUG,"Spherical harmonics precomputation completed");
    
}




// Note that the below SF is the value of the spherical function at a SINGLE POINT, and not the value x the AREA for a point. 
// (Read above for details, and see the two sh2sf functions in image_operators, where one of them uses the above basis function and the other uses the precomputed values as below.)
float SH::toSF(float *sh, float *dir) {
    
    float* unit_dir = new float[3];
    unit_dir[0]     = dir[0];
    unit_dir[1]     = dir[1];
    unit_dir[2]     = dir[2];

    auto getPhiIndex = [&]()->std::size_t {
        return  coeffCount*((std::size_t)((unit_dir[0]+1)*scalingFactor_phi)*numberOfSamples_phi + (std::size_t)((unit_dir[1]+1)*scalingFactor_phi));
    };

    auto getThetaIndex = [&]()->std::size_t {
        return  (int)((unit_dir[2]+1)*scalingFactor_theta)*coeffCount;
    };
    
    orderDirections(unit_dir);
    verifyUnitRange(unit_dir);
	
    float *phiComp 		= precomputedPhiComponent   +   getPhiIndex();
    float *thetaComp  	= precomputedThetaComponent + getThetaIndex();
    delete[] unit_dir;
    
    float amp = 0;
    for (int i=0; i<coeffCount; i++)
        amp += sh[i]*phiComp[i]*thetaComp[i];

    if (amp>0) 	return amp;
    else 		return 0;

}
