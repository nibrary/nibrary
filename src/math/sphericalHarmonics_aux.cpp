#include "base/verbose.h"
#include "sphericalHarmonics_aux.h"

// Returns number of spherical harmonics coefficients for given order and whether odd coefficients are ignored
int NIBR::getNumberOfSHCoeffs(int order, bool ignoreOddCoeffs) {

    if (order > 33) {
        disp(MSG_ERROR,"Order %d is too large. Max allowed order is 33.", order);
        return 0;
    }

    if (order%2 == 1 && ignoreOddCoeffs) order--;

    int numberOfCoeffs = 0;

    if (ignoreOddCoeffs) {
        numberOfCoeffs = (order+3)*order/2+1;
    } else {
        numberOfCoeffs = (order+1)*(order+1);
    }

    disp(MSG_DETAIL,"Number of spherical harmonics coefficients for order %d and ignoreOddCoeffs %s is %d", order, ignoreOddCoeffs ? "true" : "false", numberOfCoeffs);

    return numberOfCoeffs;

}

// Returns spherical harmonics order and whether odd coefficients are ignored
std::tuple<int, bool> NIBR::getSHOrderFromNumberOfCoeffs(int numberOfCoefficients) {

    int  shOrder         = sqrt(numberOfCoefficients);

    bool ignoreOddCoeffs = false;

    if ((shOrder*shOrder)!=numberOfCoefficients) {
        shOrder         = (sqrt(8*numberOfCoefficients+1)-3)/2;
        ignoreOddCoeffs = true;
    } else {
        shOrder         = shOrder-1;
        ignoreOddCoeffs = false;
    }

    if (shOrder > 33) {
        disp(MSG_ERROR,"Order %d is too large. Max allowed order is 33.", shOrder);
        return std::make_tuple(0, false);
    }

    disp(MSG_DETAIL,"Spherical harmonics order of input is %d, ignoreOddCoeffs: %s", shOrder, ignoreOddCoeffs ? "true" : "false");

    return std::make_tuple(shOrder, ignoreOddCoeffs);

}


// x is cos(theta)
void NIBR::computeLegendrePolynomials(double *plm, double x, int order, bool ignoreOddCoeffs) {

	plm[0] = 1.0/std::sqrt(4.0*M_PI);

    auto sphPlmInd = [&] (int l,int m) -> int { 
        if (ignoreOddCoeffs) {
            return (l*(l+1))/2+m;
        } else {
            return l*l+l+m;
        }
    };

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

// SH = Ylm  * SF
// SF = Ylm' * SH
// Note that in the above SF is the value of the spherical function at a SINGLE POINT, and not the value x the AREA for a point! This distinction is important!
// For example, consider the case, when we want to convert a spherical functions, i.e. all the values on the sphere, and their directions,
// in the discrete case, to compute the total energy using a spherical integral, we need to multiply the values by (4*pi)/(number of directions), and add them.
// This is needed to ensure that the energy in the SH domain is same as that one in the SF domain.
// This is the reason why sh2sf function in image_operators makes that scaling.

void NIBR::SH_basis(std::vector<std::vector<float>>& Ylm, std::vector<Point3D>& inpCoords, int order, bool ignoreOddCoeffs) {
    
    int coeffCount = getNumberOfSHCoeffs(order, ignoreOddCoeffs);

    double *plm    = new double[coeffCount];
       
    for (std::size_t i=0; i<inpCoords.size(); i++) {

        std::vector<float> basis;
        
        double phi   = std::atan2(inpCoords[i][1],inpCoords[i][0]);
        double theta = inpCoords[i][2]; //for plm we need acos theta
        
        computeLegendrePolynomials(plm, theta, order, ignoreOddCoeffs);
        
        basis.push_back(plm[0]);
        
        if (ignoreOddCoeffs) {
			for(int l = 2; l <= order; l+=2) {
                for (int m=-l; m<=l; m++) {
                    double ang = (fabs((double)m))*phi;
                    if (m<0)  		basis.push_back(SQRT2*std::sin(ang)*plm[(l*(l+1))/2-m]);
                    else if (m==0) 	basis.push_back(plm[(l*(l+1))/2]);
                    else  			basis.push_back(SQRT2*std::cos(ang)*plm[(l*(l+1))/2+m]);
                }
			}
			
		} else {
			for(int l = 1; l <= order; l+=1) {
                for (int m=-l; m<=l; m++) {
                    double ang = (fabs((double)m))*phi;
                    if (m<0)  		basis.push_back(SQRT2*std::sin(ang)*plm[l*l+l-m]);
                    else if (m==0) 	basis.push_back(plm[l*l+l]);
                    else  			basis.push_back(SQRT2*std::cos(ang)*plm[l*l+l+m]);
                }
			}
		}

        Ylm.push_back(basis);
        
    }
    
    delete[] plm;

}