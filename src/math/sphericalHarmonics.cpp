#include "base/config.h"
#include "sphericalHarmonics.h"

namespace NIBR {

    namespace SH {

        int 	sphericalHarmonicOrder             		= 0;
        bool    isEven                                  = true;
        int 	numberOfSphericalHarmonicCoefficients  	= 0;

        size_t  numberOfSamples_phi 					= 0;
        size_t  numberOfSamples_theta 					= 0;

        float   scalingFactor_phi 						= 0;
        float   scalingFactor_theta 				    = 0;

        float  *precomputedPhiComponent 			    = NULL;
        float  *precomputedThetaComponent            	= NULL;

        OrderOfDirections orderOfDirs                   = ORDEROFDIRECTIONS_NOTSET;

        std::function<void(float *)>* orderDirections   = NULL;

    }

}


using namespace NIBR;
using namespace NIBR::SH;

int NIBR::SH::getNumberOfSHCoeffs() {return NIBR::SH::numberOfSphericalHarmonicCoefficients;}

NIBR::OrderOfDirections NIBR::convertOrderOfDirections(std::string ood) {
        
    if      (ood=="")    return XYZ;
    else if (ood=="XYZ") return XYZ;
    else if (ood=="XYz") return XYz;
    else if (ood=="XyZ") return XyZ;
    else if (ood=="Xyz") return Xyz;
    else if (ood=="xYZ") return xYZ;
    else if (ood=="xYz") return xYz;
    else if (ood=="xyZ") return xyZ;
    else if (ood=="xyz") return xyz;
    else if (ood=="XZY") return XZY;
    else if (ood=="XZy") return XZy;
    else if (ood=="XzY") return XzY;
    else if (ood=="Xzy") return Xzy;
    else if (ood=="xZY") return xZY;
    else if (ood=="xZy") return xZy;
    else if (ood=="xzY") return xzY;
    else if (ood=="xzy") return xzy;
    else if (ood=="YXZ") return YXZ;
    else if (ood=="YXz") return YXz;
    else if (ood=="YxZ") return YxZ;
    else if (ood=="Yxz") return Yxz;
    else if (ood=="yXZ") return yXZ;
    else if (ood=="yXz") return yXz;
    else if (ood=="yxZ") return yxZ;
    else if (ood=="yxz") return yxz;
    else if (ood=="YZX") return YZX;
    else if (ood=="YZx") return YZx;
    else if (ood=="YzX") return YzX;
    else if (ood=="Yzx") return Yzx;
    else if (ood=="yZX") return yZX;
    else if (ood=="yZx") return yZx;
    else if (ood=="yzX") return yzX;
    else if (ood=="yzx") return yzx;
    else if (ood=="ZYX") return ZYX;
    else if (ood=="ZYx") return ZYx;
    else if (ood=="ZyX") return ZyX;
    else if (ood=="Zyx") return Zyx;
    else if (ood=="zYX") return zYX;
    else if (ood=="zYx") return zYx;
    else if (ood=="zyX") return zyX;
    else if (ood=="zyx") return zyx;
    else if (ood=="ZXY") return ZXY;
    else if (ood=="ZXy") return ZXy;
    else if (ood=="ZxY") return ZxY;
    else if (ood=="Zxy") return Zxy;
    else if (ood=="zXY") return zXY;
    else if (ood=="zXy") return zXy;
    else if (ood=="zxY") return zxY;
    else if (ood=="zxy") return zxy;
    else {
		NIBR::disp(MSG_FATAL,"Unknown order of directions. Acceptable formats are like \"Xyz\"");
        return XYZ;
	}
    
}

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

void orderDirections_XYZ(float*  ) {return;}
void orderDirections_XYz(float* d) {float r[3]={d[0],d[1],d[2]}; d[0] =  r[0]; d[1] =  r[1]; d[2] = -r[2];}
void orderDirections_XyZ(float* d) {float r[3]={d[0],d[1],d[2]}; d[0] =  r[0]; d[1] = -r[1]; d[2] =  r[2];}
void orderDirections_Xyz(float* d) {float r[3]={d[0],d[1],d[2]}; d[0] =  r[0]; d[1] = -r[1]; d[2] = -r[2];}
void orderDirections_xYZ(float* d) {float r[3]={d[0],d[1],d[2]}; d[0] = -r[0]; d[1] =  r[1]; d[2] =  r[2];}
void orderDirections_xYz(float* d) {float r[3]={d[0],d[1],d[2]}; d[0] = -r[0]; d[1] =  r[1]; d[2] = -r[2];}
void orderDirections_xyZ(float* d) {float r[3]={d[0],d[1],d[2]}; d[0] = -r[0]; d[1] = -r[1]; d[2] =  r[2];}
void orderDirections_xyz(float* d) {float r[3]={d[0],d[1],d[2]}; d[0] = -r[0]; d[1] = -r[1]; d[2] = -r[2];}

void orderDirections_XZY(float* d) {float r[3]={d[0],d[2],d[1]}; d[0] =  r[0]; d[1] =  r[1]; d[2] =  r[2];}
void orderDirections_XZy(float* d) {float r[3]={d[0],d[2],d[1]}; d[0] =  r[0]; d[1] =  r[1]; d[2] = -r[2];}
void orderDirections_XzY(float* d) {float r[3]={d[0],d[2],d[1]}; d[0] =  r[0]; d[1] = -r[1]; d[2] =  r[2];}
void orderDirections_Xzy(float* d) {float r[3]={d[0],d[2],d[1]}; d[0] =  r[0]; d[1] = -r[1]; d[2] = -r[2];}
void orderDirections_xZY(float* d) {float r[3]={d[0],d[2],d[1]}; d[0] = -r[0]; d[1] =  r[1]; d[2] =  r[2];}
void orderDirections_xZy(float* d) {float r[3]={d[0],d[2],d[1]}; d[0] = -r[0]; d[1] =  r[1]; d[2] = -r[2];}
void orderDirections_xzY(float* d) {float r[3]={d[0],d[2],d[1]}; d[0] = -r[0]; d[1] = -r[1]; d[2] =  r[2];}
void orderDirections_xzy(float* d) {float r[3]={d[0],d[2],d[1]}; d[0] = -r[0]; d[1] = -r[1]; d[2] = -r[2];}

void orderDirections_YXZ(float* d) {float r[3]={d[1],d[0],d[2]}; d[0] =  r[0]; d[1] =  r[1]; d[2] =  r[2];}
void orderDirections_YXz(float* d) {float r[3]={d[1],d[0],d[2]}; d[0] =  r[0]; d[1] =  r[1]; d[2] = -r[2];}
void orderDirections_YxZ(float* d) {float r[3]={d[1],d[0],d[2]}; d[0] =  r[0]; d[1] = -r[1]; d[2] =  r[2];}
void orderDirections_Yxz(float* d) {float r[3]={d[1],d[0],d[2]}; d[0] =  r[0]; d[1] = -r[1]; d[2] = -r[2];}
void orderDirections_yXZ(float* d) {float r[3]={d[1],d[0],d[2]}; d[0] = -r[0]; d[1] =  r[1]; d[2] =  r[2];}
void orderDirections_yXz(float* d) {float r[3]={d[1],d[0],d[2]}; d[0] = -r[0]; d[1] =  r[1]; d[2] = -r[2];}
void orderDirections_yxZ(float* d) {float r[3]={d[1],d[0],d[2]}; d[0] = -r[0]; d[1] = -r[1]; d[2] =  r[2];}
void orderDirections_yxz(float* d) {float r[3]={d[1],d[0],d[2]}; d[0] = -r[0]; d[1] = -r[1]; d[2] = -r[2];}

void orderDirections_YZX(float* d) {float r[3]={d[1],d[2],d[0]}; d[0] =  r[0]; d[1] =  r[1]; d[2] =  r[2];}
void orderDirections_YZx(float* d) {float r[3]={d[1],d[2],d[0]}; d[0] =  r[0]; d[1] =  r[1]; d[2] = -r[2];}
void orderDirections_YzX(float* d) {float r[3]={d[1],d[2],d[0]}; d[0] =  r[0]; d[1] = -r[1]; d[2] =  r[2];}
void orderDirections_Yzx(float* d) {float r[3]={d[1],d[2],d[0]}; d[0] =  r[0]; d[1] = -r[1]; d[2] = -r[2];}
void orderDirections_yZX(float* d) {float r[3]={d[1],d[2],d[0]}; d[0] = -r[0]; d[1] =  r[1]; d[2] =  r[2];}
void orderDirections_yZx(float* d) {float r[3]={d[1],d[2],d[0]}; d[0] = -r[0]; d[1] =  r[1]; d[2] = -r[2];}
void orderDirections_yzX(float* d) {float r[3]={d[1],d[2],d[0]}; d[0] = -r[0]; d[1] = -r[1]; d[2] =  r[2];}
void orderDirections_yzx(float* d) {float r[3]={d[1],d[2],d[0]}; d[0] = -r[0]; d[1] = -r[1]; d[2] = -r[2];}

void orderDirections_ZYX(float* d) {float r[3]={d[2],d[1],d[0]}; d[0] =  r[0]; d[1] =  r[1]; d[2] =  r[2];}
void orderDirections_ZYx(float* d) {float r[3]={d[2],d[1],d[0]}; d[0] =  r[0]; d[1] =  r[1]; d[2] = -r[2];}
void orderDirections_ZyX(float* d) {float r[3]={d[2],d[1],d[0]}; d[0] =  r[0]; d[1] = -r[1]; d[2] =  r[2];}
void orderDirections_Zyx(float* d) {float r[3]={d[2],d[1],d[0]}; d[0] =  r[0]; d[1] = -r[1]; d[2] = -r[2];}
void orderDirections_zYX(float* d) {float r[3]={d[2],d[1],d[0]}; d[0] = -r[0]; d[1] =  r[1]; d[2] =  r[2];}
void orderDirections_zYx(float* d) {float r[3]={d[2],d[1],d[0]}; d[0] = -r[0]; d[1] =  r[1]; d[2] = -r[2];}
void orderDirections_zyX(float* d) {float r[3]={d[2],d[1],d[0]}; d[0] = -r[0]; d[1] = -r[1]; d[2] =  r[2];}
void orderDirections_zyx(float* d) {float r[3]={d[2],d[1],d[0]}; d[0] = -r[0]; d[1] = -r[1]; d[2] = -r[2];}

void orderDirections_ZXY(float* d) {float r[3]={d[2],d[0],d[1]}; d[0] =  r[0]; d[1] =  r[1]; d[2] =  r[2];}
void orderDirections_ZXy(float* d) {float r[3]={d[2],d[0],d[1]}; d[0] =  r[0]; d[1] =  r[1]; d[2] = -r[2];}
void orderDirections_ZxY(float* d) {float r[3]={d[2],d[0],d[1]}; d[0] =  r[0]; d[1] = -r[1]; d[2] =  r[2];}
void orderDirections_Zxy(float* d) {float r[3]={d[2],d[0],d[1]}; d[0] =  r[0]; d[1] = -r[1]; d[2] = -r[2];}
void orderDirections_zXY(float* d) {float r[3]={d[2],d[0],d[1]}; d[0] = -r[0]; d[1] =  r[1]; d[2] =  r[2];}
void orderDirections_zXy(float* d) {float r[3]={d[2],d[0],d[1]}; d[0] = -r[0]; d[1] =  r[1]; d[2] = -r[2];}
void orderDirections_zxY(float* d) {float r[3]={d[2],d[0],d[1]}; d[0] = -r[0]; d[1] = -r[1]; d[2] =  r[2];}
void orderDirections_zxy(float* d) {float r[3]={d[2],d[0],d[1]}; d[0] = -r[0]; d[1] = -r[1]; d[2] = -r[2];}

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

    delete NIBR::SH::orderDirections;
}


void NIBR::SH::precompute(int _sphericalHarmonicOrder, OrderOfDirections _orderOfDirs, size_t _num) {
    
    if (NIBR::SH::precomputedPhiComponent != NULL) {return;}

	NIBR::SH::numberOfSamples_phi 	    =   _num;
	NIBR::SH::numberOfSamples_theta     = 8*_num;

	NIBR::SH::sphericalHarmonicOrder    = _sphericalHarmonicOrder;
    NIBR::SH::orderOfDirs               = _orderOfDirs;

    NIBR::SH::orderDirections = new std::function<void(float *)>();
    
    switch (NIBR::SH::orderOfDirs)
    {  
        case XYZ : *NIBR::SH::orderDirections = &orderDirections_XYZ; break;
        case XYz : *NIBR::SH::orderDirections = &orderDirections_XYz; break;
        case XyZ : *NIBR::SH::orderDirections = &orderDirections_XyZ; break;
        case Xyz : *NIBR::SH::orderDirections = &orderDirections_Xyz; break;
        case xYZ : *NIBR::SH::orderDirections = &orderDirections_xYZ; break;
        case xYz : *NIBR::SH::orderDirections = &orderDirections_xYz; break;
        case xyZ : *NIBR::SH::orderDirections = &orderDirections_xyZ; break;
        case xyz : *NIBR::SH::orderDirections = &orderDirections_xyz; break;
        
        case XZY : *NIBR::SH::orderDirections = &orderDirections_XZY; break; 
        case XZy : *NIBR::SH::orderDirections = &orderDirections_XZy; break;
        case XzY : *NIBR::SH::orderDirections = &orderDirections_XzY; break; 
        case Xzy : *NIBR::SH::orderDirections = &orderDirections_Xzy; break;
        case xZY : *NIBR::SH::orderDirections = &orderDirections_xZY; break; 
        case xZy : *NIBR::SH::orderDirections = &orderDirections_xZy; break;
        case xzY : *NIBR::SH::orderDirections = &orderDirections_xzY; break; 
        case xzy : *NIBR::SH::orderDirections = &orderDirections_xzy; break;
        
        case YXZ : *NIBR::SH::orderDirections = &orderDirections_YXZ; break;
        case YXz : *NIBR::SH::orderDirections = &orderDirections_YXz; break;
        case YxZ : *NIBR::SH::orderDirections = &orderDirections_YxZ; break;
        case Yxz : *NIBR::SH::orderDirections = &orderDirections_Yxz; break;
        case yXZ : *NIBR::SH::orderDirections = &orderDirections_yXZ; break;
        case yXz : *NIBR::SH::orderDirections = &orderDirections_yXz; break;
        case yxZ : *NIBR::SH::orderDirections = &orderDirections_yxZ; break;
        case yxz : *NIBR::SH::orderDirections = &orderDirections_yxz; break;

        case YZX : *NIBR::SH::orderDirections = &orderDirections_YZX; break;
        case YZx : *NIBR::SH::orderDirections = &orderDirections_YZx; break;
        case YzX : *NIBR::SH::orderDirections = &orderDirections_YzX; break;
        case Yzx : *NIBR::SH::orderDirections = &orderDirections_Yzx; break;
        case yZX : *NIBR::SH::orderDirections = &orderDirections_yZX; break;
        case yZx : *NIBR::SH::orderDirections = &orderDirections_yZx; break;
        case yzX : *NIBR::SH::orderDirections = &orderDirections_yzX; break;
        case yzx : *NIBR::SH::orderDirections = &orderDirections_yzx; break;
        
        case ZYX : *NIBR::SH::orderDirections = &orderDirections_ZYX; break;
        case ZYx : *NIBR::SH::orderDirections = &orderDirections_ZYx; break;
        case ZyX : *NIBR::SH::orderDirections = &orderDirections_ZyX; break;
        case Zyx : *NIBR::SH::orderDirections = &orderDirections_Zyx; break;
        case zYX : *NIBR::SH::orderDirections = &orderDirections_zYX; break;
        case zYx : *NIBR::SH::orderDirections = &orderDirections_zYx; break;
        case zyX : *NIBR::SH::orderDirections = &orderDirections_zyX; break;
        case zyx : *NIBR::SH::orderDirections = &orderDirections_zyx; break;
        
        case ZXY : *NIBR::SH::orderDirections = &orderDirections_ZXY; break;
        case ZXy : *NIBR::SH::orderDirections = &orderDirections_ZXy; break;
        case ZxY : *NIBR::SH::orderDirections = &orderDirections_ZxY; break;
        case Zxy : *NIBR::SH::orderDirections = &orderDirections_Zxy; break;
        case zXY : *NIBR::SH::orderDirections = &orderDirections_zXY; break;
        case zXy : *NIBR::SH::orderDirections = &orderDirections_zXy; break;
        case zxY : *NIBR::SH::orderDirections = &orderDirections_zxY; break;
        case zxy : *NIBR::SH::orderDirections = &orderDirections_zxy; break;
        
        default: { } break;
    }
    
    
    
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
     
    auto preComputePhi = [&](NIBR::MT::TASK task)->void {

        size_t c        = task.no*NIBR::SH::numberOfSamples_phi*NIBR::SH::numberOfSphericalHarmonicCoefficients;
		double x 		= (double)(task.no)*delta_phi-1;
        
		for (size_t j=0; j<NIBR::SH::numberOfSamples_phi; j++) {

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
    
	auto preComputeTheta = [&](NIBR::MT::TASK task)->void {

        size_t c        = task.no*NIBR::SH::numberOfSphericalHarmonicCoefficients;
        
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
       
    for (size_t i=0; i<inpCoords.size(); i++) {

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

    auto getPhiIndex = [&]()->size_t {
        return  numberOfSphericalHarmonicCoefficients*((size_t)((unit_dir[0]+1)*scalingFactor_phi)*numberOfSamples_phi + (size_t)((unit_dir[1]+1)*scalingFactor_phi));
    };

    auto getThetaIndex = [&]()->size_t {
        return  (int)((unit_dir[2]+1)*scalingFactor_theta)*numberOfSphericalHarmonicCoefficients;
    };
    
    (*NIBR::SH::orderDirections)(unit_dir);
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
