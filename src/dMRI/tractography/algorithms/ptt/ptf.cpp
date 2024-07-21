#include "ptf.h"

using namespace NIBR;

PTF::PTF() {
    
	p 			= new float[3];
    F           = new float*[3];
    F[0]        = new float[3];
    F[1]        = new float[3];
    F[2]        = new float[3];
    PP          = new float[9];
    
	likelihood 	= 0.0;
    firstVal    = NAN;
    
    angularSeparation = 2.0*M_PI/float(TRACKER::params_ptt.probeCount);
    probeStepSize     = TRACKER::params_ptt.probeLength/(TRACKER::params_ptt.probeQuality-1);
    if (TRACKER::params_ptt.img_FOD->getSHorder()%2==0) {
        probeNormalizer   = 1/float(TRACKER::params_ptt.probeQuality*TRACKER::params_ptt.probeCount);
    } else {
        probeNormalizer   = 1/float((TRACKER::params_ptt.probeQuality-1)*TRACKER::params_ptt.probeCount);
    }
    
}

PTF::~PTF() {
    delete[] p;
    delete[] F[0];
    delete[] F[1];
    delete[] F[2];
    delete[] F;
    delete[] PP;
}

void PTF::setPosition(float* _p) {
    p[0] = _p[0];
    p[1] = _p[1];
    p[2] = _p[2];
    disp(MSG_DEBUG, "PTF seed coordinate: [%.2f,%.2f,%.2f]",p[0],p[1],p[2]);
}

void PTF::copy(PTF *ptf) {
    k1       =  ptf->k1;
    k2       =  ptf->k2;
    k1_cand  =  ptf->k1_cand;
    k2_cand  =  ptf->k2_cand;
    
	for (int i=0; i<3; i++) {
        p[i] = ptf->p[i];
        for (int j=0; j<3; j++) {
            F[i][j] = ptf->F[i][j];
        }
        
	}
	
	for (int i=0; i<9; i++)
        PP[i]  = ptf->PP[i];
	
	likelihood 	      = ptf->likelihood;
    firstVal          = ptf->firstVal;
    initFirstVal      = ptf->initFirstVal;
    initPosteriorMax  = ptf->initPosteriorMax;
}

void PTF::getFlippedCopy(PTF *ptf) {

    this->copy(ptf);
    
    for (int i=0; i<3; i++) {
		F[0][i]  	*= -1;
		F[1][i] 	*= -1;
	}
	k1      *= -1;
    k1_cand *= -1;

	likelihood 	= 0.0;
    firstVal    = initFirstVal; // When we flip, we now recover the firstVal, which will be used in calcDataSupport for the next candidate
}


void PTF::print() {
	std::cout << "p:  " << p[0] << " " << p[1] << " " << p[2] << std::endl;
	std::cout << "k:  " << std::sqrt(k1*k1+k2*k2) << std::endl;
	std::cout << "k1: " << k1 << std::endl;
	std::cout << "k2: " << k2 << std::endl;
	std::cout << "T: "  << F[0][0] << " " <<  F[0][1] << " " <<  F[0][2] << std::endl;
	std::cout << "K1: " << F[1][0] << " " <<  F[1][1] << " " <<  F[1][2] << std::endl;
	std::cout << "K2: " << F[2][0] << " " <<  F[2][1] << " " <<  F[2][2] << std::endl;
	std::cout << "likelihood: "       << likelihood << std::endl;
    std::cout << "initPosteriorMax: " << initPosteriorMax << std::endl;
}


float PTF::calcDataSupport() {
    
    float _p[3];
    float _F[3][3];
    float _T[3] ={0,0,0};
    float _N1[3]={0,0,0};
    float _N2[3]={0,0,0};
    
    prepPropagator(probeStepSize);

    if (isnan(firstVal)) {
        firstVal = TRACKER::params_ptt.img_FOD->getFODamp(p,F[0]);
    }

    // Copy initial _p and _F    
    for (int i=0; i<3; i++) {
        _p[i] = p[i];
        for (int j=0; j<3; j++) {
            _F[i][j] = F[i][j];
        }   
	}
    
    if (TRACKER::params_ptt.img_FOD->getSHorder()%2==0) {
    
        likelihood = firstVal;
        
        for (int q=0; q<(TRACKER::params_ptt.probeQuality-1); q++) {
            
            
            for (int i=0; i<3; i++) {
                _p[i]  += PP[0]*_F[0][i] +  PP[1]*_F[1][i]  +  PP[2]*_F[2][i];
                _T[i]   = PP[3]*_F[0][i] +  PP[4]*_F[1][i]  +  PP[5]*_F[2][i];
            }
            normalize(_T);
            
            if (q<(TRACKER::params_ptt.probeQuality-1)) {
                
                for (int i=0; i<3; i++) {
                    _N2[i]  = PP[6]*_F[0][i] +  PP[7]*_F[1][i]  +  PP[8]*_F[2][i];
                }
                
                cross(_N1,_N2,_T);
                for (int i=0; i<3; i++) {
                    _F[0][i] =  _T[i];
                    _F[1][i] = _N1[i];
                    _F[2][i] = _N2[i];
                }
                
            }
            
            
            if (TRACKER::params_ptt.probeCount==1) {
                
                float val = TRACKER::params_ptt.img_FOD->getFODamp(_p,_T);
                
                if ((TRACKER::params_ptt.checkWeakLinks==true) && (val < TRACKER::params_ptt.weakLinkThresh)) {
                    likelihood  = 0;
                    return 0;
                } else {
                    likelihood += val;
                }
                
                
            } else {
                
                float totVal = 0;
                
                if (q==(TRACKER::params_ptt.probeQuality-1)) {
                    for (int i=0; i<3; i++) {
                        _N2[i]  = PP[6]*_F[0][i] +  PP[7]*_F[1][i]  +  PP[8]*_F[2][i];
                    }
                    cross(_N1,_N2,_T);
                }
                
                
                for (float c=0; c<TRACKER::params_ptt.probeCount; c++) {
                    
                    float pp[3];
                    
                    for (int i=0; i<3; i++) {
                        pp[i] = _p[i] + _N1[i]*TRACKER::params_ptt.probeRadius*std::cos(c*angularSeparation) + _N2[i]*TRACKER::params_ptt.probeRadius*std::sin(c*angularSeparation);
                    }
                    
                    float val = TRACKER::params_ptt.img_FOD->getFODamp(pp,_T);
                    
                    if ((TRACKER::params_ptt.checkWeakLinks==true) && (val < TRACKER::params_ptt.weakLinkThresh)) {
                        likelihood    = 0;
                        return 0;
                    } else {
                        totVal += val;
                    }

                } 
                
                likelihood += totVal;
                
            }
            
            prepPropagator(probeStepSize);
        }
        
    } else {
        
        likelihood = 0;
        
        float pn[3];
        float Tb[3];
        float Te[3];
        
        for (int q=0; q<(TRACKER::params_ptt.probeQuality-1); q++) {
            
            
            for (int i=0; i<3; i++) {
                pn[i]  = _p[i] + PP[0]*_F[0][i] +  PP[1]*_F[1][i]  +  PP[2]*_F[2][i];
                _T[i]  =         PP[3]*_F[0][i] +  PP[4]*_F[1][i]  +  PP[5]*_F[2][i];
               _N2[i]  =         PP[6]*_F[0][i] +  PP[7]*_F[1][i]  +  PP[8]*_F[2][i];
                Tb[i]  = pn[i] - _p[i];
            }
            normalize(_T);
            cross(_N1,_N2,_T);
            normalize(Tb);
            for (int i=0; i<3; i++) {
                Te[i]  = -Tb[i];
            }
            
            
            if (TRACKER::params_ptt.probeCount==1) {
                
                float link = (TRACKER::params_ptt.img_FOD->getFODamp(_p,Tb) + TRACKER::params_ptt.img_FOD->getFODamp(pn,Te))/float(2.0);
                
                if ((TRACKER::params_ptt.checkWeakLinks==true) && (link < TRACKER::params_ptt.weakLinkThresh)) {
                    likelihood  = 0;
                    return 0;
                } else {
                    likelihood += link;
                }
                
            } else {
                
                for (float c=0; c<TRACKER::params_ptt.probeCount; c++) {
                    
                    float pp[3];
                    float ppn[3];
                    
                    for (int i=0; i<3; i++) {
                        pp[i]  =  _p[i] + _F[1][i]*TRACKER::params_ptt.probeRadius*std::cos(c*angularSeparation) + _F[2][i]*TRACKER::params_ptt.probeRadius*std::sin(c*angularSeparation);
                        ppn[i] =  pn[i] +   _N1[i]*TRACKER::params_ptt.probeRadius*std::cos(c*angularSeparation) + _N2[i]*TRACKER::params_ptt.probeRadius*std::sin(c*angularSeparation);
                        Tb[i]  = ppn[i] -    pp[i];
                    }
                    
                    normalize(Tb);
                    for (int i=0; i<3; i++) {
                        Te[i]  = -Tb[i];
                    }
                    
                    float link = (TRACKER::params_ptt.img_FOD->getFODamp(pp,Tb) + TRACKER::params_ptt.img_FOD->getFODamp(ppn,Te))/float(2.0);
                    
                    if ((TRACKER::params_ptt.checkWeakLinks==true) && (link < TRACKER::params_ptt.weakLinkThresh)) {
                        likelihood  = 0;
                        return 0;
                    } else {
                        likelihood += link;
                    }

                }
                
            }
            
            // Update _F here
            if (q<(TRACKER::params_ptt.probeQuality-1)) {
                for (int i=0; i<3; i++) {
                       _p[i] = pn[i];
                    _F[0][i] =  _T[i];
                    _F[1][i] = _N1[i];
                    _F[2][i] = _N2[i];
                }
            }
            
            prepPropagator(probeStepSize);
        }
        
        
        
    }

    likelihood *= probeNormalizer;

    // if (TRACKER::params_ptt.dataSupportExponent != 1)
    //     likelihood  = std::pow(likelihood,TRACKER::params_ptt.dataSupportExponent);

    return likelihood;
}

float PTF::getCandidate() {
    doRandomThings.getARandomPointWithinDisk(&k1_cand, &k2_cand, TRACKER::params_ptt.maxCurvature);
    return calcDataSupport();
}

float PTF::getInitCandidate(float *initDir) {
    
    doRandomThings.getARandomMovingFrame(F,initDir);
    doRandomThings.getARandomPointWithinDisk(&k1_cand, &k2_cand, TRACKER::params_ptt.maxCurvature);
    k1 = k1_cand;
    k2 = k2_cand;   
    
    if (TRACKER::params_ptt.img_FOD->getSHorder()%2==0) {
        
        // First part of the probe
        likelihood = 0.0;
        
        if (TRACKER::params_ptt.probeCount==1) {
            likelihood = TRACKER::params_ptt.img_FOD->getFODamp(p,F[0]);
        } else {
            
            for (float c=0; c<TRACKER::params_ptt.probeCount; c++) {
                        
                float pp[3];
                
                for (int i=0; i<3; i++) {
                    pp[i] = p[i] + F[1][i]*TRACKER::params_ptt.probeRadius*std::cos(c*angularSeparation) + F[2][i]*TRACKER::params_ptt.probeRadius*std::sin(c*angularSeparation);
                }
                
                float val = TRACKER::params_ptt.img_FOD->getFODamp(pp,F[0]);
                
                if ((TRACKER::params_ptt.checkWeakLinks==true) && (val < TRACKER::params_ptt.weakLinkThresh)) {
                    likelihood  = 0;
                    return 0;
                } else {
                    likelihood += val;
                }

            }
            
        }
        
        // At init we use only lastVal because this will be added in calcDataSupport
        initFirstVal = likelihood;
        firstVal     = initFirstVal;
    }
    
    // lastVal is not used in the asymmetric FOD case
    return calcDataSupport();

}
