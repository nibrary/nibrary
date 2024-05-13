#include "recon_transhi2015.h"

using namespace NIBR;

void TranShi2015::print() {

    disp(MSG_INFO, "Running dMRI model: Tran & Shi, IEEE TMI, 2015");
    disp(MSG_INFO, "shOrder:                    %d",    shOrder);
    disp(MSG_INFO, "deltaStep:                  %.3e",  deltaStep);
    disp(MSG_INFO, "D_inAx:                     %.5f",  init_D_inAx);
    disp(MSG_INFO, "D_trapped:                  %.5f",  init_D_trapped);
    disp(MSG_INFO, "init_D_exAx_iso:            %.5f",  init_D_exAx_iso);
    disp(MSG_INFO, "init_minNumConstraint:      %d",    init_minNumConstraint);
    disp(MSG_INFO, "constraintUpdateCount:      %d",    constraintUpdateCount);
    if (fastOptimization) {
        disp(MSG_INFO, "fastOptimization:           ON");
    } else {
        disp(MSG_INFO, "fastOptimization:           OFF");
    }
    disp(MSG_INFO, "bValLow:                    %d",  int(bValLow));
    if (bValHigh>100000000000) {
        disp(MSG_INFO, "bValHigh:                   inf");
    } else {
        disp(MSG_INFO, "bValHigh:                   %d",  int(bValHigh));
    }
    disp(MSG_INFO, "maxIter:                    %d",    maxIter);
    disp(MSG_INFO, "xi_init:                    %.2f",  xi_init);
    disp(MSG_INFO, "xi_step:                    %.2f",  xi_step);
    disp(MSG_INFO, "xi_stepCount:               %d",    xi_stepCount);
    disp(MSG_INFO, "maxCrossings:               %d",    maxCrossings);
    disp(MSG_INFO, "noiseFloor:                 %.5f",  noiseFloor);        

}
