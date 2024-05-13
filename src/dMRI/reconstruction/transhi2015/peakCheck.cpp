#include "recon_transhi2015.h"

using namespace NIBR;

std::vector<float> TranShi2015::getPeakValues(const Eigen::VectorXd& y, double peakThresh) {

    std::vector<float> peakValues;

    for(int i = 0; i < y.size(); i++) {

        if (y(i) > peakThresh) {
            
            bool isPeak = true;

            for(const int& j : this->sphere.neighboringVertices[i]) {
                if(y(i) <= y(j)) {
                    isPeak = false;
                    break;
                }
            }

            if (isPeak) {
                peakValues.push_back(y(i));
            }

        }

    }

    if (peakValues.empty()) 
        return peakValues;
    
    // This part sorts, reverses and removes every other element in the array peakValues. 
    // This is done in order to take into account the 2 directions of each peak.
    std::sort(peakValues.begin(), peakValues.end(), std::greater<>());
    
    std::vector<float> everyOtherElement;

    for (size_t i = 0; i < peakValues.size(); i += 2) {
        everyOtherElement.push_back(peakValues[i]);
    }

    return everyOtherElement;
}

// Determine the target number of peaks:
//       max(a)> 0.5:   FOD_THD = max(a)/5;
//  0.25<max(a)<=0.5:   FOD_THD = 0.1;
//       max(a)< 0.25:  FOD_THD = max(0.05,max(a)/2.5);
//       max(a)< 0.125: FOD_THD = 0.05;

// Calculates the peak threshold and target number of peaks
// Returns the peak threshold
float TranShi2015::getTargetPeakCountAndThresh(const std::vector<float>& peakValues, int &targetPeakCount)
{
    float FOD_THD = std::min(std::max(0.05, peakValues[0] / 2.5), std::max(0.1, peakValues[0] / 5.0));

    int numBigPeaks     = 0;
    int numSmallPeaks   = 0;

    for(auto i : peakValues) {
        if (i >= FOD_THD)  numBigPeaks++;
        else if (i > 0.05) numSmallPeaks++;
    }

    targetPeakCount = (numBigPeaks >= maxCrossings) ? maxCrossings : numBigPeaks + std::min(1, numSmallPeaks);
    
    return (targetPeakCount >= 1) ? std::max(0.05f, peakValues[targetPeakCount - 1]*0.5f) : 0.05;
}
