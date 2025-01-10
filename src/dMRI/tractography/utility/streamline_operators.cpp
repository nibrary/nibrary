#include "streamline_operators.h"
#include "math/core.h"
#include "math/interp1D.h"
#include <complex>

using namespace NIBR;

float NIBR::getStreamlineLength(float** inp, int len) 
{
    if (len<2)
        return 0;

    float length = 0;

    for (int i=1; i<len; i++) {
        length += dist(inp[i],inp[i-1]);
    }

    return length;
}

float NIBR::getStreamlineLength(std::vector<std::vector<float>>& inp)
{
    int len = inp.size();

    if (len<2)
        return 0;

    float length = 0;

    for (int i=1; i<len; i++) {
        length += dist(inp[i],inp[i-1]);
    }

    return length;
}

float NIBR::getStreamlineLength(std::vector<Point>& inp) {

    int len = inp.size();

    if (len<2)
        return 0;

    float length = 0;

    for (int i=1; i<len; i++) {
        length += dist(inp[i],inp[i-1]);
    }

    return length;

}

float NIBR::getStreamlineLength_hermite(float** inp, int len, int div) 
{
    if (len<2)
        return 0;

    if (len==2)
        return dist(inp[1],inp[0]);

    float length = 0;
    float T0[3], T1[3], Tn1[3], Tn2[3];

    vec3sub(T0, inp[1], inp[0]);

    for (int i=1; i<len; i++) {
        
        if (i!=(len-1)) {
            vec3sub(Tn1,inp[i],  inp[i-1]);
            vec3sub(Tn2,inp[i+1],inp[i]  );
            vec3add(T1,Tn1,Tn2);
            T1[0] *= 0.5f;
            T1[1] *= 0.5f;
            T1[2] *= 0.5f;
        } else {
            vec3sub(T1,inp[i],inp[i-1]);
        }

        length += hermiteLength(inp[i-1], T0, inp[i], T1, div);

        T0[0] = T1[0];
        T0[1] = T1[1];
        T0[2] = T1[2];
    }

    return length;
}

float NIBR::getStreamlineLength_hermite(std::vector<std::vector<float>>& inp, int div)
{

    int len = int(inp.size());

    if (len<2)
        return 0;

    if (len==2)
        return dist(inp[1],inp[0]);

    float length = 0;
    float T0[3], T1[3], Tn1[3], Tn2[3];

    vec3sub(T0, inp[1], inp[0]);

    for (int i=1; i<len; i++) {
        
        if (i!=(len-1)) {
            vec3sub(Tn1,inp[i],  inp[i-1]);
            vec3sub(Tn2,inp[i+1],inp[i]  );
            vec3add(T1,Tn1,Tn2);
            T1[0] *= 0.5f;
            T1[1] *= 0.5f;
            T1[2] *= 0.5f;
        } else {
            vec3sub(T1,inp[i],inp[i-1]);
        }

        length += hermiteLength(inp[i-1], T0, inp[i], T1, div);

        T0[0] = T1[0];
        T0[1] = T1[1];
        T0[2] = T1[2];
    }

    return length;
}

std::vector<float> NIBR::getStreamlineLength_hermiteWithSpeed(std::vector<std::vector<float>>& inp, int div)
{

    int len = int(inp.size());

    std::vector<float> out;

    if (len<2)
        return std::vector<float>();

    if (len==2) {
        float T[3];
        vec3sub(T, inp[1], inp[0]);
        return hermiteLengthWithSpeed(inp[0], T, inp[1], T, div);
    }

    float length = 0;
    float T0[3], T1[3], Tn1[3], Tn2[3];

    vec3sub(T0, inp[1], inp[0]);

    for (int i=1; i<len; i++) {
        
        if (i!=(len-1)) {
            vec3sub(Tn1,inp[i],  inp[i-1]);
            vec3sub(Tn2,inp[i+1],inp[i]  );
            vec3add(T1,Tn1,Tn2);
            T1[0] *= 0.5f;
            T1[1] *= 0.5f;
            T1[2] *= 0.5f;
        } else {
            vec3sub(T1,inp[i],inp[i-1]);
        }

        auto segCumLen = hermiteLengthWithSpeed(inp[i-1], T0, inp[i], T1, div);

        for (auto it = segCumLen.begin(); it != segCumLen.end(); ++it) {
            *it += length;
        }

        length = segCumLen.back();

        out.insert(out.end(),segCumLen.begin(),segCumLen.end());

        T0[0] = T1[0];
        T0[1] = T1[1];
        T0[2] = T1[2];
    }

    return out;
}

float** NIBR::colorStreamline(std::vector<std::vector<float>>& inp)
{

    int len = int(inp.size());

    float** colors = new float*[len];

    for (int n = 0; n < (len-1); n++) {
        float* dir = new float[3];
        vec3sub(dir,inp[n+1],inp[n]);
        normalize(dir);
        vec3abs(dir);
        colors[n] = dir;
    }
    
    float* dir = new float[3];
    dir[0] = colors[len-2][0];
    dir[1] = colors[len-2][1];
    dir[2] = colors[len-2][2];

    colors[len-1] = dir;

    return colors;

}



// Function to calculate the one-sided Hausdorff distance from trk1 to trk2
float oneSidedHausdorffDistance(const std::vector<std::vector<float>>& trk1, const std::vector<std::vector<float>>& trk2) {
    float maxDist = 0;

    for (const auto& point1 : trk1) {
        float minDist = std::numeric_limits<float>::infinity();
        for (const auto& point2 : trk2) {
            float currentDist = dist(point1, point2);
            if (currentDist < minDist) {
                minDist = currentDist;
            }
        }
        if (minDist > maxDist) {
            maxDist = minDist;
        }
    }

    return maxDist;
}

// Two-sided Hausdorff distance
float NIBR::getHausdorffDistance(std::vector<std::vector<float>>& trk1, std::vector<std::vector<float>>& trk2) {
    float hd1 = oneSidedHausdorffDistance(trk1, trk2);
    float hd2 = oneSidedHausdorffDistance(trk2, trk1);
    return std::max(hd1, hd2);
}

// Minimum average direct-flip (MDF) distance
float NIBR::getMDFDistance(std::vector<std::vector<float>>& trk1, std::vector<std::vector<float>>& trk2) {

    if (trk1.size() != trk2.size()) {
        disp(MSG_ERROR, "Number of points along streamlines must be equal.");
        return NAN;
    }

    float directDist  = 0;
    float flippedDist = 0;
    size_t n = trk1.size();

    for (size_t i = 0; i < n; ++i) {
        directDist  += dist(trk1[i], trk2[i]);
        flippedDist += dist(trk1[i], trk2[n - i - 1]);
    }

    return std::min(directDist, flippedDist) / float(n);

}