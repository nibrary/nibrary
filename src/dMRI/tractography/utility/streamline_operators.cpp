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

float NIBR::getStreamlineLength(Streamline& inp)
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

float NIBR::getStreamlineLength_hermite(Streamline& inp, int div)
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

std::vector<float> NIBR::getStreamlineLength_hermiteWithSpeed(const Streamline& inp, int div)
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

std::vector<Point3D> NIBR::colorStreamline(const Streamline& inp)
{

    int len = int(inp.size());

    std::vector<Point3D> colors;
    colors.reserve(len);

    if (len < 1) {
        return colors;
    }

    if (len == 1) {
        colors.push_back({0.0f,0.0f,0.0f});
        return colors;
    }

    for (int n = 0; n < (len-1); n++) {
        Point3D dir;
        vec3sub(dir,inp[n+1],inp[n]);
        normalize(dir);
        vec3abs(dir);
        colors.push_back(std::move(dir));
    }
    
    colors.push_back(colors.back());

    return colors;

}



// Function to calculate the one-sided Hausdorff distance from trk1 to trk2
float oneSidedHausdorffDistance(const Streamline& trk1, const Streamline& trk2) {
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
float NIBR::getHausdorffDistance(Streamline& trk1, Streamline& trk2) {
    float hd1 = oneSidedHausdorffDistance(trk1, trk2);
    float hd2 = oneSidedHausdorffDistance(trk2, trk1);
    return std::max(hd1, hd2);
}

// Minimum average direct-flip (MDF) distance
float NIBR::getMDFDistance(Streamline& trk1, Streamline& trk2) {

    if (trk1.size() != trk2.size()) {
        disp(MSG_ERROR, "Number of points along streamlines must be equal.");
        return NAN;
    }

    float directDist  = 0;
    float flippedDist = 0;
    std::size_t n = trk1.size();

    for (std::size_t i = 0; i < n; ++i) {
        directDist  += dist(trk1[i], trk2[i]);
        flippedDist += dist(trk1[i], trk2[n - i - 1]);
    }

    return std::min(directDist, flippedDist) / float(n);

}

bool NIBR::areStreamlinesIdentical(const Streamline& s1, const Streamline& s2, float tolerance) {
    if (s1.size() != s2.size()) {
        return false;
    }
    if (s1.empty()) {
        return true;
    }
    for (size_t i = 0; i < s1.size(); ++i) {
        const auto& p1 = s1[i];
        const auto& p2 = s2[i];
        if ( (std::fabs(p1[0] - p2[0]) > tolerance) ||
             (std::fabs(p1[1] - p2[1]) > tolerance) ||
             (std::fabs(p1[2] - p2[2]) > tolerance) ) {
            return false;
        }
    }
    return true;
}