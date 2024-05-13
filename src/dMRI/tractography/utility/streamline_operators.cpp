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