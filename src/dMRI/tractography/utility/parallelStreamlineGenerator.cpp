#include "parallelStreamlineGenerator.h"

using namespace NIBR;

void NIBR::getParallelStreamlines(Tractogram& out, NIBR::TractogramReader* tractogram, float radius, int ringCount, int pointCountPerRing)
{
    
    int N   = tractogram->numberOfStreamlines;
    int par = ringCount*pointCountPerRing+1;
    
    // Create output
    if (N<1)
        return;
    else
        out.resize(N*par);

    // Prepare rings
    float radStep = radius/float(ringCount);
    float angStep = TWOPI/float(pointCountPerRing);
    std::vector<float> scale_N;
    std::vector<float> scale_B;
    scale_N.push_back(0.0f);
    scale_B.push_back(0.0f);
    for (float i=1; i<(ringCount+1); i++)
        for (float j=0; j<pointCountPerRing; j++) {
            scale_N.push_back(i*radStep*std::sin(angStep*j));
            scale_B.push_back(i*radStep*std::cos(angStep*j));
        }

    // Iterate throught the whole tractogram
    tractogram->reset();

    auto genParallel = [&](const NIBR::MT::TASK& task)->void{

        auto [success,streamline,streamlineId] = tractogram->getNextStreamline();

        int len = streamline.size();

        if (len<2)
            return;

        for (int p=0; p<par; p++) {
            out[streamlineId*par+p].resize(len);
        }

        // To make cylinders around the streamlines, we will compute the parallel transport frame
        float T[3], N[3], B[3], curr_T[3];
        float R_axis[3], R_angle;
        float R[4][4];

        // Get a random initial PTF and handle the initial point
        vec3sub(T,streamline[1],streamline[0]);
        normalize(T);
        NIBR::MT::RNDM()[task.threadId]->getAUnitRandomPerpVector(N,T);
        normalize(N);
        cross(B,T,N);

        for (int p=0; p<par; p++)
            for (int i=0; i<3; i++)
                out[streamlineId*par+p][0][i] = streamline[0][i] + N[i]*scale_N[p] + + B[i]*scale_B[p];

        for (auto l=1; l<(len-1); l++) {

            vec3sub(curr_T,streamline[l+1],streamline[l]);
            normalize(curr_T);

            R_angle = std::acos(std::clamp(dot(T,curr_T),-1.0,1.0));

            curr_T[0] += EPS4;  // for numerical stability reasons
            curr_T[1] += EPS4;
            curr_T[2] += EPS4;
            
            cross(R_axis,T,curr_T);
            normalize(R_axis);

            axisangle2Rotation(R_axis, R_angle, R);
            rotate(curr_T,T,R); // curr_T is used as a temp var
            rotate(R_axis,N,R); // R_axis is used as a temp var

            T[0] = curr_T[0];
            T[1] = curr_T[1];
            T[2] = curr_T[2];
            normalize(T);

            N[0] = R_axis[0];
            N[1] = R_axis[1];
            N[2] = R_axis[2];
            normalize(N);

            cross(B,T,N);

            for (int p=0; p<par; p++)
                for (int i=0; i<3; i++)
                    out[streamlineId*par+p][l][i] = streamline[l][i] + N[i]*scale_N[p] + + B[i]*scale_B[p];

        }

        // Repeat the last N, B for the last segment
        for (int p=0; p<par; p++)
            for (int i=0; i<3; i++)
                out[streamlineId*par+p][len-1][i] = streamline[len-1][i] + N[i]*scale_N[p] + + B[i]*scale_B[p];        

    };

    NIBR::MT::MTRUN(N, "Generating parallel streamlines", genParallel);

}


void NIBR::getParallelStreamlines(Tractogram& out, NIBR::TractogramReader* tractogram, float sigma, int par)
{
    
    int N = tractogram->numberOfStreamlines;

    // Create output
    if (N<1)
        return;
    else
        out.resize(par*N);

    // Prepare points to track
    std::vector<float> scale_N;
    std::vector<float> scale_B;
    scale_N.push_back(0);
    scale_B.push_back(0);
    for (int n=0; n<par; n++) {
        scale_N.push_back(NIBR::MT::RNDM()[0]->normal_m0_s1()*sigma);
        scale_B.push_back(NIBR::MT::RNDM()[0]->normal_m0_s1()*sigma);
    }

    // Iterate throught the whole tractogram

    tractogram->reset();

    auto genParallel = [&](const NIBR::MT::TASK& task)->void{

        auto [success,streamline,streamlineId] = tractogram->getNextStreamline();

        int len = streamline.size();
        
        if (len<2)
            return;

        for (int p=0; p<par; p++) {
            out[streamlineId*par+p].resize(len);
        }

        // To make cylinders around the streamlines, we will compute the parallel transport frame
        float T[3], N[3], B[3], curr_T[3];
        float R_axis[3], R_angle;
        float R[4][4];

        // Get a random initial PTF and handle the initial point
        vec3sub(T,streamline[1],streamline[0]);
        normalize(T);
        NIBR::MT::RNDM()[task.threadId]->getAUnitRandomPerpVector(N,T);
        normalize(N);
        cross(B,T,N);

        for (int p=0; p<par; p++)
            for (int i=0; i<3; i++)
                out[streamlineId*par+p][0][i] = streamline[0][i] + N[i]*scale_N[p] + + B[i]*scale_B[p];

        for (auto l=1; l<(len-1); l++) {

            vec3sub(curr_T,streamline[l+1],streamline[l]);
            normalize(curr_T);

            R_angle = std::acos(std::clamp(dot(T,curr_T),-1.0,1.0));

            curr_T[0] += EPS4;  // for numerical stability reasons
            curr_T[1] += EPS4;
            curr_T[2] += EPS4;
            
            cross(R_axis,T,curr_T);
            normalize(R_axis);

            axisangle2Rotation(R_axis, R_angle, R);
            rotate(curr_T,T,R); // curr_T is used as a temp var
            rotate(R_axis,N,R); // R_axis is used as a temp var

            T[0] = curr_T[0];
            T[1] = curr_T[1];
            T[2] = curr_T[2];
            normalize(T);

            N[0] = R_axis[0];
            N[1] = R_axis[1];
            N[2] = R_axis[2];
            normalize(N);

            cross(B,T,N);

            for (int p=0; p<par; p++)
                for (int i=0; i<3; i++)
                    out[streamlineId*par+p][l][i] = streamline[l][i] + N[i]*scale_N[p] + + B[i]*scale_B[p];

        }

        // Repeat the last N, B for the last segment
        for (int p=0; p<par; p++)
            for (int i=0; i<3; i++)
                out[streamlineId*par+p][len-1][i] = streamline[len-1][i] + N[i]*scale_N[p] + + B[i]*scale_B[p];

    };

    NIBR::MT::MTRUN(N, "Generating parallel streamlines", genParallel);

}