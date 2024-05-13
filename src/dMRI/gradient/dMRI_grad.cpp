#include "dMRI_grad.h"
#include "base/verbose.h"
#include "math/rotation.h"

using namespace NIBR;

// REFACTOR: This function is missing the part where bvals_fname and bvecs_fname should be read from text files
std::tuple<bool,std::vector<std::vector<float>> > NIBR::readGradTable(std::string bvals_fname, std::string bvecs_fname)
{

    std::vector<float> bvals;
    std::vector<std::vector<float>> bvecs;
    std::vector<std::vector<float>> gradTable;

    bool success = true;

    // REFACTOR: Below is the missing part that should fill bvals and bvecs
    // Read the files here if anything goes wrong set success = false;
    
    bvals = readBvalues(bvals_fname, &success);
    bvecs = readBvectors(bvecs_fname, &success);

    success = (bvecs.size() == bvals.size());

    if (!success) {
        return std::make_tuple(false,gradTable);
    }

    for (size_t n = 0; n < bvals.size(); n++) {

        if (norm(bvecs[n])>EPS4) {
            normalize(bvecs[n]);
        } else {
            bvecs[n][0] = 0.0f;
            bvecs[n][1] = 0.0f;
            bvecs[n][2] = 0.0f;
        }

        std::vector<float> grad = bvecs[n];
        grad.push_back(bvals[n]);
        gradTable.push_back(grad);
    }

    return std::make_tuple(true,gradTable);
    

}

void NIBR::rotateGradTable(std::vector<std::vector<float>>* gradTable, float ijk2xyz[][4])
{
    for (std::vector<float>& g : *gradTable) {
        rotate(g,ijk2xyz);
        if (norm(g)>EPS8) {
            normalize(g);
        } else {
            g[0] = g[1] = g[2] = 0.0f;
        }
    }
}

// Reads bvalues
std::vector<float> NIBR::readBvalues(std::string bvals_fname, bool *success)
{
    std::vector<float> bvalues;
    std::ifstream bvalStream(bvals_fname);
    if (!bvalStream) {
        *success = false;
        NIBR::disp(MSG_ERROR,"Failed to read bvalues file.");
        return bvalues;
    }
    float bval;
    while(bvalStream >> bval) {
        bvalues.push_back(bval);
    }
    
    bvalStream.close();
    return bvalues;
}

// Reads bvectors
std::vector<std::vector<float>> NIBR::readBvectors(std::string bvecs_fname, bool *success)
{
    std::vector<std::vector<float>> bvectors;
    std::ifstream bvecStream(bvecs_fname);
    std::vector<float> bvecs_1D;
    if (!bvecStream) {
        *success = false;
        NIBR::disp(MSG_ERROR,"Failed to read bvectors file.");
        return bvectors;
    }
    float bvec;
    while(bvecStream >> bvec) {
        bvecs_1D.push_back(bvec);
    }
    bvecStream.close();
    if (bvecs_1D.size() % 3 != 0) {
        *success = false;
        NIBR::disp(MSG_ERROR,"Wrong Bvector data dimensions.");
        return bvectors;
    }
    int third = bvecs_1D.size() / 3;
    for(int i = 0; i < third; i++) {
        bvectors.push_back({bvecs_1D[i], bvecs_1D[i + third], bvecs_1D[i + 2 * third]});
    }
    return bvectors;
}