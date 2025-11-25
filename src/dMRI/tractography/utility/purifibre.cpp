#include "purifibre.h"
#include "image/sf_image.h"
#include "dMRI/tractography/mappers/tractogram2imageMapper.h"
#include "dMRI/tractography/mappers/gridder_4segmentLength.h"
#include "tractogram_operators.h"


using namespace NIBR;

std::vector<float> NIBR::getFico( TractogramReader* tractogram, 
                                  float trimFactor, 
                                  float voxDim, 
                                  std::tuple<float,int> anisotropicSmoothing, 
                                  float sphericalSmoothing
                                ) 
{

    int N = tractogram->numberOfStreamlines;

    if (N<1) {
        NIBR::disp(MSG_WARN,"Empty tractogram");
        return std::vector<float>();
    }

    // Validate parameters
    trimFactor          = std::clamp(trimFactor,0.0f,50.0f);
    if (voxDim<=0.0f) voxDim = 4.0f;
    voxDim              = std::max(voxDim, float(EPS6));
    sphericalSmoothing  = std::clamp(sphericalSmoothing,0.0f,180.0f);

    std::get<0>(anisotropicSmoothing) = std::max(std::get<0>(anisotropicSmoothing), 0.0f);
    std::get<1>(anisotropicSmoothing) = std::max(std::get<1>(anisotropicSmoothing), 0);

    
    // FICO values
    std::vector<float> fico;
    fico.resize(N);


    // Reset tractogram reader
    tractogram->reset();


    // Compute sTODI
    NIBR::SF::init(true,17);
    NIBR::SF_Image img;
    std::vector<float> bb = getTractogramBBox(tractogram);
    bb.push_back(-0.5);
    bb.push_back(int64_t(SF::getSFCoords().size())-0.5);
    img.createFromBoundingBox(4,bb,voxDim,false);
    disp(MSG_DEBUG,"SF image dims: %lld x %lld x %lld x %lld",img.imgDims[0],img.imgDims[1],img.imgDims[2],img.imgDims[3]);

    Tractogram2ImageMapper<float> gridder(tractogram,&img);
    gridder.anisotropicSmoothing(anisotropicSmoothing);
    allocateGrid_4segmentLength_sf(&gridder);
    gridder.run(processor_4segmentLength_sf<float>, outputCompiler_4segmentLength_sf<float>);
    deallocateGrid_4segmentLength_sf(&gridder);


    // Smooth sTODI image
    img.smooth(sphericalSmoothing);

    
    // Compute SECO and FICO
    tractogram->reset();

    auto run = [&]()->void {

        auto [success,streamline,streamlineId] = tractogram->getNextStreamline();
        
        if (!success) {
            disp(MSG_WARN,"Could not read streamline %lu",streamlineId);
            fico[streamlineId] = 0.0f;
            return;
        }
        
        int len = streamline.size();
        if (len<2) {
            fico[streamlineId] = 0.0f;
            return;
        }

        float T[3];

        int trim = int(len*trimFactor*0.5*0.01);
        if (trim>=int(len/2-1)) trim = int(len/2)-1;
        if (trim<=0) trim  = 0;

        float minSeco   = std::numeric_limits<float>::infinity();

        for (int l=trim; l<(len-trim-1); l++) {
            vec3sub(T,streamline[l+1],streamline[l]);
            normalize(T);

            float seco = img.getSFval(streamline[l].data(),T); // segment-to-bundle coupling (SECO)
            if (seco < minSeco)
                minSeco = seco;
            
        }
        float seco = img.getSFval(streamline[len-trim-1].data(),T);
        if (seco < minSeco)
            minSeco = seco;

        fico[streamlineId] = std::log(minSeco+1);

    };
    NIBR::MT::MTRUN(N, "Computing FICO", run);

    NIBR::SF::clean();

    return fico;

}


std::vector<size_t> NIBR::purify(const std::vector<float>& fico, float puriFactor) {

    int N = fico.size();

    if (N<1) {
        NIBR::disp(MSG_WARN,"Empty FICO values");
        return std::vector<size_t>();
    }

    std::vector<size_t> idx;
    int remN = std::floor(float(fico.size())*puriFactor*0.01);

    // Create index+fico pairs
    std::vector<std::tuple<size_t,float>> ficoWithIdx(N);
    for (size_t n=0; n<ficoWithIdx.size(); n++) {
        ficoWithIdx[n] = std::make_tuple(n,fico[n]);
    }

    // Sort based on fico values
    std::sort(ficoWithIdx.begin(), ficoWithIdx.end(), [](auto a, auto b) {return std::get<1>(a) < std::get<1>(b);} );

    for (int n=remN; n<N; n++)
        idx.push_back(std::get<0>(ficoWithIdx[n]));

    return idx;


}
    