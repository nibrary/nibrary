#pragma once

#include "base/config.h"

namespace NIBR {

class Trekker {

public:
    
    Trekker();
    Trekker(std::string    a);
    ~Trekker();

    void        run();      // Runs the tracker and populates TRACKER::tractogram
    void        reset();    // Clears TRACKER::tractogram. run() will generate a new TRACKER::tractogram with the same settings.


    // General options
    void        numberOfThreads(int n);
    void        runTimeLimit(int t);
    void        idleTimeLimit(int t);

    // Seeding options
    void        seed_clear();   // All seed options are deleted

    void        seed_count(long count);
    void        seed_density(double density);
    void        seed_trials(int n);

    void        seed_surface_faceDensity(std::string surf_faceDensity_fname);
    void        seed_surface_vertDensity(std::string surf_vertDensity_fname);
    void        seed_surface_fieldDensity(std::string surf_fieldDensity);
    void        seed_surface_density_fileDataType(std::string  densityFileDataType);
    void        seed_surface_useNormForDir(bool useNorm);
    void        seed_surface_useInside(bool useInside);

    // Pathway options
    bool        pathway_addSeed(std::vector<std::string> rule);
    bool        pathway_addRule(std::vector<std::string> rule);
    bool        pathway_remove(int ruleInd);
    
    bool        pathway_minLength(double len);
    bool        pathway_maxLength(double len);
    bool        pathway_stopAtMax(bool q);
    bool        pathway_oneSided (bool q);
    bool        pathway_skipSeed (bool q);
    bool        pathway_inOrder  (bool q);
    bool        pathway_noEdgeSeed(bool q);
    
    // Algorithm options
    void        algorithm_clear();    // Everything about the algorithm is deleted. All algorithm parameters can always be changed without using this function.

    // FOD options
    void        fod(std::string img_FOD_path);
    void        fodSphere(std::string fod_sphere_path);
    void        fodIsSym(bool fodIsSym);
    void        fodDiscretization(bool fodDiscretization);
    void        orderOfDirections(std::string orderOfDirectionsTextInput);

    // Tracking options
    void        stepSize(double stepSize);
    void        minRadiusOfCurvature(double minRadiusOfCurvature);
    void        minDataSupport(double minDataSupport);
    void        dataSupportExponent(double dataSupportExponent);
    void        ignoreWeakLinks(double weakLinkThresh);

    // Sampling options
    void        maxEstInterval(int maxEstInterval);
    void        initMaxEstTrials(int initMaxEstTrials);
    void        propMaxEstTrials(int propMaxEstTrials);
    void        maxSamplingPerStep(int triesPerRejectionSampling);
    void        useBestAtInit(bool useBestAtInit);
    void        useLegacySampling(bool useLegacySampling);
    void        samplingQuality(int samplingQuality);

    // Probe options
    void        probeCount(int probeCount);
    void        probeRadius(double probeRadius);
    void        probeLength(double probeLength);
    void        probeQuality(int probeQuality);

    // Output options
    void        writeStepSize(float _writeStepSize);
    void        saveFrame(bool saveFrame);
    
};

}
