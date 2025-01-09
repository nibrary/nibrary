#include "trekker.h"
// #include "algorithms/ptt/algorithm_ptt_params.h"
#include "pathway/pathway.h"
#include "pathway/parseRule.h"
#include "tracker/tracker.h"
#include "tracker/trackerThread.h"
#include <climits>
#include <cstdint>
#include <ostream>
#include <string>

using namespace NIBR;

Trekker::Trekker()
{
    TRACKER::algorithm      = PTT;
    TRACKER::nThreads       = MT::MAXNUMBEROFTHREADS();
    TRACKER::runTimeLimit   = 0;
    TRACKER::idleTimeLimit  = 0;
    TRACKER::pw.enableTracking(true);
}

Trekker::Trekker(std::string a)
{
    if (a=="ptt") TRACKER::algorithm = PTT;

    TRACKER::nThreads  = MT::MAXNUMBEROFTHREADS();
    TRACKER::runTimeLimit   = 0;
    TRACKER::idleTimeLimit  = 0;
    TRACKER::pw.enableTracking(true);
}

Trekker::~Trekker()
{
    disp(MSG_DEBUG,"Cleaning up trekker");
    seed_clear();
    disp(MSG_DEBUG,"Seed cleaned");
    algorithm_clear();
    disp(MSG_DEBUG,"Algorithm cleaned");
    reset();
    disp(MSG_DEBUG,"Tracker cleaned");
}

// General options
void Trekker::numberOfThreads(int n)      {TRACKER::nThreads      = (n>0)  ? n : 1;}
void Trekker::runTimeLimit(int t)         {TRACKER::runTimeLimit  = (t>=0) ? t : 0;}
void Trekker::idleTimeLimit(int t)        {TRACKER::idleTimeLimit = (t>=0) ? t : 0;}

// Seeding options
void Trekker::seed_clear()  {TRACKER::seed.clear();}

void Trekker::seed_count(long count) {long tmp = (count>0) ? count : INT32_MAX; TRACKER::seed.setSeedCount(tmp);}
void Trekker::seed_density(double density) {double tmp = (density>0) ? density : 1; if (tmp!=density) disp(MSG_WARN,"Density is forced to be 1 since it has to be positive.");TRACKER::seed.setSeedDensity(density);}
void Trekker::seed_trials(int n) {TRACKER::seed.trials = (n>1) ? n : 1;}

void Trekker::seed_surface_faceDensity(std::string surf_faceDensity_fname) {TRACKER::seed.surf_faceDensity_filename = surf_faceDensity_fname;}
void Trekker::seed_surface_vertDensity(std::string surf_vertDensity_fname) {TRACKER::seed.surf_vertDensity_filename = surf_vertDensity_fname;}
void Trekker::seed_surface_fieldDensity(std::string surf_fieldDensity) {TRACKER::seed.surf_density_fieldname = surf_fieldDensity;}
void Trekker::seed_surface_density_fileDataType(std::string  densityFileDataType) {TRACKER::seed.surf_densityFile_dataType = densityFileDataType;}
void Trekker::seed_surface_useNormForDir(bool useNorm) {TRACKER::seed.useSurfNorm = useNorm;}

// Pathway options
bool Trekker::pathway_addSeed(std::vector<std::string> rule) { 

    if (TRACKER::pw.hasOneSeed()) {
        TRACKER::pw.remove(TRACKER::pw.getSeedRuleInd());
    } 
    
    auto r = parseSeedInput(rule);
    return (r.src != undef_src) ? TRACKER::seed.setSeed(r) : false;
}


bool Trekker::pathway_addRule(std::vector<std::string> rule) {
    
    auto rules = parsePathwayInput(rule); 
    
    for (auto r : rules) { 
        if(!TRACKER::pw.add(r)) 
            return false;
    } 
    
    return true;
}

bool Trekker::pathway_remove(int ruleInd)                   {return TRACKER::pw.remove(ruleInd);}
bool Trekker::pathway_minLength(double len)                 {return TRACKER::pw.setMinLength(len);}
bool Trekker::pathway_maxLength(double len)                 {return TRACKER::pw.setMaxLength(len);}
bool Trekker::pathway_stopAtMax(bool q)                     {return TRACKER::pw.stopAtMax(q);}
bool Trekker::pathway_oneSided (bool q)                     {return TRACKER::pw.oneSided(q);}
bool Trekker::pathway_inOrder  (bool q)                     {return TRACKER::pw.inOrder(q);}
bool Trekker::pathway_skipSeed (bool q)                     {return TRACKER::pw.skipSeed(q);}
// bool Trekker::pathway_noEdgeSeed (bool q)                   {return TRACKER::pw.noEdgeSeed(q);}

// Tracking options
void Trekker::algorithm_clear() {
     switch (TRACKER::algorithm) {
        case PTT:
            TRACKER::params_ptt.clear();
            break;

        default:
            break;
    }
}


// FOD options

void Trekker::fod(std::string img_FOD_path) {
    switch (TRACKER::algorithm) {
        case PTT:
            TRACKER::params_ptt.img_FOD_path = img_FOD_path;
            TRACKER::params_ptt.needsUpdate();
            break;

        default:
            break;
    }
}


void Trekker::fodSphere(std::string fod_sphere_path) {
    switch (TRACKER::algorithm) {
        case PTT:
            TRACKER::params_ptt.fod_sphere_path = fod_sphere_path;
            TRACKER::params_ptt.needsUpdate();
            break;

        default:
            break;
    }
}


void Trekker::fodIsSym(bool fodIsSym) {
    switch (TRACKER::algorithm) {
        case PTT:
            TRACKER::params_ptt.fodIsSym = fodIsSym;
            TRACKER::params_ptt.needsUpdate();
            break;

        default:
            break;
    }
}


void Trekker::fodDiscretization(bool fodDiscretization) {
    switch (TRACKER::algorithm) {
        case PTT:
            TRACKER::params_ptt.fodDiscretization = fodDiscretization;
            TRACKER::params_ptt.needsUpdate();
            break;

        default:
            break;
    }
}


void Trekker::orderOfDirections(std::string orderOfDirectionsTextInput) {
    switch (TRACKER::algorithm) {
        case PTT:
            TRACKER::params_ptt.orderOfDirectionsTextInput = orderOfDirectionsTextInput;
            TRACKER::params_ptt.needsUpdate();
            break;

        default:
            break;
    }
}


// Tracking options

std::pair<double, std::string> parseParamStr(std::string paramStr) {

    std::vector<std::string> val = splitString(paramStr,',');

    if ( (val.size() < 1) && (val.size() > 2) ) {
        disp(MSG_FATAL, "Unexpected tracking parameter: %s", paramStr.c_str());
        return std::make_pair(std::numeric_limits<double>::quiet_NaN(), "");
    }

    double v = NIBR::isNumber(val[0]);

    if ((val.size() == 2) && (isnan(v)) ) {
        disp(MSG_FATAL, "First value is not numeric: %s", paramStr.c_str());
        return std::make_pair(std::numeric_limits<double>::quiet_NaN(), "");
    }

    if (val.size() == 2) 
        return std::make_pair(v,val[1]);

    if ( (val.size() == 1) && (!isnan(v)) )
        return std::make_pair(v,"");
    
    // if ( (val.size() == 1) &&  (isnan(v)) )
    
    return std::make_pair(v,val[0]);

}

void Trekker::writeStepSize(float _writeStepSize) {
    switch (TRACKER::algorithm) {
        case PTT:
            TRACKER::params_ptt.outputStep_global = _writeStepSize;
            TRACKER::params_ptt.needsUpdate();
            break;

        default:
            break;
    }
}

void Trekker::writeStepSize(std::string paramStr) {

    auto pp = parseParamStr(paramStr);

    switch (TRACKER::algorithm) {
        case PTT:
            if (!isnan(pp.first)) TRACKER::params_ptt.outputStep_global = pp.first;
            TRACKER::params_ptt.outputStep_img_path = pp.second;
            TRACKER::params_ptt.needsUpdate();
            break;

        default:
            break;
    }

}

void Trekker::stepSize(double stepSize) {
    switch (TRACKER::algorithm) {
        case PTT:
            TRACKER::params_ptt.stepSize_global = stepSize;
            TRACKER::params_ptt.needsUpdate();
            break;

        default:
            break;
    }
}

void Trekker::stepSize(std::string paramStr) {

    auto pp = parseParamStr(paramStr);

    switch (TRACKER::algorithm) {
        case PTT:
            if (!isnan(pp.first)) TRACKER::params_ptt.stepSize_global = pp.first;
            TRACKER::params_ptt.stepSize_img_path = pp.second;
            TRACKER::params_ptt.needsUpdate();
            break;

        default:
            break;
    }
}


void Trekker::minRadiusOfCurvature(double minRadiusOfCurvature) {
    switch (TRACKER::algorithm) {
        case PTT:
            TRACKER::params_ptt.minRadiusOfCurvature_global = minRadiusOfCurvature;
            TRACKER::params_ptt.needsUpdate();
            break;

        default:
            break;
    }
}

void Trekker::minRadiusOfCurvature(std::string paramStr) {

    auto pp = parseParamStr(paramStr);

    switch (TRACKER::algorithm) {
        case PTT:
            if (!isnan(pp.first)) TRACKER::params_ptt.minRadiusOfCurvature_global = pp.first;
            TRACKER::params_ptt.minRadiusOfCurvature_img_path = pp.second;
            TRACKER::params_ptt.needsUpdate();
            break;

        default:
            break;
    }
}


void Trekker::minDataSupport(double minDataSupport) {
    switch (TRACKER::algorithm) {
        case PTT:
            TRACKER::params_ptt.minDataSupport_global = minDataSupport;
            TRACKER::params_ptt.needsUpdate();
            break;

        default:
            break;
    }
}

void Trekker::minDataSupport(std::string paramStr) {

    auto pp = parseParamStr(paramStr);

    switch (TRACKER::algorithm) {
        case PTT:
            if (!isnan(pp.first)) TRACKER::params_ptt.minDataSupport_global = pp.first;
            TRACKER::params_ptt.minDataSupport_img_path = pp.second;
            TRACKER::params_ptt.needsUpdate();
            break;

        default:
            break;
    }

}


void Trekker::dataSupportExponent(double dataSupportExponent) {
    switch (TRACKER::algorithm) {
        case PTT:
            TRACKER::params_ptt.dataSupportExponent_global = dataSupportExponent;
            TRACKER::params_ptt.needsUpdate();
            break;

        default:
            break;
    }
}

void Trekker::dataSupportExponent(std::string paramStr) {
    
    auto pp = parseParamStr(paramStr);

    switch (TRACKER::algorithm) {
        case PTT:
            if (!isnan(pp.first)) TRACKER::params_ptt.dataSupportExponent_global = pp.first;
            TRACKER::params_ptt.dataSupportExponent_img_path = pp.second;
            TRACKER::params_ptt.needsUpdate();
            break;

        default:
            break;
    }

}


void Trekker::ignoreWeakLinks(double weakLinkThresh) {
    switch (TRACKER::algorithm) {
        case PTT:
            TRACKER::params_ptt.weakLinkThresh = weakLinkThresh;
            TRACKER::params_ptt.needsUpdate();
            break;

        default:
            break;
    }
}

void Trekker::paramImgMask(std::string param_mask_fname) {
    switch (TRACKER::algorithm) {
        case PTT:
            TRACKER::params_ptt.img_param_mask_path = param_mask_fname;
            TRACKER::params_ptt.needsUpdate();
            break;

        default:
            break;
    }
}

void Trekker::saveSeedIndex(bool q) {
    TRACKER::saveSeedIndex = q;
}


// Sampling options
void Trekker::useBestAtInit(bool useBestAtInit) {
    switch (TRACKER::algorithm) {
        case PTT:
            TRACKER::params_ptt.useBestAtInit = useBestAtInit;
            TRACKER::params_ptt.needsUpdate();
            break;

        default:
            break;
    }
}

void Trekker::useLegacySampling(bool useLegacySampling) {
    switch (TRACKER::algorithm) {
        case PTT:
            TRACKER::params_ptt.useLegacySampling = useLegacySampling;
            TRACKER::params_ptt.needsUpdate();
            break;

        default:
            break;
    }
}

void Trekker::samplingQuality(int samplingQuality) {
    switch (TRACKER::algorithm) {
        case PTT:
            TRACKER::params_ptt.samplingQuality = samplingQuality;
            TRACKER::params_ptt.needsUpdate();
            break;

        default:
            break;
    }
}


void Trekker::maxEstInterval(int maxEstInterval) {
    switch (TRACKER::algorithm) {
        case PTT:
            TRACKER::params_ptt.maxEstInterval_global = maxEstInterval;
            TRACKER::params_ptt.needsUpdate();
            break;

        default:
            break;
    }
}

void Trekker::maxEstInterval(std::string paramStr) {
    
    auto pp = parseParamStr(paramStr);

    switch (TRACKER::algorithm) {
        case PTT:
            if (!isnan(pp.first)) TRACKER::params_ptt.maxEstInterval_global = pp.first;
            TRACKER::params_ptt.maxEstInterval_img_path = pp.second;
            TRACKER::params_ptt.needsUpdate();
            break;

        default:
            break;
    }

}


void Trekker::initMaxEstTrials(int initMaxEstTrials) {
    switch (TRACKER::algorithm) {
        case PTT:
            TRACKER::params_ptt.initMaxEstTrials_global = initMaxEstTrials;
            TRACKER::params_ptt.needsUpdate();
            break;

        default:
            break;
    }
}

void Trekker::initMaxEstTrials(std::string paramStr) {
    
    auto pp = parseParamStr(paramStr);

    switch (TRACKER::algorithm) {
        case PTT:
            if (!isnan(pp.first)) TRACKER::params_ptt.initMaxEstTrials_global = pp.first;
            TRACKER::params_ptt.initMaxEstTrials_img_path = pp.second;
            TRACKER::params_ptt.needsUpdate();
            break;

        default:
            break;
    }

}


void Trekker::propMaxEstTrials(int propMaxEstTrials) {
    switch (TRACKER::algorithm) {
        case PTT:
            TRACKER::params_ptt.propMaxEstTrials_global = propMaxEstTrials;
            TRACKER::params_ptt.needsUpdate();
            break;

        default:
            break;
    }
}

void Trekker::propMaxEstTrials(std::string paramStr) {

    auto pp = parseParamStr(paramStr);

    switch (TRACKER::algorithm) {
        case PTT:
            if (!isnan(pp.first)) TRACKER::params_ptt.propMaxEstTrials_global = pp.first;
            TRACKER::params_ptt.propMaxEstTrials_img_path = pp.second;
            TRACKER::params_ptt.needsUpdate();
            break;

        default:
            break;
    }

}

void Trekker::maxSamplingPerStep(int triesPerRejectionSampling) {
    switch (TRACKER::algorithm) {
        case PTT:
            TRACKER::params_ptt.triesPerRejectionSampling_global = triesPerRejectionSampling;
            TRACKER::params_ptt.needsUpdate();
            break;

        default:
            break;
    }
}

void Trekker::maxSamplingPerStep(std::string paramStr) {
    
    auto pp = parseParamStr(paramStr);

    switch (TRACKER::algorithm) {
        case PTT:
            if (!isnan(pp.first)) TRACKER::params_ptt.triesPerRejectionSampling_global = pp.first;
            TRACKER::params_ptt.triesPerRejectionSampling_img_path = pp.second;
            TRACKER::params_ptt.needsUpdate();
            break;

        default:
            break;
    }

}

// Probe options
void Trekker::probeCount(int probeCount) {
    switch (TRACKER::algorithm) {
        case PTT:
            TRACKER::params_ptt.probeCount_global = probeCount;
            TRACKER::params_ptt.needsUpdate();
            break;

        default:
            break;
    }
}

void Trekker::probeCount(std::string paramStr) {
    
    auto pp = parseParamStr(paramStr);

    switch (TRACKER::algorithm) {
        case PTT:
            if (!isnan(pp.first)) TRACKER::params_ptt.probeCount_global = pp.first;
            TRACKER::params_ptt.probeCount_img_path = pp.second;
            TRACKER::params_ptt.needsUpdate();
            break;

        default:
            break;
    }

}


void Trekker::probeRadius(double probeRadius) {
    switch (TRACKER::algorithm) {
        case PTT:
            TRACKER::params_ptt.probeRadius_global = probeRadius;
            TRACKER::params_ptt.needsUpdate();
            break;

        default:
            break;
    }
}

void Trekker::probeRadius(std::string paramStr) {
    
    auto pp = parseParamStr(paramStr);

    switch (TRACKER::algorithm) {
        case PTT:
            if (!isnan(pp.first)) TRACKER::params_ptt.probeRadius_global = pp.first;
            TRACKER::params_ptt.probeRadius_img_path = pp.second;
            TRACKER::params_ptt.needsUpdate();
            break;

        default:
            break;
    }

}


void Trekker::probeLength(double probeLength) {
    switch (TRACKER::algorithm) {
        case PTT:
            TRACKER::params_ptt.probeLength_global = probeLength;
            TRACKER::params_ptt.needsUpdate();
            break;

        default:
            break;
    }
}

void Trekker::probeLength(std::string paramStr) {
    
    auto pp = parseParamStr(paramStr);

    switch (TRACKER::algorithm) {
        case PTT:
            if (!isnan(pp.first)) TRACKER::params_ptt.probeLength_global = pp.first;
            TRACKER::params_ptt.probeLength_img_path = pp.second;
            TRACKER::params_ptt.needsUpdate();
            break;

        default:
            break;
    }

}


void Trekker::probeQuality(int probeQuality) {
    switch (TRACKER::algorithm) {
        case PTT:
            TRACKER::params_ptt.probeQuality_global = probeQuality;
            TRACKER::params_ptt.needsUpdate();
            break;

        default:
            break;
    }
}

void Trekker::probeQuality(std::string paramStr) {
    
    auto pp = parseParamStr(paramStr);

    switch (TRACKER::algorithm) {
        case PTT:
            if (!isnan(pp.first)) TRACKER::params_ptt.probeQuality_global = pp.first;
            TRACKER::params_ptt.probeQuality_img_path = pp.second;
            TRACKER::params_ptt.needsUpdate();
            break;

        default:
            break;
    }

}

// Reset
void Trekker::reset() { TRACKER::reset();}

// Runner
void Trekker::run() {

    // Verify seed count and run time limit requirement
    if ( ((TRACKER::seed.sCount == INT64_MAX) || (TRACKER::seed.sCount <= 0)) && ((TRACKER::runTimeLimit == INT32_MAX) || (TRACKER::runTimeLimit <= 0)) ) {
        disp(MSG_ERROR, "Both seed count and run time limit can't be infinite");
        return;
    }

    // Reset tracker. This restarts the timer, seed counter and logger, allowing the run() to be called again
    disp(MSG_DETAIL, "Resetting tracker");
    this->reset();

    // Print general options
    TRACKER::print();

    // Check seed
    disp(MSG_DETAIL, "Updating seed");
    if (!TRACKER::seed.update()) {
        disp(MSG_ERROR, "Initialiation of seed failed");
        return;
    }
    TRACKER::seed.print();    

    // Check pathway
    disp(MSG_DETAIL, "Updating pathway rules");
    if (!TRACKER::pw.verify()) {
        disp(MSG_ERROR, "Verification of pathway rules failed");
        return;
    }

    // Check algorithm
    disp(MSG_DETAIL, "Updating PTT");
    if (!TRACKER::params_ptt.update()) {
        disp(MSG_ERROR, "Initialiation of PTT failed");
        return;
    }

    TRACKER::params_ptt.print();
    TRACKER::pw.print();
    
    // Run
    auto runTrekker = [&](NIBR::MT::TASK task)->bool {
        TrackingThread tracker(task.no);
        return tracker.track();        
    };

    
    // Print progress
    std::thread coreThread([&]() {
        NIBR::MT::MTRUN(TRACKER::seed.getMaxTrackerCount(), runTrekker,TRACKER::seed.getCount());
    });

    // Hide message below VERBOSE_INFO

    if (NIBR::VERBOSE()==VERBOSE_INFO) {

        std::string preamble = "\033[1;32mNIBRARY::INFO: \033[0;32m";

        size_t success = 0;
        size_t discard = 0;
        size_t fail    = 0;

        // Skip one line for readability
        std::cout << std::endl;
        int lineCount = 0;

        auto dispCount = [&](size_t& cumCnt, size_t cnt, std::string reason)->void {
            if (cnt > 0) {
                size_t spaces = 34 - reason.length();
                reason = "   " + reason + std::string(spaces,' ') + ": ";
                std::cout << preamble << reason << cnt << "\033[0m" << std::endl; 
                lineCount++;
                cumCnt += cnt;
            }
        };

        auto printReport = [&]()->void {

            std::cout << "\r\033[K" << std::flush;           // Clear the current line
            for (int n=0; n<(lineCount-1); n++)
                std::cout << "\033[A\r\033[K" << std::flush; // Clear previous lines

            lineCount = success = discard = fail = 0;

            std::cout << preamble << "Success report" << "\033[0m" << std::endl; lineCount++;
            dispCount(success,  TRACKER::log_success_REACHED_MAXLENGTH_LIMIT.load(),              "Reached max length");
            dispCount(success,  TRACKER::log_success_REACHED_MINDATASUPPORT_LIMIT.load(),         "Reached min data support");
            dispCount(success,  TRACKER::log_success_SATISFIED_PATHWAY_RULES.load(),              "Satisfied pathway rules");

            std::cout << preamble << "Discard report" << "\033[0m" << std::endl; lineCount++;
            dispCount(discard,  TRACKER::log_discard_TOO_SHORT.load(),                            "Too short");
            dispCount(discard,  TRACKER::log_discard_TOO_LONG.load(),                             "Too long");
            dispCount(discard,  TRACKER::log_discard_DISCARD_ROI_REACHED.load(),                  "Reached discard region");
            dispCount(discard,  TRACKER::log_discard_REQUIRED_ROI_NOT_MET.load(),                 "Required region not found");
            dispCount(discard,  TRACKER::log_discard_REQUIRED_ROI_ORDER_NOT_MET.load(),           "Required order not satisfied");
            dispCount(discard,  TRACKER::log_discard_CANT_MEET_STOP_CONDITION.load(),             "Can't meet stop condition");
            dispCount(discard,  TRACKER::log_discard_ENDED_INSIDE_DISCARD_ROI.load(),             "Ended inside discard region");
            dispCount(discard,  TRACKER::log_discard_REACHED_TIME_LIMIT.load(),                   "Reached time limit");

            std::cout << preamble << "Fail report" << "\033[0m" << std::endl; lineCount++;
            dispCount(fail,     TRACKER::log_failed_BY_THE_ALGORITHM_AT_INITIALIZATION.load(),    "Initialization at the seed failed");
            dispCount(fail,     TRACKER::log_failed_BY_THE_ALGORITHM.load(),                      "Algorithm failed to propagate");
            dispCount(fail,     TRACKER::log_failed_UNKNOWN_REASON.load(),                        "Unknown reason");

            std::cout << std::endl << preamble;
            std::cout << "Success: "  << success              << "      ";
            std::cout << "Discard: "  << discard              << "      ";
            std::cout << "Fail: "     << fail                 << "      ";
            std::cout << "Total: "    << success+discard+fail << "      ";
            std::cout << "Duration: " << TRACKER::runTime()   << " sec";
            std::cout << "\033[0m"    << std::flush;
            lineCount=lineCount+2;

        };


        NIBR::MT::FINISHEDTASKCOUNT() = 0;
        while ((!TRACKER::countIsReached) && (TRACKER::seed.getMaxTrackerCount()>long(MT::FINISHEDTASKCOUNT().load())) )
        {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));

            printReport();

            if(!TRACKER::isWithinTimeLimits()) MT::FINISHEDTASKCOUNT().store(TRACKER::seed.getMaxTrackerCount());
    
        }

        std::cout << std::endl;

        if (TRACKER::idletimeLimitReached) std::cout << preamble << "Reached idle time limit." << "\033[0m" << std::endl << std::flush;
        if (TRACKER::runtimeLimitReached)  std::cout << preamble << "Reached run time limit."  << "\033[0m" << std::endl << std::flush;

    }

    coreThread.join();
    
}
