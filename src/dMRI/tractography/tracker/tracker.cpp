#include "tracker.h"
#include <climits>
#include <cstdint>

using namespace NIBR;

// User parameters
Algorithm_Type  TRACKER::algorithm;
int             TRACKER::nThreads;
int             TRACKER::runTimeLimit;
int             TRACKER::idleTimeLimit;
Seed            TRACKER::seed;
bool            TRACKER::saveSeedIndex;
Pathway         TRACKER::pw;
Params_PTT      TRACKER::params_ptt;

// Derived parameters
Tractogram                              TRACKER::tractogram;
std::vector<int>                        TRACKER::streamlineLength;
std::vector<int>                        TRACKER::seedIndex;
TractogramField                         TRACKER::seedIndexField;
std::chrono::steady_clock::time_point   TRACKER::initTime;
std::chrono::steady_clock::time_point   TRACKER::lastTime;
bool                                    TRACKER::runtimeLimitReached;
bool                                    TRACKER::idletimeLimitReached;
std::size_t                             TRACKER::currentCount;
bool                                    TRACKER::countIsReached;
int                                     TRACKER::ready_thread_id;
std::mutex                              TRACKER::trackKeeper;

// Loggers
std::atomic<size_t> TRACKER::log_success_REACHED_MAXLENGTH_LIMIT;
std::atomic<size_t> TRACKER::log_success_REACHED_MINDATASUPPORT_LIMIT;
std::atomic<size_t> TRACKER::log_success_SATISFIED_PATHWAY_RULES;
std::atomic<size_t> TRACKER::log_discard_TOO_SHORT;
std::atomic<size_t> TRACKER::log_discard_TOO_LONG;
std::atomic<size_t> TRACKER::log_discard_DISCARD_ROI_REACHED;
std::atomic<size_t> TRACKER::log_discard_REQUIRED_ROI_NOT_MET;
std::atomic<size_t> TRACKER::log_discard_REQUIRED_ROI_ORDER_NOT_MET;
std::atomic<size_t> TRACKER::log_discard_CANT_MEET_STOP_CONDITION;
std::atomic<size_t> TRACKER::log_discard_ENDED_INSIDE_DISCARD_ROI;
std::atomic<size_t> TRACKER::log_discard_REACHED_TIME_LIMIT;
std::atomic<size_t> TRACKER::log_failed_UNKNOWN_REASON;
std::atomic<size_t> TRACKER::log_failed_BY_THE_ALGORITHM_AT_INITIALIZATION;
std::atomic<size_t> TRACKER::log_failed_BY_THE_ALGORITHM;
std::atomic<size_t> TRACKER::log_unexpected_TRACKING_STATUS;

Tractogram& TRACKER::getTractogram() {return TRACKER::tractogram;}

TractogramField& TRACKER::getSeedIndexField() {

    // Clear existing values if any
    if (seedIndexField.data != NULL) {
        int*** toDel = reinterpret_cast<int***>(seedIndexField.data);
        for (size_t s = 0; s < currentCount; s++) {
            for (int l = 0; l < streamlineLength[s]; l++) {
                delete[] toDel[s][l];
            }
            delete[] toDel[s];
        }
        delete[] toDel;
    }

    // Prepare seed index field
    int*** idx = new int**[currentCount];

    for (size_t n = 0; n < currentCount; n++) {
        idx[n]  = new int*[streamlineLength[n]];
        for (int l = 0; l < streamlineLength[n]; l++) {
            idx[n][l]    = new int[1];
            idx[n][l][0] = (int(l) == seedIndex[n]);
        }
    }

    seedIndexField.owner      = POINT_OWNER;
    seedIndexField.name       = "seedIndex";
    seedIndexField.datatype   = INT32_DT;
    seedIndexField.dimension  = 1;
    seedIndexField.data       = reinterpret_cast<void*>(idx);

    return seedIndexField;
}

void TRACKER::reset() {

    // Reset seed
    seed.reset(); // Sets the seed counter to zero
    runTimeLimit  = (runTimeLimit<=0)  ? INT32_MAX             : runTimeLimit;
    idleTimeLimit = (idleTimeLimit<=0) ? DEFAULT_IDLETIMELIMIT : idleTimeLimit;
    
    // Reset seedIndex and seedIndexField
    seedIndex.clear();
    streamlineLength.clear();
    
    if (seedIndexField.data != NULL) {
        int*** toDel = reinterpret_cast<int***>(seedIndexField.data);
        for (size_t s = 0; s < currentCount; s++) {
            for (int l = 0; l < streamlineLength[s]; l++) {
                delete[] toDel[s][l];
            }
            delete[] toDel[s];
        }
        delete[] toDel;
    }   

    // Reset derived parameters
    tractogram.clear();
    initTime               = std::chrono::steady_clock::now();
    lastTime               = initTime;
    runtimeLimitReached    = false;
    currentCount           = 0;
    countIsReached         = false;
    ready_thread_id        = 0;

    // Reset loggers
    log_success_REACHED_MAXLENGTH_LIMIT             = 0;
    log_success_REACHED_MINDATASUPPORT_LIMIT        = 0;
    log_success_SATISFIED_PATHWAY_RULES             = 0;
    log_discard_TOO_SHORT                           = 0;
    log_discard_TOO_LONG                            = 0;
    log_discard_DISCARD_ROI_REACHED                 = 0;
    log_discard_REQUIRED_ROI_NOT_MET                = 0;
    log_discard_REQUIRED_ROI_ORDER_NOT_MET          = 0;
    log_discard_CANT_MEET_STOP_CONDITION            = 0;
    log_discard_ENDED_INSIDE_DISCARD_ROI            = 0;
    log_discard_REACHED_TIME_LIMIT                  = 0;
    log_failed_UNKNOWN_REASON                       = 0;
    log_failed_BY_THE_ALGORITHM_AT_INITIALIZATION   = 0;
    log_failed_BY_THE_ALGORITHM                     = 0;
    log_unexpected_TRACKING_STATUS                  = 0;
}

void TRACKER::print() {

    disp(MSG_INFO,"GENERAL OPTIONS");

    std::cout << "\033[32m";

    std::cout << "numberOfThreads : " << nThreads  << std::endl;
    if (runTimeLimit==INT32_MAX) {
	    std::cout << "runTimeLimit    : infinite " << std::endl;
    } else {
        std::cout << "runTimeLimit    : " << runTimeLimit/60 << " min"<< std::endl;
    }

    std::cout << "idleTimeLimit   : " << idleTimeLimit/60 << " min" << std::endl;

    std::cout << "\033[0m";
}

int TRACKER::runTime() {
    return int(std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now()-TRACKER::initTime).count());
}

int TRACKER::idleTime() {
    return int(std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now()-TRACKER::lastTime).count());
}

bool TRACKER::isWithinTimeLimits() {
    TRACKER::idletimeLimitReached = (TRACKER::idleTime() > TRACKER::idleTimeLimit);
    TRACKER::runtimeLimitReached  = (TRACKER::runTime()  > TRACKER::runTimeLimit );

    return !(TRACKER::idletimeLimitReached||TRACKER::runtimeLimitReached);
}

