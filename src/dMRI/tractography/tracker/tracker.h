#pragma once

#include "base/nibr.h"
#include "../algorithms/tractographyAlgorithm.h"
#include "../algorithms/ptt/algorithm_ptt_params.h"
#include "../seed/seed.h"
#include "../pathway/pathway.h"
#include "../io/tractogramField.h"

#define DEFAULT_IDLETIMELIMIT 120 // If not streamline is computed in the last 120 sec idleTimeLimitReached will be set to true

namespace NIBR {

namespace TRACKER 
{

// User parameters
extern Algorithm_Type  algorithm;
extern int             nThreads;
extern int             runTimeLimit;
extern int             idleTimeLimit;
extern Seed            seed;
extern bool            saveSeedIndex;
extern Pathway 		   pw;
extern Params_PTT      params_ptt;

// Derived parameters
extern Tractogram                                   tractogram;
extern std::vector<int>                             streamlineLength;
extern std::vector<int>                             seedIndex;
extern TractogramField                              seedIndexField;
extern std::chrono::steady_clock::time_point        initTime;
extern std::chrono::steady_clock::time_point        lastTime;
extern bool                                         runtimeLimitReached;  // time passed since the TRACKER was initialized or reset
extern bool                                         idletimeLimitReached; // time passed since the last successful streamline was computed and appended on the tractogram
extern std::size_t                                  currentCount;
extern bool                                         countIsReached;
extern int                                          ready_thread_id;
extern std::mutex                                   trackKeeper;

// Logger
struct Logger {

    std::atomic<std::size_t> log_success_REACHED_MAXLENGTH_LIMIT{0};
    std::atomic<std::size_t> log_success_REACHED_MINDATASUPPORT_LIMIT{0};
    std::atomic<std::size_t> log_success_SATISFIED_PATHWAY_RULES{0};

    std::atomic<std::size_t> log_discard_TOO_SHORT{0};
    std::atomic<std::size_t> log_discard_TOO_LONG{0};
    std::atomic<std::size_t> log_discard_DISCARD_ROI_REACHED{0};
    std::atomic<std::size_t> log_discard_REQUIRED_ROI_NOT_MET{0};
    std::atomic<std::size_t> log_discard_REQUIRED_ROI_ORDER_NOT_MET{0};
    std::atomic<std::size_t> log_discard_CANT_MEET_STOP_CONDITION{0};
    std::atomic<std::size_t> log_discard_ENDED_INSIDE_DISCARD_ROI{0};
    std::atomic<std::size_t> log_discard_REACHED_TIME_LIMIT{0};

    std::atomic<std::size_t> log_failed_UNKNOWN_REASON{0};
    std::atomic<std::size_t> log_failed_BY_THE_ALGORITHM_AT_INITIALIZATION{0};
    std::atomic<std::size_t> log_failed_BY_THE_ALGORITHM{0};

    std::atomic<std::size_t> log_unexpected_TRACKING_STATUS{0};
};

extern Logger trackerLogger;



Tractogram&         getTractogram();
Logger&             getLogger();
TractogramField&    getSeedIndexField();

int  runTime();
int  idleTime();
bool isWithinTimeLimits();

void reset();
void print();

}

}