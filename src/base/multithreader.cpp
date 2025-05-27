#include "multithreader.h"
#include "config.h"
#include "verbose.h"
#include <cstddef>
#include <atomic>

#include <geogram/basic/common.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/logger.h>

namespace NIBR 
{
    namespace MT 
    {
        std::mutex                                  prox_mx;
        int                                         maxNumberOfThreads;
        std::atomic<std::size_t>                    finishedTaskCount;
        std::atomic<std::size_t>                    finishedTaskCountToStop;
        uint16_t                                    finishedThreadId;
        std::vector<std::unique_ptr<RandomDoer>>    rndm;
        bool                                        dispRange_ext;
        std::size_t                                 dispRange_start;
        std::size_t                                 dispRange_total;
    }
}

using namespace NIBR;
using namespace NIBR::MT;

std::mutex&                                 NIBR::MT::PROC_MX()                     {return NIBR::MT::prox_mx;}
int&                                        NIBR::MT::MAXNUMBEROFTHREADS()          {return NIBR::MT::maxNumberOfThreads;}
std::atomic<std::size_t>&                   NIBR::MT::FINISHEDTASKCOUNT()           {return NIBR::MT::finishedTaskCount;}
std::atomic<std::size_t>&                   NIBR::MT::FINISHEDTASKCOUNTTOSTOP()     {return NIBR::MT::finishedTaskCountToStop;}
uint16_t&                                   NIBR::MT::FINISHEDTHREADID()            {return NIBR::MT::finishedThreadId;}
std::vector<std::unique_ptr<RandomDoer>>&   NIBR::MT::RNDM()                        {return NIBR::MT::rndm;}
void                                        NIBR::MT::SET_DISP_RANGE(std::size_t start, std::size_t total) {NIBR::MT::dispRange_start = start; NIBR::MT::dispRange_total = total; NIBR::MT::dispRange_ext = true;}


void NIBR::MT::SETMAXNUMBEROFTHREADS(int n)  {

    NIBR::MT::maxNumberOfThreads = ((n > 0) && (n <= NIBR::MT::maxNumberOfThreads)) ? n : NIBR::MT::maxNumberOfThreads;

    // Enabling geogram multithreading
    disableTerminalOutput();
    GEO::Process::enable_multithreading(NIBR::MT::maxNumberOfThreads > 1);
    GEO::Process::set_max_threads(NIBR::MT::maxNumberOfThreads); 
    enableTerminalOutput();

}

void NIBR::MT::MTINIT() 
{
    // Ensure initialization happens only once
    static std::once_flag MT_init_flag;
    std::call_once(MT_init_flag, []() {

        #ifdef _WIN32
            SYSTEM_INFO sysinfo;
            GetSystemInfo(&sysinfo);
            NIBR::MT::maxNumberOfThreads = sysinfo.dwNumberOfProcessors;
        #elif defined(__APPLE__)
            NIBR::MT::maxNumberOfThreads = sysconf(_SC_NPROCESSORS_ONLN);
        #else
            FILE* pipe = popen("nproc", "r");
            if (!pipe) {
                NIBR::disp(MSG_WARN,"Cannot determine number of processing units (nproc). 1 thread will be used.");
                NIBR::MT::maxNumberOfThreads = 1; // default value
            };
            char buffer[128];
            std::string result = "";
            while (!feof(pipe)) {
                if (fgets(buffer, 128, pipe) != nullptr)
                    result += buffer;
            }
            pclose(pipe);
            std::stringstream ss(result);
            ss >> NIBR::MT::maxNumberOfThreads;
            if (ss.fail()) {
                NIBR::disp(MSG_WARN,"Error parsing number of processing units (nproc). 1 thread will be used.");
                NIBR::MT::maxNumberOfThreads = 1; // default value
            }
        #endif

        if (NIBR::MT::maxNumberOfThreads < 1) {
            NIBR::MT::maxNumberOfThreads = 1; // default value
        }

        NIBR::MT::rndm.clear();
        for (int i=0; i<NIBR::MT::maxNumberOfThreads; i++)
            NIBR::MT::rndm.push_back(std::make_unique<RandomDoer>());

    });
}


void NIBR::MT::MTRUN(std::size_t range, int numberOfThreads, std::function<void(const TASK&, Barrier&)> f) 
{

    if (range == 0) return;
    
    if (numberOfThreads == 0) numberOfThreads = NIBR::MT::maxNumberOfThreads;
    numberOfThreads = std::min(numberOfThreads, static_cast<int>(range));

    NIBR::MT::finishedTaskCount = 0;

    NIBR::MT::Barrier taskBarrier(numberOfThreads);

    auto taskRunner = [&](uint16_t threadNo) {

        while (true) {

            std::size_t taskIndex = finishedTaskCount.fetch_add(1); // assigns the previous value before adding 1

            if (taskIndex >= range) break;  // No more tasks left

            try {
                TASK task = {taskIndex, threadNo};
                f(task,taskBarrier);
            } catch (const std::exception& e) {
                disp(MSG_FATAL, "Failed executing task %d", taskIndex);
            }

        }

    };

    std::vector<std::thread> threadPool;
    threadPool.reserve(numberOfThreads);

    for (int i = 0; i < numberOfThreads; i++) {
        threadPool.emplace_back(taskRunner, i);
    }

    for (auto& t : threadPool) {
        if (t.joinable()) {
            t.join();
        }
    }

    return;

}



void NIBR::MT::MTRUN(std::size_t range, int numberOfThreads, std::function<bool(const TASK&, Barrier&)> f, std::size_t stopLim) 
{

    if (range == 0) return;
    
    if (numberOfThreads == 0) numberOfThreads = NIBR::MT::maxNumberOfThreads;
    numberOfThreads = std::min(numberOfThreads, static_cast<int>(range));

    NIBR::MT::finishedTaskCount       = 0;
    NIBR::MT::finishedTaskCountToStop = 0;

    NIBR::MT::Barrier taskBarrier(numberOfThreads);

    auto taskRunner = [&](uint16_t threadNo) {

        while (true) {

            std::size_t taskIndex = finishedTaskCount.fetch_add(1); // assigns the previous value before adding 1

            if ( (finishedTaskCount >= range) || (finishedTaskCountToStop >= stopLim) ) break;  // No more tasks left

            try {
                TASK task = {taskIndex, threadNo};
                if (f(task,taskBarrier)==true) {finishedTaskCountToStop.fetch_add(1);}
            } catch (const std::exception& e) {
                disp(MSG_FATAL, "Failed executing task %d", taskIndex);
            }

        }

    };

    std::vector<std::thread> threadPool;
    threadPool.reserve(numberOfThreads);

    for (int i = 0; i < numberOfThreads; i++) {
        threadPool.emplace_back(taskRunner, i);
    }

    for (auto& t : threadPool) {
        if (t.joinable()) {
            t.join();
        }
    }

    return;

}


void NIBR::MT::MTRUN(std::size_t range, int numberOfThreads, std::string message, std::function<void(const TASK&, Barrier&)> f)
{

    std::thread coreThread([&]() {
        MTRUN(range, numberOfThreads, f);
    });

    // Hide message below VERBOSE_INFO
    if (NIBR::VERBOSE()>=VERBOSE_INFO) {

        std::size_t totCount = (dispRange_ext) ? dispRange_total : range;
        std::size_t begCount = (dispRange_ext) ? dispRange_start : 0.0f;

        float progressScaler = 100.0f / float(totCount);

        std::string preamble = "\033[1;32mNIBRARY::INFO: \033[0;32m";

        // Display initial message and progress
        if (begCount > 0) {
            std::cout << "\r\033[K" << std::flush; // Clear the current line
        }
        std::cout << preamble << message << ": " << std::fixed << std::setprecision(2) << float(begCount) * progressScaler << "%" << "\033[0m" << std::flush;

        finishedTaskCount = 0;
        while (range>finishedTaskCount)
        {
            float curCount = finishedTaskCount + begCount;
            std::cout << "\r\033[K" << std::flush; // Clear the current line
            std::cout << preamble << message << ": " << std::fixed << std::setprecision(2) << curCount * progressScaler << "%" << "\033[0m" << std::flush;
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }

        if ((begCount + range) == totCount) {
            std::cout << "\r\033[K" << preamble << message << ": 100%" << std::endl;
        }

    }

    coreThread.join();

    dispRange_ext = false;
    
    return;

}


void NIBR::MT::MTRUN(std::size_t range, int numberOfThreads, std::string message, std::function<bool(const TASK&, Barrier&)> f, std::size_t stopLim)
{

    std::thread coreThread([&]() {
        MTRUN(range, numberOfThreads, f, stopLim);
    });

    // Hide message below VERBOSE_INFO
    if (NIBR::VERBOSE()>=VERBOSE_INFO) {

        std::string preamble = "\033[1;32mNIBRARY::INFO: \033[0;32m";

        // Display initial message and progress
        std::cout << preamble << message << " (success) : 0%" << "\033[0m" << std::endl;
        std::cout << preamble << message << " (total)   : 0%" << "\033[0m" << std::flush;


        finishedTaskCount           = 0;
        finishedTaskCountToStop     = 0;
        float localTaskCount        = 0;
        float localTaskCountToStop  = 0;
        float progressScaler        = 100.0f/float(range);
        float progressScalerToStop  = 100.0f/float(stopLim);


        while ( (range>localTaskCount) && (stopLim>localTaskCountToStop))
        {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));

            std::cout << "\r\033[K" << std::flush;       // Clear the current line
            std::cout << "\033[A\r\033[K" << std::flush; // Clear the line above the previous line

            std::cout << preamble << message << " (success) : " << std::fixed << std::setprecision(2) << progressScalerToStop*localTaskCountToStop << "%" << "\033[0m" << std::endl;
            std::cout << preamble << message << " (total)   : " << std::fixed << std::setprecision(2) << progressScaler*localTaskCount             << "%" << "\033[0m" << std::flush;

            localTaskCount       = finishedTaskCount;
            localTaskCountToStop = finishedTaskCountToStop;

        }

        localTaskCount       = finishedTaskCount;
        localTaskCountToStop = finishedTaskCountToStop;

        if (localTaskCount>range)         localTaskCount       = range;
        if (localTaskCountToStop>stopLim) localTaskCountToStop = stopLim;

        std::cout << "\r\033[K" << std::flush;       // Clear the current line
        std::cout << "\033[A\r\033[K" << std::flush; // Clear the line above the previous line

        std::cout << preamble << message << " (success) : " << std::fixed << std::setprecision(2) << progressScalerToStop*localTaskCountToStop << "%" << "\033[0m" << std::endl;
        std::cout << preamble << message << " (total)   : " << std::fixed << std::setprecision(2) << progressScaler*localTaskCount             << "%" << "\033[0m" << std::flush;

    }

    std::cout << std::endl;

    coreThread.join();
    
    return;

}

std::vector<std::pair<int,int>> NIBR::MT::createTaskRange(int taskCount, int workerCount) 
{

    std::vector<std::pair<int,int>> taskRange;

    int baseTask  = taskCount / workerCount;
    int remainder = taskCount % workerCount;

    int begin = 0;

    for (int i = 0; i < workerCount; i++) {

        // Calculate the ending point based on the baseTask and whether
        // this worker should handle an extra task from the remainder
        
        int end = begin + baseTask - 1; 
        if (remainder > 0) {
            end++;
            remainder--;
        }

        taskRange.emplace_back(std::make_pair(begin, end));

        // Set the starting point for the next iteration
        begin = end + 1;
    }
    
    return taskRange;

}





// Overrides
void NIBR::MT::MTRUN(std::size_t range, std::function<void(const TASK&, Barrier&)> f)
{
    auto wrapper_f = [&](const TASK& mttask, Barrier& barrier) {
        return f(mttask,barrier);
    };
    NIBR::MT::MTRUN(range, maxNumberOfThreads, wrapper_f);
}


void NIBR::MT::MTRUN(std::size_t range, std::function<bool(const TASK&, Barrier&)> f, std::size_t stopLim)
{
    auto wrapper_f = [&](const TASK& mttask, Barrier& barrier) {
        return f(mttask,barrier);
    };
    NIBR::MT::MTRUN(range, maxNumberOfThreads, wrapper_f, stopLim);
}


void NIBR::MT::MTRUN(std::size_t range, std::string message, std::function<void(const TASK&, Barrier&)> f)
{
    auto wrapper_f = [&](const TASK& mttask, Barrier& barrier) {
        return f(mttask,barrier);
    };
    NIBR::MT::MTRUN(range, maxNumberOfThreads, message, wrapper_f);
}


void NIBR::MT::MTRUN(std::size_t range, std::string message, std::function<bool(const TASK&, Barrier&)> f, std::size_t stopLim)
{
    auto wrapper_f = [&](const TASK& mttask, Barrier& barrier) {
        return f(mttask,barrier);
    };
    NIBR::MT::MTRUN(range, maxNumberOfThreads, message, wrapper_f, stopLim);
}


void NIBR::MT::MTRUN(std::size_t range, int numberOfThreads, std::function<void()> f)
{
    auto wrapper_f = [&](const TASK&, Barrier&) {
        return f();
    };
    NIBR::MT::MTRUN(range, numberOfThreads, wrapper_f);
}

void NIBR::MT::MTRUN(std::size_t range, int numberOfThreads, std::function<void(const NIBR::MT::TASK&)> f)
{
    auto wrapper_f = [&](const TASK& mttask, Barrier&) {
        return f(mttask);
    };
    NIBR::MT::MTRUN(range, numberOfThreads, wrapper_f);
}

void NIBR::MT::MTRUN(std::size_t range, int numberOfThreads, std::function<bool()> f, std::size_t stopLim)
{
    auto wrapper_f = [&](const TASK&, Barrier&) {
        return f();
    };
    NIBR::MT::MTRUN(range, numberOfThreads, wrapper_f, stopLim);
}


void NIBR::MT::MTRUN(std::size_t range, int numberOfThreads, std::function<bool(const NIBR::MT::TASK&)> f, std::size_t stopLim)
{
    auto wrapper_f = [&](const TASK& mttask, Barrier&) {
        return f(mttask);
    };
    NIBR::MT::MTRUN(range, numberOfThreads, wrapper_f, stopLim);
}


void NIBR::MT::MTRUN(std::size_t range, int numberOfThreads, std::string message, std::function<void()> f)
{
    auto wrapper_f = [&](const TASK&, Barrier&) {
        return f();
    };
    NIBR::MT::MTRUN(range, numberOfThreads, message, wrapper_f);
}


void NIBR::MT::MTRUN(std::size_t range, int numberOfThreads, std::string message, std::function<void(const NIBR::MT::TASK&)> f)
{
    auto wrapper_f = [&](const TASK& mttask, Barrier&) {
        return f(mttask);
    };
    NIBR::MT::MTRUN(range, numberOfThreads, message, wrapper_f);
}



void NIBR::MT::MTRUN(std::size_t range, int numberOfThreads, std::string message, std::function<bool()> f, std::size_t stopLim)
{
    auto wrapper_f = [&](const TASK&, Barrier&) {
        return f();
    };
    NIBR::MT::MTRUN(range, numberOfThreads, message, wrapper_f, stopLim);
}


void NIBR::MT::MTRUN(std::size_t range, int numberOfThreads, std::string message, std::function<bool(const NIBR::MT::TASK&)> f, std::size_t stopLim)
{
    auto wrapper_f = [&](const TASK& mttask, Barrier&) {
        return f(mttask);
    };
    NIBR::MT::MTRUN(range, numberOfThreads, message, wrapper_f, stopLim);
}


void NIBR::MT::MTRUN(std::size_t range, std::function<void()> f)
{
    auto wrapper_f = [&](const TASK&, Barrier&) {
        return f();
    };
    NIBR::MT::MTRUN(range, wrapper_f);
}



void NIBR::MT::MTRUN(std::size_t range, std::function<void(const NIBR::MT::TASK&)> f)
{
    auto wrapper_f = [&](const TASK& mttask, Barrier&) {
        return f(mttask);
    };
    NIBR::MT::MTRUN(range, wrapper_f);
}



void NIBR::MT::MTRUN(std::size_t range, std::function<bool()> f, std::size_t stopLim)
{
    auto wrapper_f = [&](const TASK&, Barrier&) {
        return f();
    };
    NIBR::MT::MTRUN(range, wrapper_f, stopLim);
}


void NIBR::MT::MTRUN(std::size_t range, std::function<bool(const NIBR::MT::TASK&)> f, std::size_t stopLim)
{
    auto wrapper_f = [&](const TASK& mttask, Barrier&) {
        return f(mttask);
    };
    NIBR::MT::MTRUN(range, wrapper_f, stopLim);
}

void NIBR::MT::MTRUN(std::size_t range, std::string message, std::function<void()> f)
{
    auto wrapper_f = [&](const TASK&, Barrier&) {
        return f();
    };
    NIBR::MT::MTRUN(range, message, wrapper_f);
}


void NIBR::MT::MTRUN(std::size_t range, std::string message, std::function<void(const NIBR::MT::TASK&)> f)
{
    auto wrapper_f = [&](const TASK& mttask, Barrier&) {
        return f(mttask);
    };
    NIBR::MT::MTRUN(range, message, wrapper_f);
}

void NIBR::MT::MTRUN(std::size_t range, std::string message, std::function<bool()> f, std::size_t stopLim)
{
    auto wrapper_f = [&](const TASK&, Barrier&) {
        return f();
    };
    NIBR::MT::MTRUN(range, message, wrapper_f, stopLim);
}


void NIBR::MT::MTRUN(std::size_t range, std::string message, std::function<bool(const NIBR::MT::TASK&)> f, std::size_t stopLim)
{
    auto wrapper_f = [&](const TASK& mttask, Barrier&) {
        return f(mttask);
    };
    NIBR::MT::MTRUN(range, message, wrapper_f, stopLim);
}