#pragma once

#ifdef _WIN32
#include <windows.h>
#include <io.h>
#undef max
#else
#include <unistd.h>
#endif

#include <iostream>
#include <string>
#include <iomanip>
#include <thread>
#include <vector>
#include <memory>
#include <functional> 
#include <mutex>
#include <condition_variable>

#if defined(_WIN32)
#include <tlhelp32.h>
#elif defined(__APPLE__)
#include <libproc.h>
#else
#include <filesystem>
#endif

#include "base/nibr.h"
#include "math/randomThings.h"

namespace NIBR 
{

    namespace MT 
    {

        struct TASK {
            std::size_t   no;
            uint16_t threadId;
        };

        int&                                             MAXNUMBEROFTHREADS();
        int                                              GETNUMBEROFRUNNINGTHREADS();
        void                                             SETMAXNUMBEROFTHREADS(int n);
        std::atomic<std::size_t>&                        FINISHEDTASKCOUNT();
        std::atomic<std::size_t>&                        FINISHEDTASKCOUNTTOSTOP();
        uint16_t&                                        FINISHEDTHREADID();
        void                                             SET_DISP_RANGE(std::size_t start, std::size_t total);
        std::vector<std::unique_ptr<NIBR::RandomDoer>>&  RNDM();        

        // C++17 compatible barrier class (This already exists in C++20)
        // Allows a fixed number of threads to wait until all threads reach the barrier.
        class Barrier {
            public:
                Barrier(std::size_t _threshold) { 
                    threshold  = _threshold;
                    count      = _threshold;
                    generation = 0;
                }
    
                // Deleted copy constructor and assignment operator
                Barrier(const Barrier&) = delete;
                Barrier& operator=(const Barrier&) = delete;
    
                // Function for threads to call when they reach the barrier point.
                // Blocks until all threads have called arrive_and_wait for the current synchronization cycle (generation).
                void arrive_and_wait() {
                    std::unique_lock<std::mutex> lock(mute);
                    auto current_gen = generation;
                    count--; // Decrement count for the current generation
    
                    if (count == 0) {
                        // Last thread arrived for this generation
                        generation++;           // Move to the next generation
                        count = threshold;      // Reset count for the next use
                        cv.notify_all();        // Wake up all waiting threads
                    } else {
                        // Not the last thread, wait until the generation changes
                        cv.wait(lock, [this, current_gen] { return generation != current_gen; });
                    }
                    // All threads proceed past this point only after the last thread
                    // has arrived and notified the others for the current generation.
                }
    
            private:
                std::mutex mute;                // Mutex to protect internal state
                std::condition_variable cv;     // Condition variable for waiting
                std::size_t threshold;          // Total number of threads expected at the barrier
                std::size_t count{0};           // Number of threads still expected in the current generation
                std::size_t generation{0};      // Tracks the current barrier cycle/generation
        };


        void MTINIT();

        // Core functions
        void MTRUN(std::size_t range, int numberOfThreads, std::function<void(const NIBR::MT::TASK&, NIBR::MT::Barrier&)> f);
        void MTRUN(std::size_t range, int numberOfThreads, std::function<bool(const NIBR::MT::TASK&, NIBR::MT::Barrier&)> f, std::size_t stopLim);
        void MTRUN(std::size_t range, int numberOfThreads, std::string message, std::function<void(const NIBR::MT::TASK&, NIBR::MT::Barrier&)> f);
        void MTRUN(std::size_t range, int numberOfThreads, std::string message, std::function<bool(const NIBR::MT::TASK&, NIBR::MT::Barrier&)> f, std::size_t stopLim);

        std::vector<std::pair<int,int>> createTaskRange(int taskCount, int workerCount);


        // Overrides
        void MTRUN(std::size_t range, std::function<void(const NIBR::MT::TASK&, NIBR::MT::Barrier&)> f);
        void MTRUN(std::size_t range, std::function<bool(const NIBR::MT::TASK&, NIBR::MT::Barrier&)> f, std::size_t stopLim);
        
        void MTRUN(std::size_t range, std::string message, std::function<void(const NIBR::MT::TASK&, NIBR::MT::Barrier&)> f);
        void MTRUN(std::size_t range, std::string message, std::function<bool(const NIBR::MT::TASK&, NIBR::MT::Barrier&)> f, std::size_t stopLim);

        void MTRUN(std::size_t range, int numberOfThreads, std::function<void()> f);
        void MTRUN(std::size_t range, int numberOfThreads, std::function<void(const NIBR::MT::TASK&)> f);

        void MTRUN(std::size_t range, int numberOfThreads, std::function<bool()> f, std::size_t stopLim);
        void MTRUN(std::size_t range, int numberOfThreads, std::function<bool(const NIBR::MT::TASK&)> f, std::size_t stopLim);

        void MTRUN(std::size_t range, int numberOfThreads, std::string message, std::function<void()> f);
        void MTRUN(std::size_t range, int numberOfThreads, std::string message, std::function<void(const NIBR::MT::TASK&)> f);

        void MTRUN(std::size_t range, int numberOfThreads, std::string message, std::function<bool()> f, std::size_t stopLim);
        void MTRUN(std::size_t range, int numberOfThreads, std::string message, std::function<bool(const NIBR::MT::TASK&)> f, std::size_t stopLim);

        void MTRUN(std::size_t range, std::function<void()> f);
        void MTRUN(std::size_t range, std::function<void(const NIBR::MT::TASK&)> f);

        void MTRUN(std::size_t range, std::function<bool()> f, std::size_t stopLim);
        void MTRUN(std::size_t range, std::function<bool(const NIBR::MT::TASK&)> f, std::size_t stopLim);

        void MTRUN(std::size_t range, std::string message, std::function<void()> f);
        void MTRUN(std::size_t range, std::string message, std::function<void(const NIBR::MT::TASK&)> f);

        void MTRUN(std::size_t range, std::string message, std::function<bool()> f, std::size_t stopLim);
        void MTRUN(std::size_t range, std::string message, std::function<bool(const NIBR::MT::TASK&)> f, std::size_t stopLim);

    }

}