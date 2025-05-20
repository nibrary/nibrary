#include "pathFilter.h"
#include "streamline_operators.h"
#include "dMRI/tractography/io/tractogramWriter.h"
#include <cstdint>
#include <tuple>
#include <chrono>
#include <atomic>
#include "base/vectorOperations.h"

using namespace NIBR;

std::tuple<std::vector<size_t>,std::vector<float>,std::vector<float>> NIBR::pathFilter(NIBR::TractogramReader* _tractogram, Pathway* pw, int numberOfThreads, int stopLim) {

    if ( (_tractogram->numberOfStreamlines<1) || (pw->isverified() == false) )
        return std::tuple<std::vector<size_t>,std::vector<float>,std::vector<float>>();

    NIBR::TractogramReader *tractogram = new NIBR::TractogramReader[numberOfThreads]();
    for (int t = 0; t < numberOfThreads; t++)
        tractogram[t].copyFrom(*_tractogram);

    std::vector<std::vector<int>>   no;
    std::vector<std::vector<float>> begInd;
    std::vector<std::vector<float>> endInd;  
    int N = tractogram[0].numberOfStreamlines;

    for (int i = 0; i < numberOfThreads; i++) { 
        no.push_back(std::vector<int>());
        begInd.push_back(std::vector<float>());
        endInd.push_back(std::vector<float>());
    }

    std::atomic<int> log_success_REACHED_MAXLENGTH_LIMIT(0);
    std::atomic<int> log_success_REACHED_MINDATASUPPORT_LIMIT(0);
    std::atomic<int> log_success_SATISFIED_PATHWAY_RULES(0);

    std::atomic<int> log_discard_DISCARDINGREASON_NOTSET(0);
    std::atomic<int> log_discard_TOO_SHORT(0);
    std::atomic<int> log_discard_TOO_LONG(0);
    std::atomic<int> log_discard_DISCARD_ROI_REACHED(0);
    std::atomic<int> log_discard_REQUIRED_ROI_NOT_MET(0);
    std::atomic<int> log_discard_REQUIRED_ROI_ORDER_NOT_MET(0);
    std::atomic<int> log_discard_CANT_MEET_STOP_CONDITION(0);
    std::atomic<int> log_discard_ENDED_INSIDE_DISCARD_ROI(0);
    std::atomic<int> log_discard_REACHED_TIME_LIMIT(0);
    std::atomic<int> log_discard_SEED_NOT_FOUND(0);
    std::atomic<int> log_discard_DISCARD_SEED(0);
    std::atomic<int> log_discard_IMPROPER_SEED(0);


    bool filteringFinished = false;

    auto applyPathwayRule = [&](const NIBR::MT::TASK& task) -> bool {

        std::vector<Point> streamline = tractogram[task.threadId].readStreamlinePoints(task.no);

        disp(MSG_DEBUG,"Streamline %d, len: %d", task.no, streamline.size());

        Walker *walker = pw->createWalker(&streamline,task.no);

        disp(MSG_DEBUG,"Created");

        if (pw->hasSeed()) {
            pw->seededProcess(walker);
        } else {
            pw->seedlessProcess(walker);
        }

        disp(MSG_DEBUG,"Walked");

        if (walker->action == DISCARD) {

			switch (walker->discardingReason) {
                case DISCARDINGREASON_NOTSET:
                    disp(MSG_DEBUG, "Streamline: %d - DISCARDINGREASON_NOTSET", walker->ind);
                    log_discard_DISCARDINGREASON_NOTSET.fetch_add(1);
                    break;
                case TOO_SHORT:
                    disp(MSG_DEBUG, "Streamline: %d - TOO_SHORT", walker->ind);
                    log_discard_TOO_SHORT.fetch_add(1);
                    break;
                case TOO_LONG:
                    disp(MSG_DEBUG, "Streamline: %d - TOO_LONG", walker->ind);
                    log_discard_TOO_LONG.fetch_add(1);
                    break;
                case DISCARD_REGION_REACHED:
                    disp(MSG_DEBUG, "Streamline: %d - DISCARD_REGION_REACHED", walker->ind);
                    log_discard_DISCARD_ROI_REACHED.fetch_add(1);
                    break;
                case REQUIRED_ROI_NOT_MET:
                    disp(MSG_DEBUG, "Streamline: %d - REQUIRED_ROI_NOT_MET", walker->ind);
                    log_discard_REQUIRED_ROI_NOT_MET.fetch_add(1);
                    break;
                case REQUIRED_ORDER_NOT_MET:
                    disp(MSG_DEBUG, "Streamline: %d - REQUIRED_ROI_ORDER_NOT_MET", walker->ind);
                    log_discard_REQUIRED_ROI_ORDER_NOT_MET.fetch_add(1);
                    break;
                case CANT_MEET_STOP_CONDITION:
                    disp(MSG_DEBUG, "Streamline: %d - CANT_MEET_STOP_CONDITION", walker->ind);
                    log_discard_CANT_MEET_STOP_CONDITION.fetch_add(1);
                    break;
                case ENDED_INSIDE_DISCARD_ROI:
                    disp(MSG_DEBUG, "Streamline: %d - ENDED_INSIDE_DISCARD_ROI", walker->ind);
                    log_discard_ENDED_INSIDE_DISCARD_ROI.fetch_add(1);
                    break;
                case REACHED_TIME_LIMIT:
                    disp(MSG_DEBUG, "Streamline: %d - REACHED_TIME_LIMIT", walker->ind);
                    log_discard_REACHED_TIME_LIMIT.fetch_add(1);
                    break;
                case SEED_NOT_FOUND:
                    disp(MSG_DEBUG, "Streamline: %d - SEED_NOT_FOUND", walker->ind);
                    log_discard_SEED_NOT_FOUND.fetch_add(1);
                    break;
                case DISCARD_SEED:
                    disp(MSG_DEBUG, "Streamline: %d - DISCARD_SEED", walker->ind);
                    log_discard_DISCARD_SEED.fetch_add(1);
                    break;
                case IMPROPER_SEED:
                    disp(MSG_DEBUG, "Streamline: %d - IMPROPER_SEED", walker->ind);
                    log_discard_IMPROPER_SEED.fetch_add(1);
                    break;
                default:
                    break;
            }

		}

        disp(MSG_DEBUG,"Updated");

        if (walker->action == KEEP) {

            if      ((walker->terminationReasonSideA==MAX_LENGTH_REACHED)      || (walker->terminationReasonSideB==MAX_LENGTH_REACHED)     )
				log_success_REACHED_MAXLENGTH_LIMIT.fetch_add(1);
			else if ((walker->terminationReasonSideA==MIN_DATASUPPORT_REACHED) || (walker->terminationReasonSideB==MIN_DATASUPPORT_REACHED))
				log_success_REACHED_MINDATASUPPORT_LIMIT.fetch_add(1);
			else
				log_success_SATISFIED_PATHWAY_RULES.fetch_add(1);


            no[task.threadId].push_back(task.no);
            begInd[task.threadId].push_back(walker->begInd);
            endInd[task.threadId].push_back(walker->endInd);
            disp(MSG_DEBUG,"pathFilter::[begInd,endInd]: %.2f-%.2f",walker->begInd,walker->endInd);
            delete walker;
            disp(MSG_DEBUG,"KEPT");
            return true;
        } else {   
            delete walker;
            disp(MSG_DEBUG,"DISCARDED");
            return false;
        }   

    };

    std::thread coreThread([&]() {
        if (stopLim==-1) NIBR::MT::MTRUN(N, applyPathwayRule);
        if (stopLim > 0) NIBR::MT::MTRUN(N, applyPathwayRule, stopLim);
        filteringFinished = true;
    });


    auto initTime = std::chrono::steady_clock::now();
    auto runTime = [&]()->int {
        return int(std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - initTime).count());
    };    


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
                reason = "   " + reason + std::string(spaces,' ') + ":";
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
            dispCount(success,  log_success_REACHED_MAXLENGTH_LIMIT.load(),              "Reached max length");
            dispCount(success,  log_success_REACHED_MINDATASUPPORT_LIMIT.load(),         "Reached end of streamline");
            dispCount(success,  log_success_SATISFIED_PATHWAY_RULES.load(),              "Satisfied pathway rules");

            std::cout << preamble << "Discard report" << "\033[0m" << std::endl; lineCount++;
            dispCount(discard,  log_discard_DISCARDINGREASON_NOTSET.load(),              "Reason not set");
            dispCount(discard,  log_discard_TOO_SHORT.load(),                            "Too short");
            dispCount(discard,  log_discard_TOO_LONG.load(),                             "Too long");
            dispCount(discard,  log_discard_DISCARD_ROI_REACHED.load(),                  "Reached discard region");
            dispCount(discard,  log_discard_REQUIRED_ROI_NOT_MET.load(),                 "Required region not found");
            dispCount(discard,  log_discard_REQUIRED_ROI_ORDER_NOT_MET.load(),           "Required order not satisfied");
            dispCount(discard,  log_discard_CANT_MEET_STOP_CONDITION.load(),             "Can't meet stop condition");
            dispCount(discard,  log_discard_ENDED_INSIDE_DISCARD_ROI.load(),             "Ended inside discard region");
            dispCount(discard,  log_discard_REACHED_TIME_LIMIT.load(),                   "Reached time limit");
            dispCount(discard,  log_discard_SEED_NOT_FOUND.load(),                       "Seed not found");
            dispCount(discard,  log_discard_DISCARD_SEED.load(),                         "Seed was discarded");
            dispCount(discard,  log_discard_IMPROPER_SEED.load(),                        "Improper seed");


            std::cout << std::endl << preamble;
            std::cout << "Success: "  << success              << "    ";
            std::cout << "Discard: "  << discard              << "    ";
            std::cout << "Total: "    << success+discard+fail << "    ";
            std::cout << "Progress: " << 100 * (success+discard+fail) / float(tractogram[0].numberOfStreamlines)  << " %    " << std::fixed << std::setprecision(2);
            std::cout << "Duration: " << runTime()   << " sec";
            std::cout << "\033[0m"    << std::flush;
            lineCount=lineCount+2;

        };


        while (!filteringFinished)
        {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
            printReport();
        }

        std::cout << std::endl;

    }

    coreThread.join();


    for (int t = 0; t < numberOfThreads; t++) 
        tractogram[t].destroyCopy();
    delete[] tractogram;  
    
    std::vector<size_t> idx; 
    std::vector<float>  begIdx; 
    std::vector<float>  endIdx;  

    for (int i = 0; i < numberOfThreads; i++) {  
        idx.insert(idx.end(), no[i].begin(), no[i].end());
        begIdx.insert(begIdx.end(), begInd[i].begin(), begInd[i].end());
        endIdx.insert(endIdx.end(), endInd[i].begin(), endInd[i].end());
    }

    if ((stopLim!=-1) && (int(idx.size()) >= stopLim)) {
        idx.resize(stopLim);
        begIdx.resize(stopLim);
        endIdx.resize(stopLim);
    }

    return std::make_tuple(idx, begIdx, endIdx);

}