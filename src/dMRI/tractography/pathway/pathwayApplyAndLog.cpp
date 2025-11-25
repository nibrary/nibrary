#include "pathway.h"

using namespace NIBR;


// Starts the interactive logging process in a separate thread.
void NIBR::Pathway::startLogger(int N)
{
    pathwayLog.reset(); 
    pathwayLog.numberOfStreamlinesToProcess = N;
    isLogging = true;
}

// Stops the interactive logging thread and prints a final summary.
void NIBR::Pathway::stopLogger()
{ 
    isLogging = false;

    if (isShowing.load()) { 
        hideLogger();
    }
}

void NIBR::Pathway::showLogger()
{
    isShowing = true;
    logDisplayer = std::thread([&]() {this->printLogger();});
}

void NIBR::Pathway::hideLogger()
{
    isShowing = false;
    if (logDisplayer.joinable()) {
        logDisplayer.join();
    }
}

void NIBR::Pathway::printLogger()
{

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

    auto runTime = [&]()->int {
        return int(std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - pathwayLog.initTime).count());
    };

    pathwayLog.initTime = std::chrono::steady_clock::now();

    auto printReport = [&]()->void {

        if (pathwayLog.numberOfStreamlinesProcessed.load() == 0) {
            pathwayLog.initTime = std::chrono::steady_clock::now();
            return;
        }

        std::cout << "\r\033[K" << std::flush;           // Clear the current line
        for (int n=0; n<(lineCount-1); n++)
            std::cout << "\033[A\r\033[K" << std::flush; // Clear previous lines

        lineCount = success = discard = fail = 0;

        std::cout << preamble << "Success report" << "\033[0m" << std::endl; lineCount++;
        dispCount(success,  pathwayLog.log_success_REACHED_MAXLENGTH_LIMIT.load(),              "Reached max length");
        dispCount(success,  pathwayLog.log_success_REACHED_MINDATASUPPORT_LIMIT.load(),         "Reached end of streamline");
        dispCount(success,  pathwayLog.log_success_SATISFIED_PATHWAY_RULES.load(),              "Satisfied pathway rules");

        std::cout << preamble << "Discard report" << "\033[0m" << std::endl; lineCount++;
        dispCount(discard,  pathwayLog.log_discard_DISCARDINGREASON_NOTSET.load(),              "Reason not set");
        dispCount(discard,  pathwayLog.log_discard_TOO_SHORT.load(),                            "Too short");
        dispCount(discard,  pathwayLog.log_discard_TOO_LONG.load(),                             "Too long");
        dispCount(discard,  pathwayLog.log_discard_DISCARD_ROI_REACHED.load(),                  "Reached discard region");
        dispCount(discard,  pathwayLog.log_discard_REQUIRED_ROI_NOT_MET.load(),                 "Required region not found");
        dispCount(discard,  pathwayLog.log_discard_REQUIRED_ROI_ORDER_NOT_MET.load(),           "Required order not satisfied");
        dispCount(discard,  pathwayLog.log_discard_CANT_MEET_STOP_CONDITION.load(),             "Can't meet stop condition");
        dispCount(discard,  pathwayLog.log_discard_ENDED_INSIDE_DISCARD_ROI.load(),             "Ended inside discard region");
        dispCount(discard,  pathwayLog.log_discard_REACHED_TIME_LIMIT.load(),                   "Reached time limit");
        dispCount(discard,  pathwayLog.log_discard_SEED_NOT_FOUND.load(),                       "Seed not found");
        dispCount(discard,  pathwayLog.log_discard_DISCARD_SEED.load(),                         "Seed was discarded");
        dispCount(discard,  pathwayLog.log_discard_IMPROPER_SEED.load(),                        "Improper seed");


        std::cout << std::endl << preamble;
        std::cout << "Success: "  << success              << "    ";
        std::cout << "Discard: "  << discard              << "    ";
        std::cout << "Total: "    << success+discard+fail << "    ";
        std::cout << "Progress: " << 100 * (success+discard+fail) / float(pathwayLog.numberOfStreamlinesToProcess)  << " %    " << std::fixed << std::setprecision(2);
        std::cout << "Duration: " << runTime()   << " sec";
        std::cout << "\033[0m"    << std::flush;
        lineCount=lineCount+2;
    };

    while (isLogging.load())
    {
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        printReport();
    }

    printReport();
    std::cout << std::endl;

    return;

}

// Returns false and empty streamline if discarded, true and the output streamline if kept
std::tuple<bool,NIBR::Streamline> NIBR::Pathway::apply(NIBR::Streamline& streamline)
{
    disp(MSG_DEBUG,"Applying filter on streamline");

    // disp(MSG_DEBUG,"streamline size before: %d", streamline.size());
    // for (auto& p : streamline) {
    //     disp(MSG_DEBUG, "%.2f,%.2f,%.2f", p[0],p[1],p[2]);
    // }

    Walker* walker = createWalker(&streamline);

    disp(MSG_DEBUG,"Created walker");

    if (hasSeed()) {
        seededProcess(walker);
        // for (auto& p : streamline) {
        //     disp(MSG_DEBUG, "%.2f,%.2f,%.2f", p[0],p[1],p[2]);
        // }
        if (walker->seedInserted) {
            streamline.erase(streamline.begin() + walker->seedInd);
        }
    } else {
        seedlessProcess(walker);
    }

    // disp(MSG_DEBUG,"streamline size after (possibly with added seed): %d", streamline.size());
    // for (auto& p : streamline) {
    //     disp(MSG_DEBUG, "%.2f,%.2f,%.2f", p[0],p[1],p[2]);
    // }

    disp(MSG_DEBUG,"Walked the pathway");

    // Update log
    if (isLogging.load()) {

        disp(MSG_DEBUG,"Updating log");
        pathwayLog.numberOfStreamlinesProcessed.fetch_add(1);

        if (walker->action == DISCARD) {

            switch (walker->discardingReason) {
                case DISCARDINGREASON_NOTSET:
                    disp(MSG_DEBUG, "Streamline: %d - DISCARDINGREASON_NOTSET", walker->ind);
                    pathwayLog.log_discard_DISCARDINGREASON_NOTSET.fetch_add(1);
                    break;
                case TOO_SHORT:
                    disp(MSG_DEBUG, "Streamline: %d - TOO_SHORT", walker->ind);
                    pathwayLog.log_discard_TOO_SHORT.fetch_add(1);
                    break;
                case TOO_LONG:
                    disp(MSG_DEBUG, "Streamline: %d - TOO_LONG", walker->ind);
                    pathwayLog.log_discard_TOO_LONG.fetch_add(1);
                    break;
                case DISCARD_REGION_REACHED:
                    disp(MSG_DEBUG, "Streamline: %d - DISCARD_REGION_REACHED", walker->ind);
                    pathwayLog.log_discard_DISCARD_ROI_REACHED.fetch_add(1);
                    break;
                case REQUIRED_ROI_NOT_MET:
                    disp(MSG_DEBUG, "Streamline: %d - REQUIRED_ROI_NOT_MET", walker->ind);
                    pathwayLog.log_discard_REQUIRED_ROI_NOT_MET.fetch_add(1);
                    break;
                case REQUIRED_ORDER_NOT_MET:
                    disp(MSG_DEBUG, "Streamline: %d - REQUIRED_ROI_ORDER_NOT_MET", walker->ind);
                    pathwayLog.log_discard_REQUIRED_ROI_ORDER_NOT_MET.fetch_add(1);
                    break;
                case CANT_MEET_STOP_CONDITION:
                    disp(MSG_DEBUG, "Streamline: %d - CANT_MEET_STOP_CONDITION", walker->ind);
                    pathwayLog.log_discard_CANT_MEET_STOP_CONDITION.fetch_add(1);
                    break;
                case ENDED_INSIDE_DISCARD_ROI:
                    disp(MSG_DEBUG, "Streamline: %d - ENDED_INSIDE_DISCARD_ROI", walker->ind);
                    pathwayLog.log_discard_ENDED_INSIDE_DISCARD_ROI.fetch_add(1);
                    break;
                case REACHED_TIME_LIMIT:
                    disp(MSG_DEBUG, "Streamline: %d - REACHED_TIME_LIMIT", walker->ind);
                    pathwayLog.log_discard_REACHED_TIME_LIMIT.fetch_add(1);
                    break;
                case SEED_NOT_FOUND:
                    disp(MSG_DEBUG, "Streamline: %d - SEED_NOT_FOUND", walker->ind);
                    pathwayLog.log_discard_SEED_NOT_FOUND.fetch_add(1);
                    break;
                case DISCARD_SEED:
                    disp(MSG_DEBUG, "Streamline: %d - DISCARD_SEED", walker->ind);
                    pathwayLog.log_discard_DISCARD_SEED.fetch_add(1);
                    break;
                case IMPROPER_SEED:
                    disp(MSG_DEBUG, "Streamline: %d - IMPROPER_SEED", walker->ind);
                    pathwayLog.log_discard_IMPROPER_SEED.fetch_add(1);
                    break;
                default:
                    break;
            }

            disp(MSG_DEBUG,"DISCARDED");

        }

        if (walker->action == KEEP) {

            if      ((walker->terminationReasonSideA==MAX_LENGTH_REACHED)      || (walker->terminationReasonSideB==MAX_LENGTH_REACHED)     )
                pathwayLog.log_success_REACHED_MAXLENGTH_LIMIT.fetch_add(1);
            else if ((walker->terminationReasonSideA==MIN_DATASUPPORT_REACHED) || (walker->terminationReasonSideB==MIN_DATASUPPORT_REACHED))
                pathwayLog.log_success_REACHED_MINDATASUPPORT_LIMIT.fetch_add(1);
            else
                pathwayLog.log_success_SATISFIED_PATHWAY_RULES.fetch_add(1);
            disp(MSG_DEBUG,"[begInd,endInd]: %.2f-%.2f",walker->begInd,walker->endInd);
            disp(MSG_DEBUG,"KEPT");
        } else {   
            
        }
        
        disp(MSG_DEBUG,"Updated log");

    }

    // Prepare output

    // Return the original streamline if discarded
    if (walker->action == DISCARD) {
        delete walker;
        return {false,Streamline()};
    }

    // We will keep this streamline

    // Reverse the begin and end indexes since walker preserves streamline direction
    if(walker->endInd < walker->begInd) std::swap(walker->begInd,walker->endInd);

    // Copy and return input streamline if no crop is needed
    if ((walker->begInd == 0) && (walker->endInd == (streamline.size()-1)) ) {
        delete walker;
        return {true,streamline}; // same as input
    }

    // Crop and return input streamline
    Streamline out;

    // Reverse the begin and end indexes
    if(walker->endInd < walker->begInd) std::swap(walker->begInd,walker->endInd);

    int   intBeg   = (int)walker->begInd;
    float fractBeg = walker->begInd-intBeg;

    int   intEnd   = (int)walker->endInd;
    float fractEnd = walker->endInd-intEnd;

    // disp(MSG_DEBUG,"wbeg - wend: %.6f - %.6f", walker->begInd, walker->endInd);
    // disp(MSG_DEBUG,"ibeg - iend: %d - %d", intBeg, intEnd);
    // disp(MSG_DEBUG,"fbeg - fend: %.6f - %.6f", fractBeg, fractEnd);

    // Handle first point
    if (fractBeg>0.0f) {
        NIBR::Point3D p;
        p[0] = streamline[intBeg][0] * (1 - fractBeg) + streamline[intBeg + 1][0] * (fractBeg);
        p[1] = streamline[intBeg][1] * (1 - fractBeg) + streamline[intBeg + 1][1] * (fractBeg);
        p[2] = streamline[intBeg][2] * (1 - fractBeg) + streamline[intBeg + 1][2] * (fractBeg);
        out.emplace_back(p);
        intBeg++;
    }

    // Handle middle points
    for (int j=intBeg; j<=intEnd; j++) {
        out.push_back({streamline[j][0],streamline[j][1],streamline[j][2]});
    }

    // Handle last point
    if (fractEnd>0.0f) {
        NIBR::Point3D p;
        p[0] = streamline[intEnd][0] * (1 - fractEnd) + streamline[intEnd + 1][0] * (fractEnd);
        p[1] = streamline[intEnd][1] * (1 - fractEnd) + streamline[intEnd + 1][1] * (fractEnd);
        p[2] = streamline[intEnd][2] * (1 - fractEnd) + streamline[intEnd + 1][2] * (fractEnd);
        out.emplace_back(p);
    }    

    delete walker;

    return {true,out};

}

// Returns kept streamlines, and discarded streamline if outputDiscard is true
std::tuple<NIBR::StreamlineBatch,std::vector<int>,NIBR::StreamlineBatch> NIBR::Pathway::apply(NIBR::StreamlineBatch& inpBatch, bool outputDiscarded)
{

    StreamlineBatch outKeep;
    StreamlineBatch outDiscard;

    outKeep.reserve(inpBatch.size());
    outDiscard.reserve(inpBatch.size());

    std::vector<int> outIndices;
    outIndices.reserve(inpBatch.size());

    std::mutex modifier;
    
    auto filterBatch = [&](const MT::TASK& task) {
        auto out = apply(inpBatch[task.no]);
        if (std::get<0>(out)) {
            std::lock_guard lock(modifier);
            outKeep.emplace_back(std::get<1>(out));
            outIndices.push_back(task.no);
        } else if (outputDiscarded) {
            std::lock_guard lock(modifier);
            outDiscard.emplace_back(inpBatch[task.no]);
        }
    };
    MT::MTRUN(inpBatch.size(), filterBatch);

    return {std::move(outKeep), std::move(outIndices), std::move(outDiscard)};

}
