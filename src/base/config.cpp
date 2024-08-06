#include "config.h"
#include "multithreader.h"
#include "verbose.h"
#include <string>

#include <geogram/basic/common.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/logger.h>

namespace NIBR 
{
    std::string   sgntr   = "nibrary v0.1.0";
    VERBOSE_LEVEL verbose = VERBOSE_INFO;

    std::string&   SGNTR()   {return sgntr;}
    VERBOSE_LEVEL& VERBOSE() {return verbose;}

    std::string  							DATESTAMP;
    std::chrono::steady_clock::time_point 	INITTIME;

    int TICTIME;
    int TOCTIME;

}

using namespace NIBR;

void NIBR::INITIALIZE() 
{
    // Ensure initialization happens only once
    static std::once_flag config_init_flag;
    std::call_once(config_init_flag, []() {
        time_t startDate;
        time (&startDate);
        struct tm * timeinfo;
        timeinfo = localtime(&startDate);
        std::ostringstream tmp;
        tmp << std::put_time(timeinfo, "%d %b %Y, %H:%M:%S");
        DATESTAMP = tmp.str();

        // Initialize multithreader
        INITTIME = std::chrono::steady_clock::now();

        NIBR::MT::MTINIT();

        // Initialize geogram 
        disableTerminalOutput();
        GEO::initialize();

        disableTerminalOutput();
        GEO::Process::enable_multithreading(false);

        disableTerminalOutput();
        GEO::Process::set_max_threads(1);

        disableTerminalOutput();
        GEO::CmdLine::import_arg_group("standard");

        disableTerminalOutput();
        GEO::CmdLine::import_arg_group("algo");

        disableTerminalOutput();
        GEO::Logger::instance()->set_quiet(true);
        
        enableTerminalOutput();

    });

}

void NIBR::TERMINATE()
{
    // Terminate geogram
    // GEO::terminate();
}

int NIBR::RUNTIME() 
{
	return int(std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now()-INITTIME).count());
}

int NIBR::MSECRUNTIME() 
{
	return int(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now()-INITTIME).count());
}

void NIBR::TIC() 
{
	TICTIME = MSECRUNTIME();
}

void NIBR::TOC() 
{
	TOCTIME = MSECRUNTIME();
    disp(MSG_INFO,"Elapsed time: %.3f sec", float(TOCTIME-TICTIME)/1000);
}

void NIBR::TOC(MESSAGE msg) 
{
	TOCTIME = MSECRUNTIME();
    disp(msg,"Elapsed time: %.3f sec (%s, line: %d)", float(TOCTIME-TICTIME)/1000, msg.fileName.c_str(),msg.lineNo);
}
