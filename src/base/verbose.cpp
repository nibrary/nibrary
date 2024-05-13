#include "config.h"
#include <stdlib.h>
#include <string>

using namespace NIBR;

int NIBR::disp(MESSAGE msg, const char *format, ...) 
{

    if (int(NIBR::VERBOSE())<=int(msg.type))
        return 0;

    va_list args;
    va_start(args, format);

    char buffer[4096];
    int rc = vsnprintf(buffer, sizeof(buffer), format, args);

    va_end(args);
    if (rc)
    {

        switch (msg.type) {

            case msg_fatal:
            {
                std::cout << "\033[1;31mNIBRARY::FATAL: \033[0;31m" << buffer << " [" << msg.fileName <<" line:"<< msg.lineNo<< "]" << "\033[0m\n" << std::flush;
                break;
            }

            case msg_error:
            {
                std::cout << "\033[1;31mNIBRARY::ERROR: \033[0;31m" << buffer << " [" << msg.fileName <<" line:"<< msg.lineNo<< "]" << "\033[0m\n" << std::flush;
                break;
            }

            case msg_warn:
            {
                std::cout << "\033[1;38;5;208mNIBRARY::WARNING: \033[0;38;5;208m" << buffer <<"\033[0m\n" << std::flush;
                break;
            }

            case msg_info:
            {
                std::cout << "\033[1;32mNIBRARY::INFO: \033[0;32m" << buffer <<"\033[0m\n" << std::flush;
                break;
            }

            case msg_detail:
            {
                std::cout << "\033[1;35mNIBRARY::DETAIL: \033[0;35m" << buffer <<"\033[0m\n" << std::flush;
                break;
            }

            case msg_debug:
            {
                std::cout << "\033[1;93mNIBRARY::DEBUG: \033[0;93m" << buffer << " [" << msg.fileName <<" line:"<< msg.lineNo<< "]" << "\033[0m\n" << std::flush;
                break;
            }

            default: break;
        }
            
    }
    return rc;
}

void NIBR::wait(const char *format, ...) {
    va_list args;
    va_start(args, format);

    char buffer[4096];
    vsnprintf(buffer, sizeof(buffer), format, args);

    std::cout << "\033[1;32mNIBRARY::WAIT: \033[0;32m" << buffer <<"(Press Enter to continue...) \033[0m\n" << std::flush;

    std::cin.get();
}


// Use a saved stdout_fd
int saved_stdout_fd = -1;

void NIBR::disableTerminalOutput() 
{

    #ifdef BUILD_FOR_WINDOWS
    // freopen("NUL", "w", stdout);
    #else
    if (saved_stdout_fd == -1) {
        saved_stdout_fd = dup(1);
    }
    freopen("/dev/null", "w", stdout);
    #endif
}


void NIBR::enableTerminalOutput() 
{
    #ifdef BUILD_FOR_WINDOWS
    // We consider the simple case for the console.
    // There might be need to use Windows API functions for more complex cases.
    // freopen("CON", "w", stdout);
    #else
    // freopen("/dev/tty", "w", stdout);
    if (saved_stdout_fd != -1) {
        dup2(saved_stdout_fd,1);
        close(saved_stdout_fd);
        saved_stdout_fd = -1;
    }
    #endif
}


