// project.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include "project.h"


#ifdef _WIN32
#include <windows.h>
#define getpid GetCurrentProcessId
#endif

using namespace std;

void close(Logger& logger)
{
    logger.close();
}

int main(int argc, char** argv, char** envp)
{
#ifdef _DEBUG
    std::cout << "PID: " << getpid() << std::endl;
#endif

    auto args = argparse::parse<MyArgs>(argc, argv);

    /* setup the simulation runtime parameters */
    args.init();

    cout << "arguments as follows:\n" << "-----------------------" << endl;
    args.print();

    Logger logger(args.get_input_filename(), args.output_dir);
    Simulator sim(args, logger);

    /* run the simulation and get the SNR table */


    sim.run();
    sim.print();
    sim.gui_run();
    sim.gui_print();

    close(logger);
    return 0;
}
