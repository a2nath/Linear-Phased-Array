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

/* define externs */
std::mutex graphics::queue_mutex;
std::mutex graphics::finished_mutex;
std::condition_variable graphics::consig;
Dimensions<unsigned> graphics::render_space;

cached::Cache<double, double> cached::cache_sin;
cached::Cache<double, double> cached::cache_dbm2w;
cached::Cache<double, double> cached::cache_db2lin;
cached::Cache<double, double> cached::cache_pow;


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

    if (!args.plot)
        sim.gui_run();

    sim.gui_print();

    close(logger);
    return 0;
}
