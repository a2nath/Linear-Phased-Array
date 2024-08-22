// project.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <queue>
#include "project.h"
//#include "matplot.h"

#ifdef GRAPHICS
#include <SFML/Graphics.hpp>
#endif

#define GRID_WIDTH 1000
#define GRID_HEIGHT 1000

using namespace std;

/* in case of errors, check this */

/* write to a file: loggerfile, string, debug-print */
void print(Logger& logger, const string& str, bool debug = true)
{
    logger.write(str);
    if (debug)
        cout << str << endl;
}

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
    sim.guistat(GRID_WIDTH, GRID_HEIGHT);


    //plot(mobile_station_pos);


    close(logger);
    return 0;
}
