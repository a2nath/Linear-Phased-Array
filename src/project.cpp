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
    auto args = argparse::parse<MyArgs>(argc, argv);

    cout << "arguments as follows:\n" << "-----------------------" << endl;
    args.print();

    /* setup the simulation runtime parameters */
    args.init();

    Logger logger(args.output_dir.value()  + "/output.txt");
    Simulator sim(args, logger);

    /* run the simulation and get the SNR table */
    sim.run();
    auto& perfmon = sim.get_perf();

    //plot(mobile_station_pos);

    for (auto& entrylist : perfmon.get_data())
    {
        for (auto& entry : entrylist.second)
        {
            cout << "SINR:" << entrylist.first
                << " " << entry << endl;
        }
    }

    close(logger);
    return 0;
}
