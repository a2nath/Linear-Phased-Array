# README #

## Table of Contents
-   [Problem statement](#problem-statement)
-   [Parameters](#parameters)
-   [Initial Solution](#initial-solution-in-matlab)
-   [Optimized Solution](#optimized-solution)
-   [Building and Integrating](#building-and-integrating)

## Problem statement 
There are 3 "cellular on wheels (COWs) portable base stations, each of which can drive a linear scan antenna array of up to 8 elements. The layout consists of 12 non-VIP spectator zones and 3 VIP spectator zones. The goal is to provide at least 24 dB of SINR ("SNR") of signal quality to all clients involved in the simulation. See `wiki/original_problem.pdf` for more info.

More details on the linear scan antenna array as [follows.](https://en.wikipedia.org/wiki/Phased_array)
<p align="center">
<img src="images/280px-Phased_array_animation_with_arrow_10frames_371x400px_100ms.gif?raw=true" width="400" alt="How it works">

 In the project the scan angles varies between -90 and +90 degrees and each COW is placed at an orientation of 0 degrees (theta) or facing East, and theta-C angle is used to convey the scan angle as seen in the image.

</p>
<p align="center"><img src="images/Phasearray.gif?raw=true" width="400" alt="How it works">
</p>


## Parameters
- one frequency channel of 20 MHz wide, f-center = 1900 MHz
- each COW has a mast that can mount and linearly scan up to 8 antennas driven from a common downlink modulator signal. 
- each antenna has hx = 20cm and hy = 40 cm.
- overall width of the array cannot be more than 240 cm.
- each of the modulated downlink streams can be drive with power from -30 dBm to +30 dBm in 1 dB steps between timeslots

## Solution
### Initial solution in Matlab

<details open>

Initial proposal and source code along with optimizations have been made in `wiki/possible_solution.pdf`. Some of the visual effects of using linear scan antenna array can be found below: 

<p align="center">
<img src="images/antenna-profile0.png?raw=true" width="800" alt="Log scan #1">
</p>

<p align="center">
<img src="images/antenna-profile1.png?raw=true" width="800" alt="Log scan #2">
</p>

<p align="center">
<img src="images/antenna-profile2.png?raw=true" width="800" alt="Log scan #3">
</p>



<p align="center">
<img src="images/bs3_lin.png?raw=true" width="800" alt="Linear plot with focus on MS #5">
</p>
</details>

### Optimized Solution
TODO

## Building and Integrating

Default is release build with -O3 optimization

```
mkdir build
cd build
cmake ..
cmake --build .
```

### Debug build w/ debug symbols
Change `cmake ..` to `cmake -DCMAKE_BUILD_TYPE=Debug ..`



