# README #

## Table of Contents
-   [Problem statement](#problem-statement)
-   [Parameters](#parameters)
-   [Matlab Solution](#matlab-solution)
-   [CPP Solution and Simulation Output](#cpp-solution-and-simulation-output)
-   [Building and Integrating](#building-and-integrating)

## 👥 Problem statement 
There are **3 "Cellular on Wheels (COWs)"** portable base stations, each equipped to drive a **linear scan antenna array** with up to **8 elements**. The layout consists of **12 non-VIP** spectator zones and **3 VIP** spectator zones. The goal 🎯 is to ensure all clients in the simulation receive at least **24 dB of SINR** ("SNR") signal quality or at least **3 VIP** zones are above the cut off if there exists no solution to the specific scenario. 📶

More details on the linear scan antenna array as [follows.](https://en.wikipedia.org/wiki/Phased_array)
<p align="center">
<img src="images/280px-Phased_array_animation_with_arrow_10frames_371x400px_100ms.gif?raw=true" width="400" alt="How it works">

In this project, the **scan angles** vary between **-90°** and **+90°**, with each **COW** placed at an orientation of **0° (theta)** (facing East). The **theta-C angle** is used to represent the scan angle, as shown in the image. 

</p>
<p align="center"><img src="images/Phasearray.gif?raw=true" width="400" alt="How it works">
</p>


## 🔧 Parameters
**Frequency Channel**: 20 MHz wide, **f-center = 1900 MHz**.
- Each **COW** can mount and **linearly scan** up to **8 antennas** driven by a common downlink modulator signal. 🛰️
- **Antenna Dimensions**: Each antenna has **hx = 20cm** and **hy = 40cm**. 📏
- **Array Width Limit**: The overall width of the array cannot exceed **240cm**.
- **Power Control**: Each modulated downlink stream can be driven with power ranging from **-30 dBm** to **+30 dBm**, in **1 dB** steps between timeslots. 

## 💡 Solution 
### Matlab Solution

<details open>

Initial proposal and source code along with optimizations have been made. Some of the visual effects of using linear scan antenna array can be found below: 

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

### CPP Solution and Simulation Output 🖥️

Overall antenna antenna patterns from using Linear Phase Array, 5 panels each of size 20 cm x 40 cm, scan-angle +40 degrees, frequency 1.9 GHz. CPP implementation using ImGUI and SFML using the following parameters in the `init.json` file

<p align="center">
<img src="images/overall.png?raw=true" width="800" alt="1000 meters x 1000 meters field">
</p>

#### Input 📥

| Parameter  | Value |
| :---         |          ---: |
| signal frequency  | 1.9 GHz |
| signal bandwidth  | 20 MHz | 
| symbolrate  | 3.84 Mega symbols |
| blockpersymbol  | 768 blocks/symbol |
| antenna-height  | 200 meters | 
| Grx of each client | -2 dB | 
| system-noise  | 5 dB | 
| rx count in the simulation | 15 | 
| tx count in the simulation | 3 | 
| SINR/SNR failure cutoff  | 24 dB | 
| tx theta-C direction  | 75, 90, 105 degrees | 
| antenna-counts for each tx  | 5, 3, 5 | 
| base_station_power_range  | -30 to +30 dBm | 
| base_station_scan_angle_range  | -90 to +90 degrees | 
| panel-spacing for each tx  | 30 cm, 40 cm, 30 cm | 
| antenna_dims  | 20 cm x 40 cm |

```python
arguments as follows:
-----------------------
      -o,--output_dir : C:\Users\User\Downloads\projects\Linear-Phased-Array
            -f,--file :
           --timeslot : 3
          --frequency : 1900000000.0
          --bandwidth : 20000000.0
         --symbolrate : 3840000.0
     --blockpersymbol : 768
             --height : 200
             --ms_grx : -2
       --system_noise : 5
      --base_stations : 3
    --mobile_stations : 15
          --timeslots : 1
             --slimit : 24.0
               --area : unknown
--base_station_theta_c : 75.0,90.0,105.0
--base_station_location : 250,200 500,180 750,200
--mobile_station_location : 250,500 350,500 450,500 550,500 650,500 750,500 250,600 350,600 450,600 550,600 650,600 750,600 450,325 500,325 550,325
--base_station_antenna_counts : 5,3,5
--base_station_power_range_dBm : -30,30
--base_station_scan_angle_range_deg : -90,90
    --antenna_spacing : 0.3,0.4,0.3
       --antenna_dims : 0.2,0.4
    --antenna_txpower : 30.0,29.0,28.0
         --scan_angle : 40.0,65.0,-40.0
       --ms_selection : 7,4,10
                   -g : false
           -q,--quiet : true
            -?,--help : false
```

<details open>
 
#### Output 🚀

``` 
timeslot: 3     power: 28.00    alpha: -40.00   placement: 650,600      cow_id: 2       sta_id: 10      sinr: 9.96
timeslot: 3     power: 28.00    alpha: -40.00   placement: 650,600      cow_id: 2       sta_id: 10      sinr: 9.96
timeslot: 3     power: 28.00    alpha: -40.00   placement: 650,500      cow_id: 2       sta_id: 4       sinr: 7.68
timeslot: 3     power: 28.00    alpha: -40.00   placement: 650,500      cow_id: 2       sta_id: 4       sinr: 7.68
timeslot: 3     power: 29.00    alpha: 65.00    placement: 350,600      cow_id: 1       sta_id: 7       sinr: 1.66
timeslot: 3     power: 29.00    alpha: 65.00    placement: 350,600      cow_id: 1       sta_id: 7       sinr: 1.66
timeslot: 3     power: 30.00    alpha: 40.00    placement: 350,600      cow_id: 0       sta_id: 7       sinr: -1.99
timeslot: 3     power: 30.00    alpha: 40.00    placement: 350,600      cow_id: 0       sta_id: 7       sinr: -1.99
timeslot: 3     power: 29.00    alpha: 65.00    placement: 650,500      cow_id: 1       sta_id: 4       sinr: -8.56
timeslot: 3     power: 29.00    alpha: 65.00    placement: 650,500      cow_id: 1       sta_id: 4       sinr: -8.56
timeslot: 3     power: 30.00    alpha: 40.00    placement: 650,600      cow_id: 0       sta_id: 10      sinr: -11.87
timeslot: 3     power: 30.00    alpha: 40.00    placement: 650,600      cow_id: 0       sta_id: 10      sinr: -11.87
timeslot: 3     power: 30.00    alpha: 40.00    placement: 650,500      cow_id: 0       sta_id: 4       sinr: -16.61
timeslot: 3     power: 30.00    alpha: 40.00    placement: 650,500      cow_id: 0       sta_id: 4       sinr: -16.61
timeslot: 3     power: 28.00    alpha: -40.00   placement: 350,600      cow_id: 2       sta_id: 7       sinr: -19.80
timeslot: 3     power: 28.00    alpha: -40.00   placement: 350,600      cow_id: 2       sta_id: 7       sinr: -19.80
timeslot: 3     power: 29.00    alpha: 65.00    placement: 650,600      cow_id: 1       sta_id: 10      sinr: -31.49
timeslot: 3     power: 29.00    alpha: 65.00    placement: 650,600      cow_id: 1       sta_id: 10      sinr: -31.49

```

#### Simulation Output Explanation 🔭 

Since there are **3 transmitters** and **15 stations**, using **single frequency** and no **MIMO** at any given timeslot, we split the exchange over 5 timeslots (15 stations ÷ 3 transmitters = 5 timeslots). 

At any timeslot:
- Each transmitter produces **1 line of output** due to the binding info from the `init.json` test file.
- In a specific timeslot (for example, `--timeslot: 3`), we should see **3 lines of output**—one for each transmitter.

However, you may notice more than 3 lines of output (e.g., **18 lines**)! This happens because:

1. I use **permutations** of various bindings between TX (transmitters) and RX (receivers) stations. 
2. Additionally, I permute the unique TX powers. Instead of fixed power levels (e.g., 30, 30, 30 dBm), I vary the power slightly (e.g., **30, 29, 28 dBm**) for different combinations. 

This approach generates more **permutations** of signal-to-noise ratio (SNR) data, using parameters that may not have been considered in the original driving script acting as a wrapper.

#### Why So Many Lines of Output? 🤔
These extra lines are due to the additional permutations and variations in TX power and bindings. In essence, we are exploring **more combinations** of potential scenarios between TX and RX stations to provide more comprehensive **SNR data**.

#### Ready for Deep Q-Learning Integration 🧠
This output can now be integrated into **Deep Q-Learning Networks (DQN)** or **Reinforcement Learning (RL)** models. For those interested, check out the relevant scripts inside the `src` directory. 📂

#### Key Directory:
- `src/`: Contains **Reinforcement Learning** scripts for machine learning applications using this dataset. 🤖

</details>

#### Summary
- **3 transmitters** ↔ **15 stations** divided across 5 timeslots.
- Unique TX power variations create more SNR permutations.
- **Output** is ready for **DQN** and **RL** applications. 


#### Plot GUI heatgraphs from each tranmitters

<details>
Tx1
<p align="center">
<img src="images/transmitter_0.png?raw=true" width="800" alt="1000 meters x 1000 meters field">
</p>

Tx2
<p align="center">
<img src="images/transmitter_0.png?raw=true" width="800" alt="1000 meters x 1000 meters field">
</p>

Tx3
<p align="center">
<img src="images/transmitter_0.png?raw=true" width="800" alt="1000 meters x 1000 meters field">
</p>
</details>


## Building and Integrating 🛠️

Default is release build with -O3 optimization

```
mkdir build
cd build
cmake ..
cmake --build .
```

### Debug build w/ debug symbols 🛠️
Change `cmake ..` to `cmake -DCMAKE_BUILD_TYPE=Debug ..`



