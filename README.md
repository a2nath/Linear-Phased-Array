<p><video src="https://github.com/user-attachments/assets/df660384-1153-43e5-af5e-90b11e6b7a32" width="100%" height="auto" autoplay muted playsinline></video>
</p>

# README #

## Table of Contents
-   [Introduction](#-introduction)
-   [Example Scenario](#-example-scenario)
-   [Parameters](#-parameters)
-   [Matlab Solution](#-solution)
-   [CPP Solution and Simulation Output](#cpp-solution-and-simulation-output-%EF%B8%8F)
-   [Building and Integrating](#building-and-integrating-%EF%B8%8F)


## 👥 Introduction
In the age of rapid digital transformation, reliable wireless communication is essential for connecting people, communities, and critical services worldwide. *Cellular on Wheels (COW)* 🐄 technology, also known as mobile base stations, plays a vital role in extending cellular coverage in areas where traditional, fixed base stations are insufficient, or infrastructure is lacking. This project leverages a linear phased array for enhancing the functionality and efficiency of COW units, enabling robust and flexible coverage solutions.

### 🐮 Why Cellular on Wheels (COW) Matters
COWs bring essential mobile connectivity in situations where permanent installations are impractical. Here’s why they make a meaningful difference:

* ### Emergency Response & Disaster Recovery

  * **COWs** provide rapid, temporary cellular coverage in disaster-stricken areas, supporting emergency responders, enabling critical communications, and reconnecting displaced populations.
  * **Example**: After natural disasters like hurricanes, earthquakes, or wildfires, COWs can be deployed within hours to restore cellular connectivity, allowing for immediate coordination and life-saving operations.

* ### Event & Crowd Coverage

  * Large events, such as **concerts**, **sports** such as the Olympics, or **festivals**, experience significant surges in cellular demand, often exceeding local infrastructure capacity.
  * COWs augment the network in these high-demand situations, improving connection reliability and data speeds for thousands of attendees.

* ### Expanding Coverage in Rural & Remote Areas:

  * COWs can provide temporary coverage in hard-to-reach or underserved areas, allowing for cellular service where permanent towers may be cost-prohibitive.
  * They are also essential for extending services during infrastructure projects or for seasonal agricultural activities.

* ### Technical Advantages of Using COWs with Phased Arrays:

  * **Dynamic Coverage**: Phased array technology enables COWs to dynamically steer beams, maximizing signal strength and efficiently targeting areas of demand.
  * **Enhanced Data Rates**: By focusing beams precisely, phased arrays reduce interference and improve data throughput, especially in high-density settings.
  * **Faster Deployment & Adaptability**: Lightweight and modular, COWs with phased arrays can be deployed within hours and adjusted to provide optimal coverage as needed.
  * **Improved Energy Efficiency**: Phased arrays in COW systems can direct energy precisely where needed, lowering power consumption compared to traditional omnidirectional antennas.

### 🏕️ Real-World Impact

  * **Increased Connectivity**: Studies show that even a temporary COW deployment can enhance network performance by up to 50% in high-traffic scenarios.
  * **Operational Efficiency**: Phased array COWs provide up to 30% better spectrum efficiency by reducing signal interference, crucial for rural areas with limited bandwidth.

Following this, you’ll find a specific example that demonstrates how a phased array setup can enhance the flexibility and effectiveness of COW technology, providing tangible improvements in coverage and capacity for real-world applications.


## 🔥 Example Scenario

There are **3 "Cellular on Wheels (COWs)"** portable base stations, each equipped to drive a **linear scan antenna array** with up to **8 elements**. The layout consists of **12 non-VIP** spectator zones and **3 VIP** spectator zones. The goal 🎯 is to ensure all clients in the simulation receive at least **24 dB of SINR** ("SNR") signal quality or at least **3 VIP** zones are above the cut off if there exists no solution to the specific scenario. 📶

More details on the linear scan antenna array as [follows.](https://en.wikipedia.org/wiki/Phased_array)

<p align="center">
<img src="images/280px-Phased_array_animation_with_arrow_10frames_371x400px_100ms.gif?raw=true" width="400" alt="How it works">

In this project, **theta-C° angle** represents the direction of the transmitter relative to EAST (0°), as shown in the image. This is the direction in which the transmitter is physically pointing. On the other hand, **scan angle°, θ** vary between **-90°** and **+90°**, relative to **theta-C° angle**, or from the perpendicular axis to the array of antennas.

</p>
<p align="center"><img src="images/Phasearray.gif?raw=true" width="400" alt="How it works">
</p>


## 🔧 Parameters
**Frequency Channel**: 20 MHz wide, **f-center = 1900 MHz**.
- Each **COW** can mount and **linearly scan** up to **8 antennas** driven by a common downlink modulator signal. 🛰️
- **Antenna Dimensions**: Each antenna has **hx = 20cm** and **hy = 40cm**. 📏
- **Array Width Limit**: The overall width of the array cannot exceed **240cm**.
- **Power Control**: Each modulated downlink stream can be driven with power ranging from **-30 dBm** to **+30 dBm**, in **1 dB** steps between timeslots. 


<details>

 ## 💡 Solution 
### Matlab Solution

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

</details>

### CPP Solution and Simulation Output 🖥️

Overall antenna antenna patterns from using Linear Phase Array, 5 panels each of size 20 cm x 40 cm, scan-angle +40 degrees, frequency 1.9 GHz. CPP implementation using ImGUI and SFML using the following parameters in the `init.json` file

#### Input 📥

<div align="center">
      
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

</div>

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



<details open>
 
## Plot GUI heatgraphs and example interference pattern 

### Scenario 1

_Perspective_ is from the middle transmitter, TX2. It is having the most influence on all 15 clients because the TX power is +29 dBm as shown in the control panel to the right. TX1 and TX3 are set to 0 dBm. The scan angles are 0 for all 3 TXs, therefore the radiation pattern is symmetric and can reach handsets on the left-side and handsets on the right-side equally. 

**Outcome**

Bind TX1 single stream link to STA 13, 14, or 15 across 3 timeslots for _highest_ SNR in the range of 48-72 dB. However, we can also use STA 3, 4, 8, or 11 since all of them have 21-48 dB of SNR. 

<p align="center">
<img src="images/state1.PNG?raw=true" width="800" alt="1000 meters x 1000 meters field">
</p>

### Scenario 2

_Perspective_ is from the left transmitter or TX1 with a strong signal of +30 dBm, like TX3. Control panel on the right will show the output power of the Linear Antenna Arrays, TX1, TX2 and TX3. Here TX2 is set to the lower +10 dBm output so its _influence_ is much lower on the heat graph. 

The scan angles are 0 for all 3 TXs, therefore the radiation pattern is symmetric for each transmitter looking at them individually. However, with the interference pattern the _perspective_ of TX1 is very different. 

Notice that the signal-to-noise ratio closer to TX3 is < 7 dB (blue-ish), because both TX1 and TX3 having the same energy output causes them to clash and interference is at a maximum from TX1's _perspective_. Vice-version, TX3's power closer to TX1 will have maximum interference, therefore the SNR (not shown) is lowest in the blue.

**Outcome**

Bind TX1 single stream link to STA 7, 10, or 12 across 3 timeslots for _highest_ SNR in the range of 34-62 dB. However, we can also use STA 1, 8, 9, or 11 since all of them have 21-48 dB of SNR. 

<p align="center">
<img src="images/state2.PNG?raw=true" width="800" alt="1000 meters x 1000 meters field">
</p>

### Scenario 3

_Perspective_ again, is from the left transmitter or TX1 with a strong signal of +30 dBm, like TX3. The only difference since last image is the location of TX3. Control panel on the right will show the output power of the Linear Antenna Arrays, TX1, TX2 and TX3. Here TX2 is set to the lower 10 dBm output so its _influence_ is much lower on the heat graph. 

The scan angles are 0 for all 3 TXs, therefore the radiation pattern is symmetric for each transmitter looking at them individually. In this scenario both TX1 and TX3 are placed almost at the same Y value to show the interference more clearly; see how very little places in the grid have SNR at 48+ dB (bright red)? This is because of the interference coming from TX3. Once switching to TX3's _perspective_, the heat map will be very similar.

**Outcome**

There are more STAs if not all, closer to 7 dB SNR. 

<p align="center">
<img src="images/state3.PNG?raw=true" width="800" alt="1000 meters x 1000 meters field">
</p>

### Scenario 4

_Perspective_ is from the middle transmitter, TX2 with a strong signal of +35 dBm therefore the signal is much stronger around STA 13, 14, and 15. TX1's power has been reduced to -10 dBm but TX3 still has +30 dBm output power. Control panel on the right will show the output power of the Linear Antenna Arrays, TX1, TX2 and TX3.

The scan angles are 0 for all 3 TXs, therefore the radiation pattern is symmetric for each transmitter looking at them individually. In this scenario we demonstrate the SNR as the distance increases, and a lot of the original antenna gain is also visible if it weren't for all the interference seen from other transmitters. See Matlab solution above; we call this "clown fingers".

**Outcome**

The effect of antenna [3-dB gain](https://en.wikipedia.org/wiki/Half-power_point) can be seen at STA 2, 3, 6, 7 even though the last two: 6 and 7 are closer to TX3. In the future I need to show the individual plots to make it clearer and prove this stuff better.

<p align="center">
<img src="images/state4.PNG?raw=true" width="800" alt="1000 meters x 1000 meters field">
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



