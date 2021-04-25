-------------
**RCD: Rotation Coordinate Descent**
-------------

Rotation Coordinate Descent (RCD) is a fast rotation averaging algorithm that achieves global optimality under mild noise conditions on the noise level of the measurements.

## Quick Start
This demo runs in MATLAB, with **RCD** and **RCDL** compiled in C++.
We tested under version R2020a on systems with:
- macOS Catalina
- Ubuntu 18.04

We provide 2 demonstration programs for SfM and SLAM settings:
- **demo_rcd.m** :
      - RCD on SfM camera graphs with *n* =1000, 2000 and 3000 cameras.
      - This demo finishes in about 1 minute.

- **demo_rcdl.m** : RCDL on torus3D and grid3D dataset.
      - RCDL on SLAM camera graphs (1) torus3D (2) grid3D.
      - This demo finishes in about 1 minute.

## Setup ##
Dependencies:
   1. CMake 3.0 or later **required** [cmake installation link](https://cmake.org/install/)
      - *MacOS*   
      ```brew install cmake```

      - *Ubuntu*  
      ```sudo apt-get install cmake```

   2. SuiteSparse **required**
      - *MacOS*
        ```brew install suite-sparse```

      - *Ubuntu*
        ```sudo apt-get install libsuitesparse-dev```


## Running demo ##

### Building RCD and RCDL libraries
 We provide a script build.sh to build **RCD** and **RCDL**.
   Please make sure you have installed all required dependencies (see Section 1).
   Execute
   ``` 
   chmod +x build.sh
    ./build.sh
   ```
   which will create the executables **RCD** and **RCDL** in *bin* folder.
   
### Running RCD
  Run **demo_rcd.m** in MATLAB.
  Here, we provide a demonstration on 3 graph sizes, where *number of nodes* = 1000, 2000 and 3000 and *graph density* = 0.4.
  ![Output](https://github.com/sfchng/Rotation_Coordinate_Descent/blob/b3abd9206df1e162020f12279c2cdbf5df65b5b2/n_1000.pdf)


### Run "demo_rcd.m" and "demo_rcdl.m" in MATLAB.

1. Run "demo_rcd.m" for the demonstration of RCD.
2. Run "demo_rcdl.m" for the demonstration of RCDL.
