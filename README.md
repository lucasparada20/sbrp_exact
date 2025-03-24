# The DL Shaped method for the SBRP
An implementation of the Disaggregated Integer L-shaped (DL-shaped) method for the stochastic bicycle repositioning problem and using the Cplex library.

* [Technical Report](https://www.cirrelt.ca/documentstravail/cirrelt-2024-26.pdf)

## Building the code in Linux

1. Clone the repository and add executable permission to a script that will call CMake for you:

```shell
git clone https://github.com/lucasparada20/sbrp_exact.git
cd sbrp
chmod u+x cmake_script_exact.sh
```
2. The compiler needs an absolute path to your installed Cplex library. To provide the path, go into src/CMakeLists.txt and edit the following line:

```cmake
set(CPLEX_DIR "/some/path/to/Cplex")
```

for example, in my Ubuntu environment, the absolute path is as follows:

![nano src/CMakeLists.txt](/images/image.png)

3. Build the code by typing:

```bash
./cmake_script_exact.sh
```

## Running the code

Inside the sbrp_exact directory, you will find a script run.sh with sample command line calls. The format is:

* instance_file : the instance to solve.
* instance_type : dins or pcg. dins denotes instances from [Dell'Amico et al., (2018)](https://doi.org/10.1016/j.trb.2018.10.015). pcg denotes the instances of our technical report. For the dins instances, resolution begins by a previously computed upper bound from an implementation of the Adaptative Large Neighborhood Search algorithm. The bounds are in sbrp/instances_dins/all_upper_bounds.txt and are read in the UbManager class from UbManager.h.
* epsilon, delta : the parameters for cost computations.
* opt_cuts : defines the optimality cuts used for the resolution. 1 for P&L cuts, 2 for Benders cuts, and 3 for hybrid cuts. The multi-cut algorithm can only be called with Benders cuts.
* algorithm : dl or m. dl denotes the DL-shaped method, and m represents the integer L-shaped multi-cut algorithm of [Dell'Amico et al., (2018)](https://doi.org/10.1016/j.trb.2018.10.015).

In run.sh, the first example calls the executable with the instance Chicago_20_1.txt from [Dell'Amico et al., (2018)](https://doi.org/10.1016/j.trb.2018.10.015), using epsilon = delta = 0.2, Benders optimality cuts, and the DL-shaped method. 

```bash
build/exec_exact instance_file=instances_dins/Chicago_20_1.txt epsilon=0.2 delta=0.2 opt_cuts=2 instance_type=dins algorithm=dl
```
The second example in run_exact.sh calls the executable with the instance ssbrp_30_20_0.txt from the [Technical Report](https://www.cirrelt.ca/documentstravail/cirrelt-2024-26.pdf), using epsilon = delta = 0.2, hybrid optimality cuts and the DL-shaped method. 

```bash
build/exec_exact instance_file=instances_pcg/ssbrp_30_20_0.txt epsilon=0.2 delta=0.2 opt_cuts=3  instance_type=pcg algorithm=dl
```
## Ouput
Upon calling the examples in run_exact.sh, you should see the following output.

### Chicago_20_1.txt

![Chicago_20_1.txt](/images/Chicago_20_1.png)

### ssbrp_30_20_0.txt

![ssbrp_30_20.txt](/images/ssbrp_30_20_0.jpg)
