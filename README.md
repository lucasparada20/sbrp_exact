# The-DL-Shaped-method-for-the-SBRP
An implementation of the Dissagregated Integer L-shaped (DL-shaped) method for the stochastic bicycle repositioning problem

* [Technical Report](https://www.cirrelt.ca/documentstravail/cirrelt-2024-26.pdf)

## Building the code

1. Clone the repository and add executable permission to a script that will call CMake:

```shell
git clone https://github.com/lucasparada20/sbrp.git
cd sbrp
chmod u+x cmake_script.sh
```
2. The compiler needs an absolute path to your installed Cplex library. To provide the path, go into src/CMakeLists.txt and edit the following line:

```cmake
set(CPLEX_DIR "/some/path/to/Cplex")
```

for example, in my Ubuntu environment, the absolute path is as follows:

![nano src/CMakeLists.txt](/images/image.png)

3. Build the code by typing:

```bash
./cmake_script.sh
```

## Running the code

Inside the sbrp directory you will find a script run.sh with sample command line calls. The format is:

* instance_file : the instance to solve.
* instance_type : dins or pcg. dins denotes instances from [Dell'Amico et al., (2018)](https://doi.org/10.1016/j.trb.2018.10.015). pcg denotes the instances of our technical report.
* epsilon, delta : the parameters for cost computations.
* algorithm : dl or m. dl denotes the DL-shaped method and m denotes the interger L-shaped multicut algorithm of [Dell'Amico et al., (2018)](https://doi.org/10.1016/j.trb.2018.10.015).

In run.sh, the first example calls the executable with the instance from Dell'Amico et al., (2018) Chicago_20_1.txt, using epsilon = delta = 0.2, Benders optimality cuts and the DL-shaped method. 

```bash
build/exec_exact instance_file=instances_dins/Chicago_20_1.txt epsilon=0.2 delta=0.2 opt_cuts=2 instance_type=dins algorithm=dl
```


