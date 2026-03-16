# Scripts that generate the results in SAASI manuscript

## Required packages/software:
R (Version 4.4.1)
Python (Version 3.13.7)

TreeTime(Python): https://treetime.readthedocs.io/en/latest/

PastML(Python): https://pastml.pasteur.fr/

Typical install time: less than 1min.

## R library used
tidyr
stringr
treeio
tidyverse
ggalluvial
ggpubr
irr
stats4
ape
deSolve
strap
ips
ggplot2
ggtree
reshape2
ggstance
diversitree
ggimage
patchwork
janitor
phytools
dplyr
rstatix
jsonlite
rootSolve

## Content
### Code:
Code for generating the figures in the manuscript.

#### H5N1:
Codes that conduct phylogeography analysis for H5N1.

#### Results:
RDS files for the results.

#### Simulation:
Codes that conduct simulation studies.

### Data:
H5N1 tree used in the analysis. See original paper (Nguyen et al. 2025, DOI: 10.1126/science.adq0900
) and GitHub: https://github.com/flu-crew/dairy-cattle-hpai-2024/tree/main.

### Figures:
Figures in the manuscript.

### SAASI:
The Ancestral State Inference algorithm that is described in the manuscript.

## Running Shell script:

In Code/Simulation/Shell, the user can find the shell script that generates the 1,000 simulated trees and performs Ancertal State Inference methods (ace,simmap,TreeTime,PastML, SAASI).
Before running the script, the following libraries in R need to be installed: saasi, diversitree, ape, phytools, tidytree, and desolve. Furthermore, TreeTime and PastML should also be installed in Python. 

To run the shell script, in the command line, run:
```
./xxx.sh 1000 5000 5000 folder_name
```
The first argument is the number of simulations.
The second argument is the maximum running time.
The third argument is the maximum number of tips in the present day.
The fourth argument is the name of the output folder.

The expected output is 1,000 folders, each folder contains all the Ancerstral State Inference methods (ace,simmap,PastML,TreeTime,SAASI) results, running times, and the true simualted phylogeny. The expect running time for one simulation (one folder) is expected to be between 1min and 5min.


