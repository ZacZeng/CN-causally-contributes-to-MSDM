# CN-causally-contributes-to-MSDM
## Aim
To replicate the results in paper titled as ["Distinct neural manifolds and critical roles of primate caudate nucleus in multimodal decision-making"](https://www.biorxiv.org/content/10.1101/2024.09.03.610907v1).

## Software 
We use MATLAB 2017a for data analysis.  
RNN simulations are conducted using pytorch-cuda 12.1.

## Installation
Download the [dataset](https://zenodo.org/records/13923317) and gunzip it.  
Download this code package. 

## Instructions
1. Add all MATLAB files to the default path in MATLAB.
2. Print `Group_GUI` in the Command Window to get the GUI named as **Group Data Analyzer**.
   * Click `Read Data Hub` to load the contents in file **DataHub.xlsm**.
   * After reading, click `CD Heading` to load neurophysiological data, 'Microstimulation' for electric stimulation data, or `Psychometric` for behavioral data.
   * Click the analysis protocol in the left panel under **Figures** and then specific item in the right.
   * Choose the animal and cells you wanted in the panel at the southwest corner.
   * Click `GO`. You can get the CN results in Figure 1, Figure 2, Figure 4, Figure 6, Figure S1, Figure S3, figure S7.
![Figure 1](https://github.com/ZacZeng/CN-causally-contributes-to-MSDM/blob/main/figure/Group_GUI.png)
