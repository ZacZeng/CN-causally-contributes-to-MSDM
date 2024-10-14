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
1. Add all files to the MATLAB path.
2. Print `Group_GUI` in the Command Window to get the GUI named as **Group Data Analyzer**.  **NOTE: You need to change the default directories in some scripts into the correct local position in your computer!!!**
   * Click `Read Data Hub` to load the contents in file **DataHub.xlsm**.
   * After reading, click `CD Heading` to load neurophysiological data, 'Microstimulation' for electric stimulation data, or `Psychometric` for behavioral data.
   * Click the analysis protocol in the left panel under **Figures** and then specific item in the right.
   * Choose the animal and cells you wanted in the panel at the southwest corner.
   * Click `GO`. You can get the CN results in Figure 1, Figure 2, Figure 4, Figure 6, Figure S1, Figure S3, figure S7.
     
![Group GUI](https://github.com/ZacZeng/CN-causally-contributes-to-MSDM/blob/main/figure/Group_GUI.png)

3. For population analysis in cortical regions (LIP, FEF, MSTd), you can run `PopulationAnaly_ZZ`.
4. Inactivation results in Figure 5 and S6 can be found by running `InactiveCompare` and `ColorSelectionCompar`.
5. Run `dPCA_RNN` to get the neural states of RNN units in Figure 3.
6. For GDDM fittings of behavioral data in RT-version multisensory heading discrimination task, source code is from [DrugowitschLab](https://github.com/DrugowitschLab/OptimalMultisensoryDecisionMakingwithRT). For running,
   * First make the Current Folder in MATLAB as [`dm-0.1.1`](https://github.com/ZacZeng/CN-causally-contributes-to-MSDM/tree/main/OptimalMultisensoryDecisionMakingwithRT/shared/ddm/dm-0.1.1).
   * Print `fit_model('Monkey_M','mf_sepk',100)` in the Command Window.
   * Set a Breakpoints in the script of `fit_model` at line ~ 135.
   * Print `plot_fits(subj_id,model_fn)`, you can get the fitting results in Figure 4.

## Acknowledgement
Group Analysis pipeline is established by [Han Hou](https://github.com/hanhou/Labtools).  
RNN simulations are performed by Ce Zhang.  
dPCA package is from [Here](https://github.com/machenslab/dPCA).  
GDDM package is from [DrugowitschLab](https://github.com/DrugowitschLab/OptimalMultisensoryDecisionMakingwithRT).
