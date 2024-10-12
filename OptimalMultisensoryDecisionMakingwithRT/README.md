# OptimalMultisensoryDecisionMakingwithRT

This repo constains some of the code use for model fitting in [Drugowitsch, DeAngelis, Klier, Angelaki & Pouget (2014)](https://doi.org/10.7554/eLife.03005) and [Drugowitsch, DeAngelis, Angelaki & Pouget (2015)](https://doi.org/10.7554/eLife.06678).

**Please note that this code is research code and for your information only. It will not run without additional adjustments - in particular adjustments the data loading routines. The below should provide sufficient information to perform these adjustments.**

The files have been tested with MATLAB R2022b on OS X Ventura 13.6, but should also work on older (at least back to 2014) or newer MATLAB and OS versions. It requires the MATLAB Optimization Toolbox to be installed.

## Code file structure

`mcmc_fit_LLH` - MATLAB files to fit different models to behavior

`data` - MATLAB files to load data, and the data itself (not provided)

`shared` - files shared across various scripts

## Installation

Download or clone all files in the repository, and add the following folders to the MATLAB path

```
data
shared
shared/inference
shared/ddm
```

The `shared/ddm` folder contains an (old) copy of the [Diffusion model toolset](https://github.com/DrugowitschLab/dm), which is a collection of C++ programs with MATLAB and Python interfaces. The toolset need to be complied into binaries to be used. While the repository already provide some binaries, it only does so for some operating systems. Please see the [Diffusion model toolset MATLAB interface documentation](https://github.com/DrugowitschLab/dm/tree/v0.1.1) for how to compile these binaries for your operating system, if needed.

Finally, navigate to the `mcmc_fit_LLH` folder to perform the fits.

## Usage

### Model fitting

The `mcmc_fit_LLH` folder contains all the code to sample from the model parameter posterior for each of the respective models, as well as the model code.

* `fit_model.m`: the main fitting script. It loads the subject data (using `load_condi_data.m`, which in turn calls the `data/load_data.m` script mentioned further below), initializes the chosen mode, and samples from the parameter posteriors. The posterior samples are written to the `mcmc_fit_LLH/fit_data` folder. If posterior samples are already present, then it loads these samples first, and adds to them. In Drugowitsch et al. (2014) We used 44000 samples to fit each of the individual models.

  For example, if you want to draw 10000 posterior samples for subject `subj01` using model `mf_sepk`, you would call `fit_model` by
        ```
        fit_model('subj01', 'mf_sepk', 10000)
        ```
  `fit_model` has some additional options that are described in the script's docstring.

* `mf_...m`: the different model that we tried. The only we finally used in Drugowitsch et al. (2014) is `mf_sepk.m` which fits each of the sensitivities separately. One can then predict the multimodel sensitivity from the unimodel ones, as done in the next script. Initial parameter values and bounds on those are hard-coded in each of the model scripts.

* `plot_k_relation.m`: this script might be of interest, as it was used it to generate the sensitiviy plot in Drugowitsch et al. (2014) (plot title 'all subjects, average of mode + SEM, ...' in the script). In particular, the average sensitivities predicted from the unimodal conditions are stored in `avg_kpred(1,:)`, and the sensitivities fitted in the combined condition are stored in `avg_kcomb(1,:)`. These averages are computed from the highest-likelihood fitted parameter for each participant, using the `mf_sepk` model. This script will only run if posterior samples for all subjects have been compted.

There are a bunch of other plotting scripts in this folder. They are provided as additional information for how the resulting data / model parameters can be read out and used.


### Data organization and loading

All the data lives in the `data` folder.

* `load_data.m`: that's the main workhorse of this folder, loading and pre-processing the requested dataset. It is called by `mcmc_fit_LLH/load_condi_data.m`. No datasets are provided, but the script should provide some idea of the format of the data expected by `load_condi_data.m`.
