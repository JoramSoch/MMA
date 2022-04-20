
# MMA

### Voxel-wise multi-modal analyses for MRI data using SPM

This code allows to perform [SPM-style](https://www.fil.ion.ucl.ac.uk/spm/) GLM analyses across multiple modalities (e.g. structural MRI, functional MRI and resting-state fMRI), using one of these modalities as the dependent variable and other modalities as voxel-wise covariates, possibly influencing the variable of interest. Model parameters and results are stored in the [SPM.mat](https://github.com/spm/spm12/blob/master/spm_spm.m#L157-L217) format, such that analysis outcomes can be viewed using SPM's results interrogation.


## Requirements

For using the toolbox, the following software should be installed:
- [MATLAB](https://de.mathworks.com/products/matlab.html) (e.g. R2020a)
- [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) (Revision 7771 as of 14/01/2020)
- [MACS](https://github.com/JoramSoch/MACS) (later than commit [67572cc](https://github.com/JoramSoch/MACS/commit/67572ccf95575ad5968d78566456ba666c44edfe))


## Getting started

The toolbox and its manual use the following terminology:
- dependent modality (DM) = data that should be used as the dependent variable in the multi-modal analysis
- independent modalities (IMs) = data that should be used as voxel-wise covariates in the multi-modal analysis

In order to perform a multi-modal analysis (MMA) using this toolbox, proceed as follows:
1. Run a factorial design* for images from the dependent modality (e.g. one-sample t-test, two-sample t-test, one-way ANOVA, full factorial). The results of this analysis are stored into a directory that is hereafter referred to as `DM_dir`.
2. Run the same factorial design* for images from an independent modality (to ensure that subjects are in the same order in all models). The results of this analysis are stored into a directory that is hereafter referred to as `IM_dir`.
3. Create a directory into which results of the multi-modal analysis are to be saved and that is hereafter referred to as `MMA_dir`. Then, run the following commands:

```matlab
SPM_mat_DM = strcat(DM_dir,'/','SPM.mat');
SPM_mat_IM ={strcat(DM_dir,'/','SPM.mat')};
con_mat_IM ={eye(p)}; % p = number of groups in the factorial design
mma_mma(SPM_mat_DM, SPM_mat_IM, con_mat_IM, MMA_dir);
```

4. This will create a voxel-wise general linear model (GLM) using the dependent modality as measured signal and the independent modality as group- and voxel-wise covariates. Model parameters will be estimated including [SPM's temporal and spatial corrections](https://github.com/JoramSoch/MMA/blob/main/mma_mma.m#L100-L110) and all contrasts from the original factorial designs will be transferred to the new model.
5. Open SPM, click "Results", select the `SPM.mat` in `MMA_dir` and access GLM results in the usual way.

* For more details, e.g. (i) how to integrate more than one independent modality `SPM_mat_IM` into the model, (ii) how to use the input variable `con_mat_IM` more efficiently and (iii) how independent modalities are treated when adding them to the model, see the [instructions in the main function of the toolbox](https://github.com/JoramSoch/MMA/blob/main/mma_mma.m#L35-L98).


## Additional functions

Once that a multi-modal analysis has been run, the following further actions can be performed:

```matlab
load(strcat(MMA_dir,'/','MMA.mat'));
```

1. Create and calculate additional contrasts (see [mma_cons.m](https://github.com/JoramSoch/MMA/blob/main/mma_cons.m#L6-L9) for input variable `xCon`):

```matlab
mma_cons(MMA, xCon, false);
```

2. Plot contrast estimates and 90% confidence intervals (similar to [spm_plot_con_CI.m](https://github.com/JoramSoch/spm_helper/blob/master/spm_plot_con_CI.m)):

```matlab
mma_plot_con_CI(MMA, con_ind, [x, y, z], 0.1, true);
```

3. Plot regression lines with 90% confidence bands (similar to [spm_plot_reg_CI.m](https://github.com/JoramSoch/spm_helper/blob/master/spm_plot_reg_CI.m)):

```matlab
mma_plot_reg_CI(MMA, reg_ind, [x, y, z], 0.1, true);
```


## MMA pipeline

The following drawing illustrates the pipeline of (1) creating separate factorial designs for dependent and independent modalities, (2) specifying how groups are collapsed when integrating the independent modalities as voxel-wise covariates and (3) and merging DM and IMs into a single MMA model:

<img src="https://raw.githubusercontent.com/JoramSoch/MMA/main/MMA_Manual/MMA_Pipeline.jpg" alt="MMA Pipeline" width=600>


## MMA examples

The following figures are outputs from exemplary applications of "mma_plot_con_CI.m" and "mma_plot_reg_CI.m":

* contrast estimates with confidence intervals:
<img src="https://raw.githubusercontent.com/JoramSoch/MMA/main/MMA_Examples/mma_plot_con_CI.png" alt="MMA: plot contrast estimates with CIs" width=600>

* regression lines with confidence bands:
<img src="https://raw.githubusercontent.com/JoramSoch/MMA/main/MMA_Examples/mma_plot_reg_CI.png" alt="MMA: plot regression lines with CBs" width=600>
