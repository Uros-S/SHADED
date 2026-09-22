# SHADED
SHADED is a **S**avitzky-Golay and **H**omogeneous-differentiator based **A**utomatic **DE**noising and **D**ifferentiation method for effective denoising and estimation of derivatives of a recorded signal up to an arbitrary order. It enables attractor reconstruction via differential embedding from noisy time series data.

In this repo, find the code used to generate the figures in Main Text and Supplementary Material of the paper "Automatic denoising and differentiation based on Savitzky-Golay filtering and Homogeneous Differentiators for attractor reconstruction via differential embedding" by U. Sutulovic, D. Proverbio, R. Katz, and G. Giordano.   
**ArXiv version** of the paper can be found at the link http://arxiv.org/abs/2609.18631 .  
**Supplementary material** for the paper can be found at the Zenodo link https://doi.org/10.5281/zenodo.22796837 .

### Credits
The code is released under a GNU General Public License. Please cite the original reference if you reuse the code, even partially, or its results.

## Usage

The code is fully self-contained and runs on Matlab (tested with Matlab R2023b). Download the whole repo with its subfolders to run it properly.  
Requires the `Global optimization` Matlab toolbox.

### Files and folders
- `main_theoretical_models.m`: main file. Run this one to reproduce the results.
- `utils`: folder containing companion functions. Check each file for a detailed description.
- `systems`: folder with model equations.

## Reproduce results
Run `main_theoretical_models.m` to reproduce the results and figures of the Main Text. 
Below, the steps and settings to reproduce each figure of the Main Text.  
For future works, you can refer to the settings below as blueprints.

### Figure 1 in Main Text

* In "%% Selection of models, noise and result to be obtained" insert the settings below and leave "%% Parameter settings for Differentiator, Savitzky-Golay filter and computational complexity" unchanged:
 
model = 1;

modality = 1;

noise_type = 1;

additive_noise = 1;

noise_var = 5;

n_derivatives = 2;

persistence = true;

### Figure 2 in Main Text

* In "%% Selection of models, noise and result to be obtained" insert the settings below and in "%% Parameter settings for Differentiator, Savitzky-Golay filter and computational complexity" insert n_d = 0;.
 
model = 1;

modality = 2;

noise_type = 1;

additive_noise = 1;

noise_var = 5;

n_derivatives = 2;

persistence = true;


### Figure 3 in Main Text

* In "%% Selection of models, noise and result to be obtained" insert the settings below and in "%% Parameter settings for Differentiator, Savitzky-Golay filter and computational complexity" insert n_d = 0; and L = 540;.
 
model = 1;

modality = 3;

noise_type = 1;

additive_noise = 1;

noise_var = 5;

n_derivatives = 2;

persistence = true;

### Figure 4 in Main Text

* In "%% Selection of models, noise and result to be obtained" insert the settings below with the desired model and noise intensity (see Section 3.1). In "%% Parameters setting for methodology with guess parameter values" insert n_d = 2; and L_guess and fl_guess as specified each panel (note that fl_guess needs to be an odd natural number, use formula (11) to convert the T values reported in the panels).

modality = 1;

noise_type = 1;

additive_noise = 1;

n_derivatives = 2;

persistence = true;


## License 
This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

The full text of the GNU General Public License can be found in the file "LICENSE.txt".
