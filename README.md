**README**

Developed by Suzanne Shoffner-Beck and Kade Wong at the Arnold Lab at University of Michigan as an expansion on Melissa Lemke's model (see https://github.com/melissalemke/Ab_FcR_ODE)

**Overview**

This repository contains MATLAB code for an ODE-based mathematical model investigating the simultaneous activation of FcγRIIa and FcγRIIIa in different tissue environments (blood, rectal, and penile). The model evaluates antibody-dependent cellular phagocytosis (ADCP) and antibody-dependent cellular cytotoxicity (ADCC) through various computational approaches, including sensitivity analyses and surface generation.

**Repository Structure**

_Helper Files/_ – Includes the ODE equations and parameter definitions necessary for running the model simulations.

_1DSA/_ – Contains scripts for one-dimensional sensitivity analysis (1DSA) to evaluate parameter influence on FcR activation.

_Boosting/_ – Simulates increases in IgG1 or IgG3 levels to assess their impact on FcR activation and immune responses.

_Surfaces/_ – Contains scripts for generating all surface plots presented in the manuscript, depicting steady-state complex formations.

_Global Sensitivity/_ – Implements global sensitivity analysis (GSA), used in the supplementary materials to identify key model parameters.


**Dependencies**

MATLAB (R2022b used for all simulations)

Parallel Computing Toolbox

Mac OS 14
