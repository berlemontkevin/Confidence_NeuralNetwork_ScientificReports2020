# Confidence_NeuralNetwork_ScientificReports2020

Numerical simulations codes for the paper "Nonlinear neural network dynamics accounts for human confidence in a sequence of perceptual decisions" in Scientific Reports (2020) : https://www.nature.com/articles/s41598-020-63582-8

## Contents

### raw_experimental_data folder
This folder contains the raw data of all participants
- Manip 1: Pure Block
- Manip 2: Feedback block
- Manip 3: Confidence block

### theoretical model folder
This folder contains the code to reproduce the main figures of the articles. It is written in Julia
- **decision_making.jl**: This file implements the decision-making model (mean-field attractor neural network).
- **Fig_fitDM.jl**: This script reproduces the figure that shows the fit between the model and experimental data
- **figConfidence.jl**: This script reproduces the figure showing the relationship between reaction times, accuracy and confidence
- **genManip3.jl**: Generate the data corresponding to a simulation of Manip 3 conditions
- **variousfunctions.jl**: Definition of several functions necessary for the analysis

### R_scripts
This folder contains the R script that is necessary to reproduce the LMM. The script is the file *scriptLMMDATA.Rmd*

*Please contact K. Berlemont (kevin.berlemont@gmail.com) if you have any remark or questions about these files.*
