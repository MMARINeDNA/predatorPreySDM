# predatorPreySDM
This repository contains code to examine potential prey species that co-occur with marine mammals, and quantify the effectiveness of potential prey proportional abundance in 3D spatial distribution models for marine mammals. A living version of the manuscript for the project can be found here (coming soon!).

## Goals and Hypotheses

This project aims to create 3D SDMs for marine mammals using a GAM framework and eDNA input data, with a focus on answering the following questions:

- RQ1: Do 3D eDNA SDMs arrive at the same mechanistic conclusions as 2D visual SDMs?
-   H1: Environmental predictors (e.g. SST, SSH, Chl-a) included in the best eDNA SDM for each species agree with environmental predictors included in published visual SDMs for the same species.
-   H2: Environmental predictors differ between SDMs built using the two datasets; this may correspond with the inclusion of a third spatial dimension (depth) in eDNA data.
- RQ2: Is 3D distribution of potential prey items correlated with 3D distribution of target marine mammal species?
- RQ3: Does the inclusion of potential prey items in 3D SDMs of marine mammals improve model performance?

## Approach

To analyze comparable datasets across multiple marine mammal species and ecosystems, we use eDNA detections from 3 techncal replicates of MarVer1 targeting all vertebrates in a biological sample. 

- RQ1: To address this question, we start with the 3D spatial model built in Van Cise et al. [in prep](https://mmarinedna.github.io/zDistribution/) and quantify improvement in model performance with the inclusion of biologically appropriate environmental predictors. We select canditate environmental predictors using previously published SDMs for each species, and download those datasets from MARSPEC and bioOracle using sdmpredictors in R. We test model performance using AIC, deviance explained, AUC, TSS, and a visual comparison to true detections.
- RQ2: To address this question, we use community analysis of variation (e.g., ANCOMBC package in R) to find vertebrate potential prey species that covary in space with our target marine mammal species.
- RQ3: Finaly, we add proportional eDNA data of vertebrate potential prey species to the best fit models from RQ1 and re-optimize model parameters using AIC, then quantify the improvement in model performance with the inclusion of potential prey items. We test model performance using AIC, deviance explained, AUC, TSS, and a visual comparison to true detections.
