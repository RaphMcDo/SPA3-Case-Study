# User Guide
 Model and simulations for McDonald et al., submitted to ICES Journal of Marine Science
 
 The files contained in this repository were used for the manuscript titled *Explicit incorporation of spatial variability in a biomass dynamics assessment model*, submitted to the ICES Journal of Marine Science. 
 
 The model files are **SPA3_SEBDAM.cpp** and **SPA3_TLM.cpp**. They contain the SEBDAM model, as described in the manuscript mentioned above, and the TLM model as described in the manuscript titled *Incorporating intra-annual variability in fisheries abundance data to better capture population dynamics*, submitted to Fisheries Research. 
 
### Case Study

 The RData files contain what is necessary for fitting both TLM and SEBDAM for the SPA 3 case study. The file **spa3_maps.RData** contains the sf objects necessary for plotting the output from SEBDAM. The file **spa3_unid_data.RData** contains all the survey data formatted for fitting SEBDAM, while **spa3_unid_TLM_data.RData** contains the same data formatted for TLM.  The **SPA3 Case Study.R** file is used to fit both TLM and SEBDAM to the SPA 3 data and includes the code required for the plots related to their model fit. 
 
The landings used in the manuscript are not included because they are subject to privacy laws. To generate realistic landings data both TLM and SEBDAM were initially run without any landings. For both models, the total landings in a given year were then simulated from a lognormal distribution using 10% of the modeled mean biomass in a given year and a variance of 0.2 on the log scale. For SEBDAM, the total landings were then distributed in each knot proportionally to the biomass density predicted by the model fit without landings.
 
### Simulation 
 
The file **Simulations SEBDAM.R** contains the code for the simulations used in the main body of the manuscript.
