# SPA3-Case-Study
 Model and simulations for McDonald et al., submitted to ICES Journal of Marine Science
 
 The files contained in this repository are the ones used for the paper titled *Explicit incorporation of spatial variability in a biomass dynamics assessment model*, submitted to the ICES Journal of Marine Science. 
 
 The model files are SPA3_SEBDAM.cpp and SPA3_TLM.cpp. They contain the SEBDAM model, as described in the paper mentionned above, and the TLM model as described in the paper titled *Incorporating intra-annual variability in fisheries abundance data to better capture population dynamics*, submitted to Fisheries Research. 
 
 The RData files contain what is necessary for fitting both TLM and SEBDAM for the SPA 3 case study. spa3_maps.RData contains the sf objects necessary for plotting the output from SEBDAM. spa3_unid_data.RData contains all the survey data formatted for fitting SEBDAM, while spa3_unid_TLM_data.RData contains the same data formatted for TLM. The only data that are not included are the landings used in the paper, as these are subject to privacy laws. To simulate the landings contained in these files, both TLM and SEBDAM were run without any landings. For both models, the total landings in a given year were simulated from a lognormal distribution with the mean set to 10% of the biomass in a given year and a variance of 0.2 on the log scale. For SEBDAM, the total landings were then distributed in each knot proportionally to their biomass density as predicted by the model fit without landings.
 
 The SPA3 Case Study.R file is the script used to fit both TLM and SEBDAM to the SPA 3 data and the code required for the plots related to their model fit. The simulations SEBDAM.R file contains the code for the simulations in the main body of the manuscript.
