6/28/2022 Hsiang-Yu Yuan (sean.yuan@cityu.edu.hk)

This file explains how to 
i) produce the following main figures in the manuscript about the impacts of public health interventions and weather on SARS-CoV-2 Omicron outbreak in Hong Kong;
ii) reproduce MCMC results.

Figure 1C: run plotFigure1_input.m to plot daily changes in vaccination, mobility and seasonal factors
Figure 2AB: run plotFigure2_transmission_comparisons
Figure 2CD: run plotFigure2_transmission_dynamics
Figure 3: run plotFigure3_population_immunity
Figure 4: run plotFigure4_npi_impacts

To produce these figures, most of the files need to resample from the posterior distribution, which may takes few minutes depending on the sample size.
Data scource: Stored in HK_virus.mat. The mat file includes number of PCR detected cases (cases), number of RAT detected cases (cases_rat), daily BNT booster rate (dailybnt), daily CoronaVac booster rate (dailysinovac), raw daily mobility (dailymobility), 7-day moving averaged daily mobility (mobility_7d), daily mean temerature (temperature), and daily mean relative humidity (relhumid).

To reproduce MCMC results, run main_MCMC_m1_12_151.m. The output will be stored in the 'out' folder.

When you plot the model outputs, replace the default MCMC file name and its path in the following line (e.g. line 10 in plotFigure2_transmission_dynamics).
load('mcmc/fullmod_mcmc_output');
 