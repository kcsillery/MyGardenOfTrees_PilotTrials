# MyGardenOfTrees_PilotTrials
**Data analysis scripts for the MyGardenOfTrees pilot trials: climate chamber and field trials**

The manuscript from this work is available on bioRxiv:

Cite it as:


The **data** folder contains raw and processed data files:

Abies_MLE.txt and Fagus_MLE.txt: maximum likelihood parameter estimates from the Markov model

climate_chamber_pheno_stages_seeds.csv: climate chamber germination data for survival analysis

dat_surv.csv: climate chamber germination data for plotting

dat2022.RData: field trial (micro-garden) data for mixed-effects models for germination; the dataset also contains the environmental variables

dat.RData: field trial (micro-garden) data for 3-year survival analysis

seed_dat_cols.csv and seed_provenance_details.csv: seed traits and provenance information, color code used

Table_S1_environmental_variables_LASSO.xlsx: the list of environmental variables considered for variable selection


The **scripts** folder contains the R scripts for analysis and visualization:

00_plot.r: functions for visualization of the Markov model outcome

0_CC_surv_analysis.r: raw data visualization, germination metrics, and survival analysis from the climate chamber trial data (Figure 2 of the manuscript)

1_MG_figure_rawdata.r: raw data visualization of the micro-garden data

2_MG_lasso_varselect.r: selection of the environmental variables for the mixed-effects models using LASSO regression

3_MG_mixedmodel_env.r: mixed-effects models of the germination data from the micro-gardens (**Asreml-R licence required**) (Figure 4 of the manuscript)

4_MG_plot_markovmodel.r: visualize the germination and developmental speed parameters estimated from the climate chamber and micro-garden trials (Figure 3 of the manuscript)

5_MG_3years_survival.r: 3-year survival analysis from the micro-gardens (Figure 5 of the manuscript)


The **docs** folder contains supplementary documents

Protocol_MyGardenOfTrees_Trials_2021-2026_EN.pdf: protocol for establishing the micro-gardens

XLSForms for performing the observations
