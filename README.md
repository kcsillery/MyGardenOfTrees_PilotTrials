# 🌱 MyGardenOfTrees — Pilot Trials

**Data analysis scripts for the *MyGardenOfTrees* pilot trials**, including climate chamber experiments and field (micro-garden) trials.

🌐 **Project website:** https://www.mygardenoftrees.eu

---

## 📄 Manuscript

Katalin Csillery, Justine Charlet de Sauvage, Madleina Caduff, Johannes Alt, Marjorie Bison, Mert Celik, Nicole Ponta, Daniel Wegmann. 2026. Provenance and environment jointly shape early regeneration in fir and beech: Evidence from distributed participatory trials across Europe. Accepted in **New Phytologist**.

Link to the published version will follow. In the meantime, an earlier version is available here: https://doi.org/10.64898/2026.01.29.702314 

---

## 📁 Repository Structure


### XLSforms/
XLSForms and associated media files used for standardized data collection in the *MyGardenOfTrees* citizen science project. Dataset S1 is available from this repository.

- **`Dataset_S1.zip`**
  All forms and associated media for download
  
- **`first_year_form/`**
  Observation forms used in the first year of the participatory trials to record germination phenology and survival.

- **`second_third_year_form/`**  
  Observation forms used in the second and third years of the participatory trials to record spring budbreak phenology and survival.

---

### `data/`
Raw and processed datasets used in the analyses.

- **`Abies_MLE.txt`**, **`Fagus_MLE.txt`**  
  Maximum likelihood parameter estimates from the Markov model

- **`climate_chamber_pheno_stages_seeds.csv`**  
  Climate chamber germination data for survival analysis

- **`dat_surv.csv`**  
  Climate chamber germination data for plotting

- **`dat2022.RData`**  
  Field trial (micro-garden) data for germination mixed-effects models  
  *(includes environmental variables)*

- **`dat.RData`**  
  Field trial (micro-garden) data for 3-year survival analysis

- **`seed_dat_cols.csv`**, **`seed_provenance_details.csv`**  
  Seed traits and provenance information (including color codes)

- **`Table_S1_environmental_variables_LASSO.xlsx`**  
  Environmental variables considered for variable selection (Table S2 in the accepted manuscript)
---

### `GDD/`
Growing-degree-days used in the Markov model for the climate chamber cycles and micro-gardens.

---

### `scripts/`
R scripts for statistical analysis and figure generation.

- **`00_plot.r`**  
  Visualization functions for Markov model outputs

- **`0_CC_surv_analysis.r`**  
  Climate chamber trial: raw data visualization, germination metrics, and survival analysis  
  *(Figure 2)*

- **`1_MG_figure_rawdata.r`**  
  Raw data visualization for micro-garden trials

- **`2_MG_lasso_varselect.r`**  
  Environmental variable selection using LASSO regression

- **`3_MG_mixedmodel_env.r`**  
  Mixed-effects models for micro-garden germination data  
  ⚠️ *Requires an **Asreml-R** license*  
  *(Figure 4)*

- **`4_MG_plot_markovmodel.r`**  
  Visualization of germination and developmental speed parameters from climate chamber and micro-garden trials  
  *(Figure 3)*

- **`5_MG_3years_survival.r`**  
  Three-year survival analysis from micro-garden trials  
  *(Figure 5)*
  
- **`7_MG_climate_transfer_distance.r`**  
  Germination rate as a function of geographic and environmental transfer distances between provenance and garden
  *(Figure S11)*

- **`8_MG_prov_specific_slopes.r`**  
  Testing provenance-specific slopes for environmental gradients 
  *(Non-significant. Wald-test results reported in the text)*

---
### `docs/`

Supporting Information for the accepted manuscript: Figs S1–S11, Tables S1–S5, and Methods S1–S3. Methods S1–S3 are also supplied to the journal as Supporting Information.

- **`Figures_S1-11_Table_S1-5.pdf`**  
Supporting Information for the accepted manuscript:
  - **Fig. S1:** Climate chamber trays for fir (*Abies*) and beech (*Fagus*) germination.
  - **Fig. S2:** Phenological stages used to score germination and early seedling development.
  - **Fig. S3:** Cumulative germination as a function of growing degree-days in the field experiments.
  - **Fig. S4:** Fir (*Abies*) phenological development in climate chambers and field experiments.
  - **Fig. S5:** Beech (*Fagus*) phenological development in climate chambers and field experiments.
  - **Fig. S6:** Diagnostic plots for mixed-effects models of germination.
  - **Fig. S7:** Climate chamber germination trajectories across provenances.
  - **Fig. S8:** Effect of experimental cycle length on germination metrics.
  - **Fig. S9:** Correlations among germination metrics, geography, and seed traits.
  - **Fig. S10:** LASSO-selected environmental predictors ranked by their standardized regression coefficients for fir and beech.
  - **Fig. S11:** Germination rate as a function of geographic and environmental transfer distances.
  - **Table S1:** Information about the seed stands and origin certificates.
  - **Table S2:** Environmental covariates used for variable selection in the LASSO regression.
  - **Table S3:** Germination metrics across species, provenances, treatments, and experimental cycles.
  - **Table S4:** Wald tests for fixed effects retained in the final mixed-effects germination models.
  - **Table S5:** Generalized additive models (GAMs) of germination rate as a function of geographic and environmental transfer distances.

- **`Methods_S2.pdf`**  
  Climate and soil data and variable selection using LASSO regression.
  
- **`Methods_S3.pdf`**  
  Description of a hidden Markov model of seed germination, phenological development, and mortality, including the inference and importance sampling scheme.
---

## 🧰 Requirements

- **R version:** 4.3.3 (2024-02-29) — *“Angel Food Cake”*
- Required R packages (see individual scripts)
- **Asreml-R** license for mixed-effects model analyses

---

## ⚙️ Model Implementation

The model and inference scheme described in the manuscript were implemented as a **command-line C++ program**, **`tree_growth`**, using the **`stattools`** library.

- **Source code & user manual:**  
  https://bitbucket.org/wegmannlab/tree_growth/

- **Reproducibility note:**  
  All estimations presented in the manuscript were performed using commit  
  **`a9ae485`**.
