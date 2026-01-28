# 🌱 MyGardenOfTrees — Pilot Trials

**Data analysis scripts for the *MyGardenOfTrees* pilot trials**, including climate chamber experiments and field (micro-garden) trials.

🌐 **Project website:** https://www.mygardenoftrees.eu

---

## 📄 Manuscript

The manuscript associated with this work is available on **bioRxiv**.

> *Citation information will be added here.*

---

## 📁 Repository Structure

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
  Environmental variables considered for variable selection

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

---

### `docs/`

Supporting Information files as referenced in the manuscript.

- **`Figures_S1-14_Table_S1.pdf`**  
  Supplementary figures and table, including:
  - Experimental setups for climate chamber and micro-garden trials
  - Phenological stages used to score germination and early development
  - Environmental predictor selection and regression coefficients
  - Germination trajectories in climate chambers and field experiments
  - Phenological development of fir (*Abies*) and beech (*Fagus*) seedlings
  - Model diagnostics  
  - **Table S1:** Germination metrics across species, provenances, treatments, and experimental cycles

- **`Methods_S1.pdf`**  
  Protocols for establishing and monitoring the *MyGardenOfTrees* climate chamber and micro-garden experiments.

- **`Methods_S2.pdf`**  
  Description of a hidden Markov model of seed germination, phenological development, and mortality, including the inference and importance sampling scheme.

- **`Dataset_S1.zip`**  
  XLSForms and associated media files used for standardized data collection in the *MyGardenOfTrees* citizen science project.

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
