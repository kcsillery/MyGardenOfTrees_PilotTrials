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
Supplementary documentation.

- **`Protocol_MyGardenOfTrees_Trials_2021-2026_EN.pdf`**  
  Protocol for establishing and monitoring micro-gardens

- **XLSForms**  
  Forms used for field and experimental observations

---

## 🧰 Requirements

- R (version used in the analyses recommended)
- Required R packages (see individual scripts)
- **Asreml-R** license for mixed-effects model analyses
