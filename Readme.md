# Technological Employment and Gender Inequalities

## Description

This project analyzes the impact of technological progress on female unemployment across five Western European countries: France, Spain, Italy, the Netherlands, and Belgium.
Using **spatial econometric models**, we explore how tech employment, science & technology workforce, regional GDP, higher education, and population density affect female unemployment, both directly and through spillover effects on neighboring regions.

---

## 📁 Project Structure

```
├── .git/                 → Git repository folder  
├── data/                 → Folder with raw and processed data and shapfiles 
├── data.csv              → Main dataset (cleaned and combined)  
├── functions.R           → Custom R functions  
├── notebook.Rmd          → R Markdown notebook (analysis & report)  
├── paper.pdf             → Full research paper  
├── Readme.md             → This file  
├── script.R             → Main R script (data prep, modeling, results)
```

---

## Requirements

Make sure you have R and the following libraries installed:

```R
library(sf)
library(plm)
library(splm)
library(sp)
library(spdep)
library(dplyr)
library(ggplot2)
library(lmtest)
library(modelsummary)
library(stargazer)
library(tidyr)
```

---

## Data

* **Source:** Eurostat, regional data (2012–2021)
* **Variables:**

  * Female unemployment rate
  * Tech employment share
  * HRST (Science & Technology workforce)
  * Regional GDP
  * Higher education share
  * Population density

---

## Methodology

* **Descriptive analysis** → Explore trends & regional disparities
* **Panel data models** → Fixed effects, random effects
* **Spatial econometric models** → SAR, SEM, SAC to account for spillover effects
* **Key finding:** Tech jobs reduce female unemployment, but increased STEM competition can worsen gender gaps

---

## 📈 Results

* Tech employment ↓ female unemployment (local + neighboring regions)
* STEM workforce ↑ female unemployment (competition effect)
* GDP growth → limited effect
* Higher education → clear positive effect
* Strong spatial dependencies confirmed

---

## How to Run

1. Clone the repository

   ```
   git clone <repository_link>
   ```
2. Open `script.R` or `notebook.Rmd` in RStudio
3. Run the code to reproduce the analysis and results

---

## 📚 Reference

Read the full paper → [paper.pdf](./paper.pdf)

---

## 🔗 Links

* 📄 **Post / Presentation:** [link](https://auradev.hashnode.dev/inegalites-de-genre-et-technologie-ce-que-nous-apprend-leconometrie-spatiale)
