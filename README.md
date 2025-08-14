# OEHHA PFAS PK Modelling Application

Interactive R/Shiny app for simulating and visualizing pharmacokinetics (PK) of PFAS across species, model types, and exposure scenarios. Built to support transparent parameterization, quick what-ifs, and reproducible downloads for analysis and review.

> Contact: **Scott.Coffin@oehha.ca.gov**; **NTES@oehha.ca.gov**

---

## Table of contents

- [Objectives](#objectives)
- [Key Features](#key-features)
- [Models & Sources](#models--sources)
- [Quick Start](#quick-start)
- [Using the App](#using-the-app)
- [Data & Files](#data--files)
- [Deployment Notes (shinyapps.io)](#deployment-notes-shinyappsio)
- [How It Works (under the hood)](#how-it-works-under-the-hood)
- [Known Limitations / Roadmap](#known-limitations--roadmap)
- [Repository Structure](#repository-structure)
- [Contributing](#contributing)
- [Citation & Acknowledgements](#citation--acknowledgements)
- [License](#license)
- [FAQ](#faq)

---

## Objectives

- Provide a unified, point-and-click interface for PFAS PK simulations (PBPK + simple TK).
- Make parameters fully visible and editable, with provenance and units.
- Compare modeled vs. measured concentrations at user-specified sampling times.
- Offer reliable exports (time series, summary statistics, interactive plots).
- Encourage transparent, reproducible, literature-grounded modeling.

---

## Key Features

- **Experiment table (editable):** define PFAS, species, sex, model type, dose, interval, exposure duration, and sampling time; upload CSV/XLSX (template included).
- **Model availability & nudges:** see which model types are supported for each PFAS/species/sex; app suggests a more-preferred **supported** model, but never changes your inputs.
- **Parameter editors:**
  - **Simple TK models:** single-compartment, two-compartment, biphasic; inline edits; auto-calc `t½` or `Vd` when derivable.
  - **MassTransferPBPK (Fischer 2025):** parameter table auto-loads PFAS rows in use (Mouse/Male), editable inline; edits are honored at runtime.
  - **PBPK (PFOS, Chou & Lin 2019):** species/sex-specific parameter grids appear when PBPK rows are present.
- **Results & downloads:**
  - Concentration-time curves (multi-model, multi-compartment).
  - Modeled vs. measured scatter (log–log with 1:1 and ±3× guide lines).
  - Matched table of modeled/observed concentrations at your sampling hour.
  - Summary stats: **Cmax**, **AUC** (mg·hr/L), **C<sub>TWA</sub>**.
  - One-click exports (CSV/XLSX/HTML widgets).
- **TK data explorer:** full curated toxicokinetic dataset with sources, units, and filters.
- **Loading indicators:** spinners for heavy renders and a progress bar during simulations.
- **Robust text handling:** normalizes PFAS labels (e.g., `pfhxa` → `PFHxA`) to avoid mismatches.

> **Model preference (nudges only):**  
> `MassTransferPBPK` → `PBPK` → `two-compartment` → `biphasic` → `single-compartment`.

---

## Models & Sources

### PBPK — PFOS (Chou & Lin, 2019)
- **Model source**: [Chou & Lin 2019](https://linkinghub.elsevier.com/retrieve/pii/S016041201930203X).
- **Species:** Mouse, Rat, Monkey, Human (**Male** in this app’s optimized sets).
- **Route:** Oral (PFOS model is oral-only).
- **Structure:** 4 organ compartments (Plasma, Liver, Kidney, Rest) with Bayesian-optimized parameters.
- **Implementation:** `mrgsolve` models in `models/*.RDS`; parameters transformed out of log10 space at runtime.

### MassTransferPBPK
- **Model source:** [(Fischer et al. 2025)](https://doi.org/10.1021/acs.est.5c05473?urlappend=%3Fref%3DPDF&jav=VoR&rel=cite-as)
- **Scope in app:** **Mouse / Male**.
- **Mechanisms:** permeability–surface (passive) + active transport clearances; renal filtration/reabsorption; oral dosing via gut.
- **Compartments exposed:** Blood/Plasma, Liver, Kidneys, Gut, Rest (e.g., `C_plasma`, `C_blood`, …).
- **Parameters:** `Additional files/Datasets/Fischer/pfas_parameters.csv` (editable in-app; edits propagate to the simulation function).

### Simple TK models
- **Single-compartment:** requires `Volume_of_Distribution_L_per_kg (Vd)`, `Half_Life_hr`.
- **Two-compartment:** requires `Vd`, `Half_Life_hr`, `Absorption_Coefficient_unitless (k_abs)`.
- **Biphasic:** requires `VD2`, `VDc`, `k_alpha`, `k_beta`, `k_abs` (for PFAS with α/β dynamics, e.g., GenX, PFOA, PFBS, PFOS).
- **Parameter sources:** curated OEHHA dataset (`tk_params.rds`, `tk_df.rds`), with authors/year/PubMed/DOI displayed in the UI.

---

## Quick Start

### Requirements
- R ≥ 4.2 (recommended)
- Platform packages: see install block below

### Install & run locally
```r
install.packages(c(
  "shiny","shinydashboard","plotly","DT","ggplot2","mrgsolve","dplyr","purrr",
  "rhandsontable","readr","readxl","tidyverse","cols4all","shinyjs","shinythemes",
  "shinycssloaders","openxlsx"
))
# Optional but recommended for reproducibility:
# install.packages("renv"); renv::init()

shiny::runApp("app.R")
```
## Using the App

1. **Experiment Inputs**
   - Add/edit rows for PFAS, Species, Sex, Model_Type, Dose (mg/kg), Interval (h), Exposure (days), Sampling hour, and optional observed Serum concentration (mg/L).
   - Or **Upload** CSV/XLSX (template is downloadable from the app).
   - Use **Model Availability** to confirm what’s supported and which model is preferred.

2. **Model Parameters**
   - **Simple TK:** edit directly in the table; the app computes `t½` or `Vd` when derivable.
   - **MassTransferPBPK:** edit per-PFAS parameters (Mouse/Male only); the run uses your edited table.
   - **PBPK (PFOS):** species/sex parameter grids appear for PBPK rows.
   - Watch the **validation banner**—it lists what’s missing (if anything) before you run.

3. **Run Simulation**
   - Click **Run Simulation**. A progress bar shows run status; results populate the **Results** tab.

4. **Explore Results**
   - **Time Series:** filter by Species / PFAS / Model / Compartment; hover for details; download data.
   - **Modeled vs. Measured:** points at your sampling times (±3× guide bands); download an interactive HTML widget.
   - **Tables & Summary:** matched modeled/observed values and summary metrics (Cmax, AUC, C<sub>TWA</sub>) with CSV/Excel export.

5. **Reset**
   - Restores the app to its default demo state.

---

## Data & Files

- `Additional files/Datasets/general/tk_params.rds` — canonical TK parameters by PFAS/species/sex.
- `Additional files/Datasets/general/tk_df.rds` — full TK dataset with sources (for exploration and tooltips).
- `Additional files/Datasets/general/mismatch_data_for_shiny.rds` — long-form TK rows powering the heatmap.
- `Additional files/Datasets/general/validation_data.xlsx` — validation table available for download.
- `Additional files/Datasets/Fischer/pfas_parameters.csv` — MassTransferPBPK parameter catalog.

> PFAS labels are normalized internally (trim, uppercase, common synonyms) to keep filters and joins consistent.

---

## Deployment Notes (shinyapps.io)

**Suggested (single instance) starting point:**
- **Max Worker Processes:** 1–2  
- **Max Connections:** 25–50  
- **Worker Load Factor:** 80–90%  
- **Idle Timeout:** 60 s  
- **Read Timeout:** 3600 s  
- **Package Cache:** Enabled

**Performance tips:**
- Avoid high worker counts unless you have heavy concurrent traffic; each worker has a limited memory budget.
- Prefer a single instance with a couple of workers; the app caches expensive joins and normalizations.
- Plotly outputs use spinners; simulation shows a progress bar to keep the UI responsive.

---

## How It Works (under the hood)

- **Simulation engine**
  - **PBPK (PFOS):** `mrgsolve` models per species (params in log10 space → transformed at run); event dosing with hourly integration.
  - **MassTransferPBPK (Mouse/Male):** `Fischer_PBPK_mouse()` reads the active parameter CSV (including in-app edits) and returns compartment concentrations; the app pivots to long format and harmonizes units.
  - **Simple TK:** repeated-dose accumulation/decay; two-compartment absorption; biphasic α/β implementation with dynamic α-phase cutover.

- **Validation:** checks required parameters per model; flags unsupported PFAS/species/sex for PBPK/MassTransferPBPK.

- **Modeled vs. Measured matching:** robust nearest-time join by PFAS/species/sex/model and dosing keys, with tolerances and numeric key matching; joins at plasma/serum.

---

## Known Limitations / Roadmap

- PBPK run path wired to PFOS ([Chou & Lin 2019](https://linkinghub.elsevier.com/retrieve/pii/S016041201930203X)).
- MassTransferPBPK currently limited to **Mouse / Male** (per available function/params) [(Fischer et al. 2025)](https://doi.org/10.1021/acs.est.5c05473?urlappend=%3Fref%3DPDF&jav=VoR&rel=cite-as).
- Two-compartment and biphasic models require complete parameter sets (the app derives `t½` or `Vd` when possible).
- Future: broader PBPK coverage, additional routes, expand MassTransferPBPK to other species/sex, richer validation messaging and tooltips.

---

## Repository Structure

```bash
.
├─ app.R
├─ models/
│   ├─ ratPBPK.RDS
│   ├─ micePBPK.RDS
│   ├─ humanPBPK.RDS
│   └─ monkeyPBPK.RDS
├─ Additional files/
│   ├─ Code/Fischer 2025/Fischer_PBPK.R
│   └─ Datasets/
│       ├─ general/
│       │   ├─ tk_params.rds
│       │   ├─ tk_df.rds
│       │   ├─ mismatch_data_for_shiny.rds
│       │   └─ validation_data.xlsx
│       └─ Fischer/
│           └─ pfas_parameters.csv
├─ www/
│   ├─ main_logo_drop.png
│   ├─ fig1.PNG
│   ├─ PFOA_PBPK.jpeg
│   └─ acsest5c05473.JPG
└─ deployment_history.txt
```

---

## Contributing

Issues and PRs are welcome!

- Use feature branches and submit focused PRs.
- Describe the change and include before/after screenshots for UI tweaks.
- For model changes, cite sources and include validation notes or tests.

---

## Citation & Acknowledgements

If you use this app in your work, please cite:

- [**Chou & Lin (2019)**](https://linkinghub.elsevier.com/retrieve/pii/S016041201930203X) — Bayesian PBPK for PFOS across mouse, rat, monkey, human.  
- [**Fischer et al. (2025)**](https://doi.org/10.1021/acs.est.5c05473?urlappend=%3Fref%3DPDF&jav=VoR&rel=cite-as) — Permeability/active-transport PBPK (“MassTransferPBPK”).  
- OEHHA curated TK datasets (literature sources are exposed in-app).

Thanks to the maintainers of **mrgsolve**, **shiny**, **plotly**, and the R community.

---

## License

MIT License

Copyright (c) 2025 Office of Environmental Health Hazard Assessment

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

---

## FAQ

**Q: Why is PBPK skipped for non-PFOS PFAS?**  
A: The PBPK run path is currently wired to the PFOS model set. Use a simple TK model or MassTransferPBPK (where supported) for other PFAS.

**Q: Why is MassTransferPBPK flagged as unsupported?**  
A: It requires **Species = Mouse**, **Sex = Male**, and a PFAS present in the Fischer parameter file. Add an experiment row matching those constraints or extend the parameter file/model.

**Q: Can I reproduce an analysis?**  
A: Yes—export the time series and summary tables; include your edited parameter downloads. Consider using `renv` to lock package versions for exact reproducibility.

