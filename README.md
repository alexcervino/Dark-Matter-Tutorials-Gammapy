# 🌌 Dark Matter Modeling and Simulations with Gammapy

Welcome to this repository of tutorials dedicated to modeling and simulating Dark Matter (DM) scenarios. 

This project provides a step-by-step workflow to understand the physical signatures of Weakly Interacting Massive Particles (WIMPs) and how to simulate their observation using the [Gammapy](https://gammapy.org/) framework.

---

## 📂 Repository Structure: The Tutorials

This repository is structured as a series of educational notebooks, taking you from theoretical spatial/spectral modeling to full observational simulations.

### 📓 Notebook 1: Dark Matter Indirect Detection with Gammapy: Basics [tutorial_basics.ipynb] 
This notebook provides the foundational methods for theoretical Dark Matter modeling. It introduces the basic concepts of defining the spatial distribution and spectral signatures of a dark matter source, guiding you through the necessary steps to compute the expected gamma-ray flux .

### 📓 Notebook 2: Dark Matter Data Handling with Gammapy  [tutorial_observation_simulation.ipynb]
This notebook focuses on the observational side of a dark matter analysis: translating a theoretical model into realistic data. It introduces Instrument Response Functions (IRFs) and shows how to use them to generate synthetic observations, as well as how to build equivalent datasets from real observations, comparing both cases. It also introduces the 1D vs. 3D analysis approaches used in the next notebook and the Asimov and Monte carlo simulation scenarios.


### 📓 Notebook 3: Dark Matter indirect search analysis with Gammapy  [tutorial_complete_analysis.ipynb]
This notebook covers a full dark matter indirect detection analysis pipeline, using simulated observations. It includes testing for signal detection, deriving upper limits or confidence intervals depending on the outcome, scanning over mass to build an exclusion curve, and computing the expected sensitivity bands (1σ/2σ) — the "Brazilian plot" — for both the annihilation and decay scenarios.

---

## ⚙️ Installation and Setup

To run these notebooks, you will need Gammapy v2.1

### Environmental Variables
The theoretical spectra generation relies on the PPPC4DMID or Cosmixs tables. You must download the Gammapy datasets and set the following environment variable on your machine pointing to that folder so Gammapy can find the `AtProduction_gammas.dat` file:

```bash
export GAMMAPY_DATA="/path/to/your/gammapy-data/"
