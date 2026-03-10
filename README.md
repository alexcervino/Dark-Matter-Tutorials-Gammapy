# 🌌 Dark Matter Modeling and Simulations with Gammapy

Welcome to this repository of tutorials dedicated to modeling and simulating Dark Matter (DM) scenarios. 

This project provides a step-by-step workflow to understand the physical signatures of Weakly Interacting Massive Particles (WIMPs) and how to simulate their observation using the [Gammapy](https://gammapy.org/) framework and **CTA (Cherenkov Telescope Array)** Instrument Response Functions (IRFs).

---

## 📂 Repository Structure: The Tutorials

This repository is structured as a series of educational notebooks, taking you from theoretical spatial/spectral modeling to full observational simulations.

### 📓 Notebook 1: Basic Modeling of Dark Matter [tutorial_basics.ipynb] 
This notebook provides the foundational methods for theoretical Dark Matter modeling. It introduces the basic concepts of defining the spatial distribution and spectral signatures of a dark matter source, guiding you through the necessary steps to compute the expected, pristine gamma-ray flux arriving at Earth independent of any specific telescope.

### 📓 Notebook 2: Observational Simulation with Gammapy  [tutorial_observation_simulation.ipynb]
This notebook focuses on the observational side, demonstrating how to translate a theoretical physical model into realistic mock data. It introduces the use of Instrument Response Functions (IRFs) to simulate how a gamma-ray telescope actually observes the source, allowing you to generate synthetic datasets and apply statistical tools to evaluate signal detectability.

### 📓 Folder: cross_section_limits
In this folder you cna find the instructions to compute the annihilation or decay limits or your DM case.
---

## ⚙️ Installation and Setup

### Prerequisites
To run these notebooks, you will need a standard scientific Python environment with the following core libraries:
* `gammapy`
* `astropy`
* `scipy`
* `numpy`
* `matplotlib`

### Environmental Variables
The theoretical spectra generation relies on the PPPC4DMID or Cosmixs tables. You must download the Gammapy datasets and set the following environment variable on your machine pointing to that folder so Gammapy can find the `AtProduction_gammas.dat` file:

```bash
export GAMMAPY_DATA="/path/to/your/gammapy-data/"