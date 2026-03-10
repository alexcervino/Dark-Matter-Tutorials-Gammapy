# Dark Matter limits calculator
This folder contains the core Python scripts used to derive upper limits on the Dark Matter (DM) properties: the **thermally averaged annihilation cross-section** ($\langle \sigma v \rangle$) in `limits_annihilation.py`, and the **inverse decay lifetime** ($1/\tau$) or **decay width** ($\Gamma$) in `limits_decay.py`.

The methodology, scripts, and mock datasets (i.e., JFactors and DFactors computed with clumpy) used to generate these results are directly derived from the official [CTA-Padova KSP dSph material repository](https://gitlab.cta-observatory.org/cta-consortium/aswg/sandbox/cta-padova/ksp_dsph_material/-/tree/main?ref_type=heads). This repository collects the core analytical framework utilized in the following CTA Consortium publication:

>  **Paper:** Abe, K., et al. (2025). Prospects for dark matter observations in dwarf spheroidal galaxies with the Cherenkov Telescope Array Observatory. Monthly Notices of the Royal Astronomical Society, 544(3), 2946–2986. **DOI:** [10.1093/mnras/staf1798](https://doi.org/10.1093/mnras/staf1798)

The tool supports both single-source and **stacked joint-likelihood analyses**, incorporating the astrophysical J-factor as a Gaussian nuisance parameter to account for uncertainties.


## 🚀 Usage

The workflow of the limits obtainment must be:

### 1. Obtain the limits
You can use limits_annihilation.py or limits_decay.py to calculate the limits. Both scripts share the same execution logic and command-line interface. Run them as follows:

**Important note:** For running these scripts you need to use the spectra.py file in this folder, since it is a duplicate of the spectra.py of Gammapy but with some modification made by the authors of the aforementioned KSP paper. Changing the spectra.py file of Gammapy in your environment with the one in the folder would be enough.

Run the script via the command line using the following arguments:

```bash
python limits_annihilation.py -s RetII DraI -c tau -p Einasto -t 100 -o tau_limits.csv
```
## ⚙️ Command-Line Arguments

The script is executed via the command line and offers the following parameters to customize your analysis:

| Short Flag | Long Flag | Description | Default Value |
| :---: | :--- | :--- | :--- |
| `-s` | `--source` | Target source or list of sources to analyze (space-separated). Example: `RetII DraI`. | `['RetII']` |
| `-c` | `--channel` | Dark Matter annihilation channel (e.g., `b`, `tau`, `mu`, `W`). | `b` |
| `-m` | `--mass` | List of Dark Matter masses to simulate, in GeV (space-separated). | `[100]` |
| `-p` | `--profile` | Dark Matter density profile to use (e.g., `Einasto`, `Burkert`). | `Einasto` |
| `-t` | `--time` | Simulated observation time in hours. | `100` |
| `-n` | `--nruns` | Number of independent mock datasets to generate per mass and channel. | `1` |
| `-b` | `--basepath` | Base directory path containing the required input materials. | `.` (Current directory) |
| `-o` | `--output` | Filename and path for the resulting CSV output. | `limits_output.csv` |
```
### 2. Quantiles
Once you have the limits computed, you can obtain the quantiles neccesary to plot the limits later with the notebook ``obtain_quantiles.ipynb``.

### 3. Visualization
In the notebook ``limits_cross_section_plot.ipynb`` you can plot all your limits results with the known Brazil plot pattern.
