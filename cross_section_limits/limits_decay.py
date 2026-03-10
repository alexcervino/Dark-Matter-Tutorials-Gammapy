#!/usr/bin/env python
# coding: utf-8

"""
CTA Dark Matter Limits Calculator
This script simulates CTA observations of dwarf spheroidal galaxies and computes 
the 95% C.L. upper limits on the Dark Matter Decay cross-section.
It supports both single-source and stacked joint-likelihood analyses.
"""

import argparse
import time
import os
import logging
import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.interpolate import interp1d

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits

from gammapy.irf import load_irf_dict_from_file
from gammapy.maps import WcsGeom, MapAxis
from gammapy.modeling.models import (
    TemplateSpatialModel,
    SkyModel,
    Models,
    FoVBackgroundModel
)
from gammapy.makers import MapDatasetMaker, SafeMaskMaker
from gammapy.modeling import Fit
from gammapy.data import Observation
from gammapy.datasets import MapDataset, Datasets
from gammapy.astro.darkmatter import DarkMatterDecaySpectralModel, PrimaryFlux

log = logging.getLogger(__name__)

# ==========================================
# Argument Parsing
# ==========================================
def getArgs():
   parser = argparse.ArgumentParser(
      description="Calculates upper limits for DM Decay cross-section."
   )
   
   parser.add_argument('-s', '--source', default=['RetII'], nargs='+',
                     help='Name of the source(s) (e.g., RetII DraI)')
   parser.add_argument('-c', '--channel', default='b',
                     help='Dark matter channel (b, mu, tau, W, etc.)')
   parser.add_argument('-p', '--profile', default='Einasto',
                     help='Dark matter density profile (Einasto, Burkert...)')
   parser.add_argument('-t', '--time', default="100", type=int,
                     help='Observation time in hours')
   parser.add_argument('-n', '--nruns', default=1, type=int,
                     help='Number of runs per DM mass and channel')
   parser.add_argument('-o', '--output', default='limits_output.csv',
                     help='Output CSV filename')
   parser.add_argument('-b', '--basepath', default='.',
                        help='Base directory path containing the input materials (IRFs, Dfactors, etc.)')
   parser.add_argument('-m','--mass', default=[100] , nargs='+',
                  help='DM masses in GeV to process')

   return parser.parse_args()

# ==========================================
# Custom MapDataset Class
# ==========================================
class DMMapDataset(MapDataset):
    """Dark matter template dataset with nuisance parameters and likelihood."""

    def __init__(self, nuisance=None, **kwargs):
        super().__init__(**kwargs)
        self.nuisance = nuisance

    @property
    def nuisance(self):
        return self._nuisance

    @nuisance.setter
    def nuisance(self, nuisance):
        try:
            if nuisance:
                assert (
                    nuisance["d"].unit == self.models.parameters["dfactor"].unit
                    and nuisance["dobs"].unit == self.models.parameters["dfactor"].unit
                    and nuisance["sigmad"].unit == self.models.parameters["dfactor"].unit
                ), "Different units in J factor"
                if "width" not in nuisance:
                    nuisance["width"] = 5
            self._nuisance = nuisance
        except Exception as ex:
            log.error(ex)
            self._nuisance = None
            raise ValueError("Nuisance parameters cannot be set")

    def stat_sum(self):
        wstat = super().stat_sum()
        liketotal = wstat
        if self.nuisance:
            liketotal += self.jnuisance()
        return liketotal

    def jnuisance(self):
        lg_dfactor = np.log10(self.models.parameters["dfactor"].value)
        jobs = self.nuisance["dobs"].value
        lg_jobs = np.log10(jobs)
        sigma = self.nuisance["sigmad"].value
        
        ee = np.power(lg_dfactor - lg_jobs, 2) / (2 * sigma * sigma)
        norm = np.log(10) * jobs * np.sqrt(2 * np.pi) * sigma
        
        res = np.exp(-ee) / norm
        return -2 * np.log(res)

# ==========================================
# Helper Functions
# ==========================================
def dfactor_inter1d(path_j, name, profile):
    """Interpolates the J-factor based on the alpha value."""
    filename = os.path.join(path_j, profile, f'rvir_{name}_MCMC.dat.Dalphaint_cls.output')
    df_j = pd.read_csv(filename, sep=r'\s+', header=None, comment='#',
                       names=['alpha', 'dfactor', 'd1', 'd2', 'd3', 'd4', 'd5', 'd6', 'd7'])
    return interp1d(df_j.alpha.values, df_j.dfactor.values, kind="quadratic")

def setModel(name, mass, ch, d_factor, spatial_model, time_decay):
    """Builds the SkyModel combining spectral DM physics and spatial templates."""
    DarkMatterDecaySpectralModel.TIME_DECAY = time_decay
    
    flux_model_fit = DarkMatterDecaySpectralModel(mass=mass, channel=ch, dfactor=d_factor)
    model_fit = SkyModel(spatial_model=spatial_model, spectral_model=flux_model_fit, name=name)
    bkg_model = FoVBackgroundModel(dataset_name=name)

    models_fit = Models([model_fit, bkg_model])

    models_fit.parameters["norm"].frozen = False
    models_fit.parameters["tilt"].frozen = False
    models_fit.parameters["norm"].value = 1
    models_fit.parameters["tilt"].value = 0

    models_fit.parameters["tau"].frozen = False
    models_fit.parameters["tau"].value = 1
    models_fit.parameters["tau"].min = 1e-6
    models_fit.parameters["dfactor"].min = 1e16 

    return models_fit

# ==========================================
# Main Execution
# ==========================================
if __name__ == "__main__":
    args = getArgs()

    # Define Paths
    mainpath = args.basepath
    path_j = os.path.join(mainpath, "cross_section_limits", "clumpy", "Dfactors", "profile")
    path_irfs = os.path.join(mainpath, "cross_section_limits", "IRFs", "cta", "prod5-v0.1", "bcf")
    path_2d = os.path.join(mainpath, "cross_section_limits", "clumpy", "Dfactors", "2Dmaps")
    
    PrimaryFlux.table_filename = os.path.join('/Users/alexcervino/Desktop/DARKMATTER/gammapy-data/dark_matter_spectra', "AtProduction_gammas.dat")
    input_csv_path = os.path.join(mainpath, "cross_section_limits", "inputs_combined_analisys.csv")

    df_input = pd.read_csv(input_csv_path)
    df_input = df_input.loc[df_input['src'].isin(args.source)]

    dSphNames = df_input['src'].values
    dfac_sigma = df_input['dfac_sigma'].values
    offset = df_input['r_tidal'].values * u.Unit("deg")
    alpha = df_input['r_tidal'].values
    irf_list = df_input['irf'].values

    # Pre-load spatial models and J-factors
    dfac = []
    hdul = []
    spatial_model = []
    
    for i, src in enumerate(dSphNames):
        j_val = dfactor_inter1d(path_j, name=src, profile=args.profile)(alpha[i]) * u.Unit("GeV2 cm-5")
        dfac.append(j_val)
        
        dfactor_filename = os.path.join(path_2d, args.profile, f'decay_{src}2D_FOVdiameter10.0deg_nside1024-DFACTOR-Jtot-image.fits')
        hdul.append(fits.open(dfactor_filename))
        spatial_model.append(TemplateSpatialModel.read(dfactor_filename))

    GLON = hdul[0][0].header['PSI_0'] * u.Unit("deg")
    GLAT = hdul[0][0].header['THETA_0'] * u.Unit("deg")
    src_pos = SkyCoord(GLON, GLAT, frame="galactic")

    # Geometry & Constants
    time_decay = 1e30 * u.Unit("s")
    CL_VALUE = 4.00  # For a 95% CL (1 DOF)
    
    axis = MapAxis.from_edges(np.logspace(np.log10(0.03), np.log10(100), 31), unit="TeV", name="energy", interp="log")
    geom = WcsGeom.create(skydir=src_pos, binsz=0.02, width=(2,2), frame="galactic", axes=[axis])

    maker = MapDatasetMaker(selection=["exposure", "background", "psf", "edisp"])
    
    masses = np.array(args.mass, dtype=float) * u.GeV    
    channels = [args.channel]

    # Setup Output Table
    columns = ['profile', 'ch', 'mass', 'tau']
    for i, name in enumerate(dSphNames):
        columns.insert(i, f'src{i+1}')
    final_table = pd.DataFrame(columns=columns)

    # Main Loop
    for run in range(args.nruns):  
        dataset = []
        
        # 1. Generate Fake Datasets
        for index, name in enumerate(dSphNames):
            filename = os.path.join(path_irfs, irf_list[index], 'irf_file.fits.gz')
            irfs = load_irf_dict_from_file(filename)
            obs_ = Observation.create(pointing=src_pos, livetime=args.time * u.hour, irfs=irfs)
            maker_safe_mask_ = SafeMaskMaker(methods=["offset-max"], offset_max=offset[index])
            
            d_fake = maker.run(DMMapDataset.create(geom, name=f"dataset_{name}"), obs_)
            d_fake = maker_safe_mask_.run(d_fake, obs_)
            d_fake.fake(int(time.time()))
            dataset.append(d_fake)

        # 2. Iterate over Channels and Masses
        for ch in channels:
            for mass in masses:
                print(f"--- Processing Run {run+1}/{args.nruns} | Mass: {mass} | Channel: {ch} ---")

                for i, src in enumerate(dSphNames):
                    dataset[i].models = setModel(f'dataset_{src}', mass, ch, dfac[i], spatial_model[i], time_decay)
                    dataset[i].nuisance = dict(j=dfac[i], jobs=dfac[i], sigmad=dfac_sigma[i] * u.Unit("GeV2 cm-5"), sigmatau=0.01)

                # Link sv parameter for stacking
                for i in range(len(dSphNames) - 1):
                    dataset[i+1].models[0].spectral_model.tau = dataset[0].models[0].spectral_model.tau

                datasets = Datasets(dataset)
                fit = Fit()

                try:
                    # Fit 1: Background only
                    datasets.models.parameters["tau"].frozen = True
                    fit.run(datasets=datasets)
                    lg_MLE_null = datasets.stat_sum()

                    # Fit 2: Unfreeze SV
                    datasets.models.parameters["tau"].frozen = False
                    fit.run(datasets=datasets)

                    # Fit 3: Unfreeze J-factor (Joint Fit)
                    for i, src in enumerate(dSphNames):
                        datasets.models[f'dataset_{src}'].parameters['dfactor'].frozen = False
                        datasets.models[f'dataset_{src}-bkg'].parameters['norm'].frozen = True
                        datasets.models[f'dataset_{src}-bkg'].parameters['tilt'].frozen = True

                    result = fit.run(datasets=datasets)
                    if not result.success:
                        print("Fit failed, skipping...")
                        continue

                    # 3. Profile Likelihood Scan
                    scale_min = datasets.models.parameters["tau"].value
                    scale_error = datasets.models.parameters["tau"].error  
                    scale_max = scale_min + 10 * scale_error

                    if scale_max > 1e5 or scale_min < 1e-5:
                        scale_min, scale_max = 1e-5, 10000

                    # Find upper bound crossing TS threshold
                    counter = 0
                    while True:
                        datasets.models.parameters['tau'].scan_n_values = 1
                        datasets.models.parameters['tau'].scan_min = scale_max
                        datasets.models.parameters['tau'].scan_max = scale_max
                        
                        ts_diff = np.max(fit.stat_profile(reoptimize=True, datasets=datasets, parameter='tau')['stat_scan'] - lg_MLE_null)
                        
                        if ts_diff >= (CL_VALUE + 1) or counter >= 100 or scale_max >= 1e10:
                            break
                        
                        scale_max *= 2
                        counter += 1

                    # Precise Scan and Interpolation
                    datasets.models.parameters['tau'].scan_n_values = 100
                    datasets.models.parameters['tau'].scan_min = scale_min
                    datasets.models.parameters['tau'].scan_max = 2 * scale_max
                    
                    profile = fit.stat_profile(reoptimize=True, datasets=datasets, parameter='tau')
                    xvals = profile[f"dataset_{dSphNames[0]}.spectral.tau_scan"]
                    yvals = profile["stat_scan"] - lg_MLE_null - CL_VALUE

                    scale_found = brentq(interp1d(xvals, yvals, kind="quadratic"), scale_min * 1.01, scale_max, maxiter=100, rtol=1e-5)
                    tau = time_decay / scale_found
                    print(f"Calculated Limit for {mass}: {tau:.2e}")

                except Exception as ex:
                    print(f"Error during fitting: {ex}")
                    continue

                # Store Results
                tmp_row = {f'src{i+1}': name for i, name in enumerate(dSphNames)}
                tmp_row.update({'profile': args.profile, 'ch': ch, 'mass': mass.value, 'tau': tau})
                
                final_table = pd.concat([final_table, pd.DataFrame([tmp_row])], ignore_index=True)

    # Clean and Save
    final_table['mass'] = final_table['mass'].astype(int)
    final_table['tau'] = final_table['tau'].astype(float)
    final_table.to_csv(args.output, index=False)
    print(f"Results saved to {args.output}")