import numpy as np
import pandas as pd
from isofit import grid, fitter
from Cluster import cluster

isochrone_filepath = '/Users/admin/Desktop/astro_research/orion/orion_populations/parsec_colibri_isochrones/output43.dat'
df = pd.read_csv('/Users/admin/Desktop/astro_research/orion/orion_populations/apogee_gaia_sample.csv')

my_grid = grid(isochrone_filepath)

lambda_orion_df = df.loc[df.label == 4]
lambda_orion = cluster(lambda_orion_df,'Lambda Orion')
lambda_orion.photometric_cuts()


my_fitter = fitter(lambda_orion, my_grid)
my_fitter.plot_best_fit()

