import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import gaussian_kde
import isofit

df = pd.read_csv('/Users/admin/Desktop/astro_research/orion/orion_populations/apogee_gaia_sample.csv')
isochrone_file = np.genfromtxt('/Users/admin/Desktop/astro_research/orion/orion_populations/parsec_colibri_isochrones/output43.dat')
isochrone_data = isofit.tabulate_isochrone_data(isochrone_file,[1,2,-1,-2,-3],['MH','logAge','G_RP','G_BP','Gmag'])
lambda_orion = df.loc[df.label == 4]
isofit.one_cluster(lambda_orion,isochrone_data)
