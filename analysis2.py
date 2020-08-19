import numpy as np
import pandas as pd
from isofit import grid, fitter
from Cluster import cluster
import matplotlib.pyplot as plt

cluster_output_directory = '/Users/cam/Desktop/astro_research/orion/orion_populations/age_analysis_groups/'
isochrone_filepath = '/Users/cam/Desktop/astro_research/orion/orion_populations/parsec_colibri_isochrones/second_analysis/'

def cluster_age_analysis(df, cluster_label : int, cluster_name : str, cost_function : str):
    
    df_group = df.loc[df['label'] == cluster_label]
    group_cluster = cluster(df_group, cluster_name, error_type = cost_function)
    #group_cluster.photometric_cuts()
#    group_cluster.output_cluster(directory = cluster_output_directory)

    iso_grid = grid.large_grid(isochrone_filepath)
    cluster_fitter = fitter(group_cluster, iso_grid)
    age = cluster_fitter.return_best_fit(type_fit = 'chi2', plot = True, save = True)

    #for g in iso_grid:
    #    g.plot_grid(stellar_pop = group_cluster)
    return age, df_group


df = pd.read_csv('/Users/cam/Desktop/astro_research/orion/apogee_gaia_sample4.csv')
df = df.loc[df['probability'] > .1]
#df = df.loc[df['phot_bp_mean_mag'] - df['phot_rp_mean_mag'] < 2.2]

#age_ngc1980 = cluster_age_analysis(df, 5, 'NGC 1980', 'astrometric')
#age_obpNear = cluster_age_analysis(df, 3, 'OBP-near', 'astrometric_photometric')

for i in range(1, 25):
    cluster_age_analysis(df, i, 'cluster', 'astrometric_photometric')
