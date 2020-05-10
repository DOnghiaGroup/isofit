'''
Cameren Swiggum (swiggum2@wisc.edu)

Classes to handle the fitting of isochrones
to a stellar population
'''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

from Cluster import cluster
from Isochrone import isochrone

class grid:
    '''
    Grid of isochrones.
    '''
    def __init__(self, filepath):
        self.filepath = filepath
        self.grid_list = []
        
        # indices [1,2,-1,-2,-3] correspond to: MH, logAge, G_RP, G_BP, Gmag
        isochrones_grid = np.genfromtxt(filepath)
        grid_colors = isochrones_grid[:,-2] - isochrones_grid[:,-1]
        grid_abs_mags = isochrones_grid[:,-3]
        grid_ages = isochrones_grid[:,2]
        grid_metallicity = isochrones_grid[:,1]

        # DataFrame representation of grid
        self.grid_df = pd.DataFrame({'color':grid_colors, 'abs_mag':grid_abs_mags, 'age':grid_ages, 'metallicity':grid_metallicity}) 

        # Populate list of Isochrone objects to represent grid
        grouped_isochrones = self.grid_df.groupby(['age'])
        for name, ind_isochrone in grouped_isochrones:
            self.grid_list.append(isochrone(ind_isochrone['color'].values, ind_isochrone['abs_mag'], ind_isochrone['age'].iloc[0], ind_isochrone['metallicity'].iloc[0]))

    def __str__(self):
        return 'Grid of {} PARSEC-COLIBRI isochrones with median log age of {}'.format(len(self.grid_list),np.median(self.grid_df['age']))


class fitter:
    '''
    Fits isochrone to stellar population
    '''

    def __init__(self, stellar_pop : cluster, isochrones : grid):

        self.stellar_pop = stellar_pop
        stellar_pop_photometry = self.stellar_pop.photometry
        self.x = stellar_pop_photometry['color']
        self.y = stellar_pop_photometry['abs_mag']

        self.isochrones = isochrones

        self.isochrone_best_fit
    
    def interpolate(self) -> np.ndarray :

        x_iso = isochrone.color
        y_iso = isochrone.abs_mag
        x_iso = x_iso[np.where((x_iso >= min(self.x)) & (x_iso <= max(self.x)))]
        y_iso = y_iso[np.where((x_iso >= min(self.x)) & (x_iso <= max(self.x)))]
        x_iso = x_iso[np.where(y_iso >= min(self.y))]
        y_iso = y_iso[np.where(y_iso >= min(self.y))]

        f_y = interp.interp1d(x_iso, y_iso, fill_value = 'extrapolate')
        y_iso_interp = f_y(x)
        return(y_iso_interp)
            
    def calc_square_error(self):
        
        for isochrone in self.isochrones:
            y_iso = interpolate(self)
            isochrone.error = np.sum((y - y_iso)**2)

    def min_square_error(self):
        
        calc_square_error(self)
        for i, isochrone in enumerate(self.isochrones):
            if i == 0:
                isochrone_min = isochrone
            elif isochrone.error <= isochrone_min.error:
                isochrone_min = isochrone
        isochrone_min.best_fit = True
        self.isochrone_best_fit = isochrone_min 
