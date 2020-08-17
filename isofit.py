'''
isofit.py

Cameren Swiggum (swiggum2@wisc.edu)

Classes to handle the fitting of isochrones
to a stellar population
'''
import os 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import scipy.interpolate as interp

from Cluster import cluster
from Isochrone import isochrone

np.seterr(divide='ignore')

class grid:
    '''
    Class to contain functions and attributes for a grid of isochrones

    Attributes:
        filepath -> str: path to isochrone grid file (downloaded from http://stev.oapd.inaf.it/cgi-bin/cmd_3.3)
        grid_list -> list: List of Isochrone objects
        grid_df -> pd.DataFrame: DataFrame representation of isochrone grid
    '''
    def __init__(self, filepath):
        '''
        init method that parses the isochrone grid file to define list of isochrones
        '''
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

        if len(set(grid_metallicity)) == 1: # Checks to see if metallicity values are varrying
            grouped_isochrones = self.grid_df.groupby(['age'])
            self.grid_list = populate_isochrone_list(grouped_isochrones) 

        else:
            grouped_isochrones = self.grid_df.groupby(['metallicity', 'age'])
            self.grid_list = populate_isochrone_list(grouped_isochrones) 
          

    @classmethod
    def large_grid(cls, path_to_grid_directory):
        '''
        Takes into account a grid of Isochrones with varying metallicity
        and extinction (Av) values. Directory should be formated as multiple
        subdirectories representing different extinction values and PARSEC files
        should contain metallicity (along with age) variations, themselves.
        '''
        directory = os.fsencode(path_to_grid_directory)

        # Each index of isochrone_large_grid will correspond to an isochrone grid with a single Av value and varying metallicities.
        isochrone_large_grid = []
        for parsec_file in os.listdir(directory):

            parsec_filename = os.fsdecode(parsec_file)
            if parsec_filename.endswith('.txt'):

                isochrone_grid = cls(
                        filepath = os.path.join(path_to_grid_directory, parsec_filename)
                        ) 

                isochrone_large_grid.append(isochrone_grid)

        return isochrone_large_grid
        
            
    def __str__(self):
        '''
        str representation of grid
        '''
        return 'Grid of {} PARSEC-COLIBRI isochrones with median log age of {}'.format(len(self.grid_list),np.median(self.grid_df['age']))
    
    def plot_grid(self, stellar_pop = None):
        '''
        Plots grid of isochrones and can also plot Cluster data
        @param stellar_pop -> cluster : Cluster object to be plotted with isochrones
        '''
        if stellar_pop != None:
            fig, ax = stellar_pop.plot_cluster_cmd(show=False)
        else:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        for isochrone in self.grid_list:
            ax.set_title('[M/H] = {}'.format(isochrone.metallicity))
            ax.plot(isochrone.color, isochrone.abs_mag, c='black', alpha = .6, linewidth=.2)
        plt.show()




class fitter:
    '''
    Fits isochrone to stellar population

    Attributes:
        stellar_pop -> cluster: data to be fitted to
        x -> np.ndarray: cluster's color
        y -> np.ndarray: cluster's absolute mag
        isochrones -> grid: grid of isochrones
    '''

    def __init__(self, stellar_pop : cluster, isochrones_grid : list):
        '''
        init method for fitter
        @param stellar_pop -> cluster, stellar grouping data cluster object
        @param isochrones -> list, list of isochrone objects
        '''

        self.stellar_pop = stellar_pop
        stellar_pop_photometry = self.stellar_pop.photometry
        self.x = stellar_pop_photometry['color']
        self.y = stellar_pop_photometry['abs_mag']
        self.y_error = stellar_pop_photometry['abs_mag_err']

        self.isochrones_grid = isochrones_grid
        

    def interpolate(self,isochrone) -> np.ndarray :
        '''
        Interpolates the parsec-colibri models
        so that their x-axis matches with the 
        data to calculate the residuals
        
        @param isochrone
        @return y_iso_interp -> np.ndarray
        '''

        x_iso = isochrone.color
        y_iso = isochrone.abs_mag

        # Limit the isochrone interpolation to the ranges of the cluster data
        x_iso = x_iso[np.where((x_iso >= min(self.x)) & (x_iso <= max(self.x)))]
        y_iso = y_iso[np.where((x_iso >= min(self.x)) & (x_iso <= max(self.x)))]
        x_iso = x_iso[np.where(y_iso >= min(self.y))]
        y_iso = y_iso[np.where(y_iso >= min(self.y))]

        f_y = interp.interp1d(x_iso, y_iso, fill_value = 'extrapolate')
        y_iso_interp = f_y(self.x)
        return(y_iso_interp)
   

    def calc_chi_square(self):
        '''
        Calculates the chi^2 value for each isochrone relative to data and 
        its errors.
        '''

        n_redder = []
        n_bluer = []
        
        for x in self.x:
            n_redder.append(len(np.where(self.x > x)[0]))
            n_bluer.append(len(np.where(self.x < x)[0]))

        n_redder = np.array(n_redder)
        n_bluer = np.array(n_bluer)
        weights = np.sqrt(n_redder/(n_bluer + 1))

        for isochrones in self.isochrones_grid:
            for isochrone in isochrones.grid_list:

                y_iso = self.interpolate(isochrone)
                chi2 = np.sum(((self.y - y_iso)/(self.y_error))**2)
                likelihood = -(1/2)*np.sum((chi2)**weights)
                isochrone.likelihood = likelihood
                isochrone.error = chi2


                
#                error.append(isochrone.error)
#        error = np.array(error)
#        error = error[np.isfinite(error)]

        
    def calc_square_error(self):
        '''
        Calculates error between model
        and data assigned it to the __error 
        attribute of the isochrone
        '''
        
        for isochrones in self.isochrones_grid:
            for isochrone in isochrones.grid_list:
                y_iso = self.interpolate(isochrone)
                isochrone.error = np.sum((self.y - y_iso)**2)

    def min_error(self,type_fit):
        '''
        Finds isochrone with minimum
        square error and returns it
        Sets this isochrone's best_fit
        attribute to true
        '''
            
        if type_fit == 'least squares':
            self.calc_square_error()

        elif type_fit == 'chi2':
            self.calc_chi_square()
        
        for isochrones in self.isochrones_grid:

            for i, isochrone in enumerate(isochrones.grid_list):

                if i == 0:
                    isochrone_min = isochrone

                elif isochrone.likelihood < isochrone_min.error:
                    isochrone_min = isochrone


        isochrone_min.best_fit = True
        isochrone_best_fit = isochrone_min 
        cluster.age = isochrone_best_fit.age
        return isochrone_best_fit

    def return_best_fit(self, type_fit = 'chi2', plot = True, save = False, save_path = None):
        '''
        Plots the best fit isochrone on 
        top of the data
        '''
        
        isochrone = self.min_error(type_fit)
        if save == True:
            isochrone.output_isochrone(self.stellar_pop.name)
        if plot == True:
            fig, ax = self.stellar_pop.plot_cluster_cmd(show=False)
            ax.set_title('[M/H] = {}'.format(isochrone.metallicity))
            ax.plot(isochrone.color,
                    isochrone.abs_mag,
                    c='black', 
                    linewidth=1, 
                    label=r'$\tau = $' + str(round(10**isochrone.age/1e6,3)) + 'Myr')
            plt.legend()
            plt.show()
        return isochrone.age


# Helper functions
def populate_isochrone_list(grouped_isochrones):

    grid_list = []
    for name, ind_isochrone in grouped_isochrones:

        color = ind_isochrone['color'].values
        abs_mag = ind_isochrone['abs_mag'].values
        age = ind_isochrone['age'].iloc[0]
        metallicity = ind_isochrone['metallicity'].iloc[0]

        grid_list.append(isochrone(color=color, abs_mag=abs_mag, age=age, metallicity=metallicity))
    return grid_list
    
        
