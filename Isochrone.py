'''
Cameren Swiggum (swiggum2@wisc.edu)

Class to represent data from Parsec/Colibri 
Isochroens. Found at http://stev.oapd.inaf.it/cgi-bin/cmd
'''
import numpy as np
import pandas as pd

class isochrone:

    def __init__(self, color : np.ndarray, abs_mag : np.ndarray, age : float, metallicity : float):
        """__init__.

        Class to contain attributes for a given isochrone, including 
        its associated error when it is found by the Fitter class

        Parameters
        ----------
        color : np.ndarray
            the "x-axis", G_BP - G_RP color
        abs_mag : np.ndarray
            the "y-axis", M_G absolute magnitude
        metallicity : 
            M_H metallicty of the isochrone
        age : float
            age of isochrone
        """
        
        self.color = color 
        self.abs_mag = abs_mag 
        self.age = age 
        self.metallicity = metallicity

        self.best_fit = False # Set to `True` if Fitter calculates Isochrone to be best fit
        self.__error = None 
        self.chi_square = None
#
#    @property
#    def chi_square(self):
#        return self.chi_square
#    
#    @chi_square.setter
#    def chi_square(self, value):
#        if self.chi_square == None:
#            self.chi_square = value
#        else:
#            raise TypeError("Isochrone chi square attribute already had a value")
#
        
    @property
    def error(self):
        return self.__error
    
    @error.setter
    def error(self, value):
        if self.__error == None:
            self.__error = value
        else:
            raise TypeError("Isochrone error attribute already had a value")
    
    def __str__(self):

        '''
        String representation of Isochrone
        '''
        out_str = 'Isochrone with log age {}, metallicity {}, and is'.format(self.age, self.metallicity)

        if self.best_fit:
            out_str += ' THE BEST FIT!'
        elif not self.best_fit:
            out_str += ' not the best fit.'

        return out_str

    def output_isochrone(self, file_to_name):
        """output_isochrone.

        Writes isochrone as file. Used to save
        best fit isochrone.

        Parameters
        ----------
        file_to_name : str
            Name of file
        """
        N = len(self.color) # length of data points

        color = self.color
        abs_mag = self.abs_mag
        metallicity = self.metallicity
        best_fit = int(self.best_fit)*np.ones(N)
        age = self.age*np.ones(N)

        df_out = pd.DataFrame({'color' : color, 'abs_mag' : abs_mag, 'metallicity': metallicity, 'best_fit' : best_fit, 'age' : age})

        # TODO: Allow user to set their own directory
        df_out.to_csv('/Users/cam/Desktop/astro_research/orion/orion_populations/best_fit_isochrones/' + file_to_name)

