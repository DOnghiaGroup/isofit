'''
Cameren Swiggum (swiggum2@wisc.edu)

Class to represent data from Parsec/Colibri 
Isochroens. Found at http://stev.oapd.inaf.it/cgi-bin/cmd
'''
import numpy as np
import pandas as pd

class isochrone:
    '''
    Class to contain attributes for a given isochrone, including 
    its associated error when it is found by the Fitter class

    Attributes:
        color: the "x-axis", G_BP - G_RP color
        abs_mag: the "y-axis", M_G absolute magnitude
        metallicty: M_H metallicty of the isochrone
        age: age of isochrone
        best_fit: Set to True if Fitter calculates Isochrone to be best fit
        __error: Error found when Fitter calculates square error
    '''

    def __init__(self, color : np.ndarray, abs_mag : np.ndarray, age : float, metallicity : float):
        '''
        init method that defines attributes and sets error to null 
        ''' 
        
        self.color = color 
        self.abs_mag = abs_mag 
        self.age = age 
        self.metallicity = metallicity
        self.best_fit = False 
        self.__error = None 
        
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
    

