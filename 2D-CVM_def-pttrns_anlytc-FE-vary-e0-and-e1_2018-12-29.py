# -*- coding: utf-8 -*-
####################################################################################################
# Alianna J. Maren
# Computing configuration variables for the Cluster Variation Method
####################################################################################################
# Import the following Python packages

import random
import itertools
import numpy as np
import pylab
import matplotlib
from math import exp
from math import log
from matplotlib import pyplot as plt
from random import randrange, uniform #(not sure this is needed, since I'm importing random)


####################################################################################################
####################################################################################################
#
# Detailed code documentation is JUST ABOVE main(), at the very end of this program. 
#
####################################################################################################
####################################################################################################
#
# Code Update Log:
#   12/12/2018: Tweaked prints and plots. 
#     - Added a new pattern - a single big round of x1 in the middle of an x2 sea.
#
#   11/25/2018: Substantial updates to style and structure, but no real change to core functionality
#     -  Changed variable names and procedure / function names to accord with Python conventions
#     -  Re-organized code for improved readability
#     -  Tighted up documentation; created better interfaces
#
#   1/20/2018: Initial code development for Perturbation Experiments 1.2 (folloiwng 1.0 and 1.1)
#     - Original code included basic computation for configuration variables but MISSING the
#       computation for the VERTICAL z(i); however, the vertical w(i) were computed, as well as all
#       the y(i). 
#     - Code results for configuration variables should be approximately correct, based on 
#       horizontal-only values for the z(i).
#     - Code included computation of thermodynamic variables, and original version included
#       a select-and-test approach for free energy minimization given an initial (randomly-generated)
#       starting point. 
#
# This specific version of the code computes a randomly-generated distribution of x1 / x2 values, 
#   dependent on the h-parameter. Then, it computes the configuration variables for the grid. Based on 
#   these, it then computes the entropy, enthalpy, and free energy values. 
#
# The crucial equations are as follows (taken from AJM's 2014 paper, "The Cluster Variation Method II: 
#   2-D Grid of Zigzag Chains":
#   h = exp(beta*epsilon/4) & lambda = 0 (Beginning of Appendix B, replicating Eqn. 2-16.)
#   We can set beta = Boltzmann's constant = 1. 
#   Thus, eps1 = epsilon = 4*log(h)
#   For the equilibrium case (which is where we have an analytic solution), eps0 = 0. 
# Thus, x1 = x2 = 0; h controls the distribution among the z, w, & y values.
# At equilibrium, when eps1 = 0, z1 = z6 = z3 = z4 = 0.125; z2 = z5 = 0.25 (due to degeneracy). 
#
#   y1 = z1 + 0.5*(0.5 - z1 - z3)
#   y2 = z3 + 0.5*(0.5 - z1 - z3)
#   z3 = (h*h - 3.0)*(h*h + 1.0)/(8.0*(h*h*h*h - 6.0*h*h + 1.0))      App. B, Eqn. 29
#   z1 = (1.0 - 3.0*h*h)*(h*h + 1.0)/(8.0*(h*h*h*h - 6.0*h*h + 1.0))  App. B, Eqn. 30
#
#
####################################################################################################
####################################################################################################
#
# Procedure to welcome the user and identify the code
#
####################################################################################################
####################################################################################################

def welcome ():

    print()    
    print()    
    print()
    print()
    print()    
    print()
    print( '******************************************************************************')
    print()
    print( 'Welcome to the 2-D Cluster Variation Method')
    print( 'Version 1.2, 01/07/2018, A.J. Maren')
    print( '  and updated 12/27/2018, by A.J. Maren')
    print( 'This version works with a randomly-generated array to identify the')
    print( '  configuration variables,')
    print( 'and finds the various thermodynamic values associated with different')
    print( '  user-specified values for epsilon0 and epsilon1.' )
    print() 
    print( 'By changing parameters in the main code, the user can select:' )
    print(  '  (O) Randomly generating (and then improving) an array, or' )
    print(  '  (1 .. N) Selecting a pre-stored array' )
    print() 
    print( 'For comments, questions, or bug-fixes, contact: alianna.maren@northwestern.edu')
    print( 'Alternate email address: alianna@aliannajmaren.com')
    print()
    print( '  NOTE: In these calculations, x1 = A (units are at value 1),')
    print( '                           and x2 = B (units are at value 0).')
    print()
    print( '******************************************************************************')
    print()
    return()


####################################################################################################
####################################################################################################
#
# Function to obtain the array size specifications (currently DEFINED for the user; not a choice)
#
# Note: The code is ONLY set up to work with a grid consisting of an EVEN number of rows
#
####################################################################################################
####################################################################################################

def obtain_array_size_specs ():
    
#    x = input('Enter array_length: ')
#    array_length = int(x)
#    print 'array_length is', array_length  
          
#    x = input('Enter layers: ')
#    layers = int(x)
#    print 'layers is', layers
 
#   NOTE: The system is designed to work with an even number of rows, e.g. layers must be an even number
    
    
    array_length = 16
    layers = 16
            
                
    array_size_list = (array_length, layers)  
    return (array_size_list)  

# ************************************************************************************************ #
#
# Pattern Storage 
#
#   This program allows the user to access various pre-stored 16x16 patterns, exemplifying:
#     - Scale-free
#     - Rich club
#     - and potentially other topologies. 
#
#   I'm going to allow the user to select a specified pattern from a pattern-selection module 
#     (still to be written) that will be called from __main__
# 
#   Since grid size is pre-determined (16x16), each pattern is called by specifiying individual rows.
#
# ************************************************************************************************ #


####################################################################################################
####################################################################################################
#
# Function to obtain the choice of a randomly-generated pattern or select a prestored pattern
#
####################################################################################################
####################################################################################################

def obtain_pattern_selection():

    pattern_select = 5
    return(pattern_select)

####################################################################################################
####################################################################################################
#
# Function to obtain a row of 2-D CVM data - PRESTORED pattern (currently part of a 16x16 grid)
#
####################################################################################################
####################################################################################################


def obtainGridRow (rowNum, patternSelect):


# Note: This is some vestigial data, from eary development stages
#       This 4x8 pattern corresponds to an illustration in an early paper. 
#
# Note: 4 rows of 8 units each - this is the small-scale, 2-D equilibrium test case 
#    rowArray0 =  [1,1,1,0,1,1,1,0] # Row 0 - top row 
#    rowArray1 =  [1,0,0,0,1,0,0,0] # Row 1 - second row (counting down from the top)
#    rowArray2 =  [1,1,1,0,1,1,1,0] # Row 2 - third row (counting down from the top) 
#    rowArray3 =  [1,0,0,0,1,0,0,0] # Row 3 - fourth row (counting down from the top)

# Note: 16 rows of 16 units each - this is the 2-D scale-free equilibrium test case 
# Equilibrium scale-free; rows 0 - 7 
    if patternSelect == 1:
        rowArray0 =   [1,0,0,0,1,0,0,1, 0,0,1,1,1,1,0,0] # Row 0 - top row 
        rowArray1 =   [1,0,1,0,0,0,1,1, 0,1,1,1,0,0,1,0] # Row 1 - second row (counting down from the top)
        rowArray2 =   [0,0,1,1,0,1,0,1, 1,0,1,0,0,1,1,1] # Row 2 - third row (counting down from the top) 
        rowArray3 =   [1,0,1,0,1,1,0,0, 1,0,0,0,1,1,1,0] # Row 3 - fourth row (counting down from the top)
        rowArray4 =   [1,1,0,0,0,1,0,1, 0,0,0,1,1,1,1,0] # Row 4 - fifth row (counting down from the top)
        rowArray5 =   [1,1,0,1,0,0,1,1, 1,0,0,1,1,1,0,0] # Row 5 - sixth row (counting down from the top)
        rowArray6 =   [1,1,0,1,0,0,0,1, 1,0,1,0,1,1,1,0] # Row 6 - seventh row (counting down from the top) 
        rowArray7 =   [1,0,1,1,0,1,0,0, 0,1,1,0,0,0,0,0] # Row 7 - eighth row (counting down from the top)  
    
# Rich Club
    if patternSelect == 2:
        rowArray0 =   [1,1,1,1,0,0,0,0, 0,0,0,1,1,1,1,1] # Row 0 - top row 
        rowArray1 =   [1,1,0,0,0,0,0,0, 0,1,1,1,1,1,1,1] # Row 1 - second row (counting down from the top)
        rowArray2 =   [1,1,1,0,0,0,0,0, 0,1,1,1,1,1,1,1] # Row 2 - third row (counting down from the top) 
        rowArray3 =   [1,1,0,0,0,0,0,0, 0,1,1,1,1,1,1,1] # Row 3 - fourth row (counting down from the top)
        rowArray4 =   [1,1,0,0,0,0,0,0, 0,0,0,1,1,1,1,1] # Row 4 - fifth row (counting down from the top)
        rowArray5 =   [1,1,0,0,0,0,0,0, 0,0,0,1,1,1,1,1] # Row 5 - sixth row (counting down from the top)
        rowArray6 =   [1,1,1,0,0,0,0,0, 0,0,0,0,1,1,1,1] # Row 6 - seventh row (counting down from the top) 
        rowArray7 =   [1,1,1,0,0,0,0,0, 0,0,0,0,0,1,1,1] # Row 7 - eighth row (counting down from the top)                         
                                                                
# Rich club with x1 < 0.5; rows 0 - 7 
    if patternSelect == 3:
        rowArray0 =   [0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0] # Row 0 - top row 
        rowArray1 =   [0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0] # Row 1 - second row (counting down from the top)
        rowArray2 =   [0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0] # Row 2 - third row (counting down from the top) 
        rowArray3 =   [0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0] # Row 3 - fourth row (counting down from the top)
        rowArray4 =   [0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0] # Row 4 - fifth row (counting down from the top)
        rowArray5 =   [0,0,0,0,0,0,1,1, 1,0,0,0,0,0,0,0] # Row 5 - sixth row (counting down from the top)
        rowArray6 =   [0,0,0,0,0,1,1,1, 1,1,0,0,0,0,0,0] # Row 6 - seventh row (counting down from the top) 
        rowArray7 =   [0,0,0,0,1,1,1,1, 1,1,1,0,0,0,0,0] # Row 7 - eighth row (counting down from the top)

# Stretched-out Rich club with x1 < 0.5; rows 0 - 7 
    if patternSelect == 4:
        rowArray0 =   [0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0] # Row 0 - top row 
        rowArray1 =   [0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0] # Row 1 - second row (counting down from the top)
        rowArray2 =   [0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0] # Row 2 - third row (counting down from the top) 
        rowArray3 =   [0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0] # Row 3 - fourth row (counting down from the top)
        rowArray4 =   [0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0] # Row 4 - fifth row (counting down from the top)
        rowArray5 =   [0,0,0,1,1,0,0,0, 0,0,1,0,0,0,0,0] # Row 5 - sixth row (counting down from the top)
        rowArray6 =   [0,0,0,1,1,1,0,0, 0,1,1,0,0,0,0,0] # Row 6 - seventh row (counting down from the top) 
        rowArray7 =   [0,0,0,0,1,1,1,1, 1,1,1,0,0,0,0,0] # Row 7 - eighth row (counting down from the top)

# Rich Club
    if patternSelect == 5:
        rowArray0 =   [1,1,1,0,0,0,0,0, 0,0,0,1,1,1,1,1] # Row 0 - top row 
        rowArray1 =   [1,0,0,0,0,0,0,0, 0,1,1,1,1,1,1,1] # Row 1 - second row (counting down from the top)
        rowArray2 =   [1,1,0,0,0,0,0,0, 0,1,1,1,1,1,1,1] # Row 2 - third row (counting down from the top) 
        rowArray3 =   [1,1,0,0,0,0,0,0, 0,1,1,1,1,1,1,1] # Row 3 - fourth row (counting down from the top)
        rowArray4 =   [1,1,0,0,0,0,0,0, 0,0,0,1,1,1,1,1] # Row 4 - fifth row (counting down from the top)
        rowArray5 =   [1,1,0,0,0,0,0,0, 0,0,0,1,1,1,1,1] # Row 5 - sixth row (counting down from the top)
        rowArray6 =   [1,1,0,0,0,0,0,0, 0,0,0,0,1,1,1,1] # Row 6 - seventh row (counting down from the top) 
        rowArray7 =   [1,1,0,0,0,0,0,0, 0,0,0,0,0,1,1,1] # Row 7 - eighth row (counting down from the top)                         
                                      
# Non-equilibrium scale-free; rows 1 - 8 (9 - 16 in graph)
#    rowArray0 =   [1,0,0,0,1,0,0,0, 0,0,1,1,1,1,0,0] # Row 0 - top row 
#    rowArray1 =   [1,0,1,0,0,0,0,0, 0,1,1,1,0,0,0,0] # Row 1 - second row (counting down from the top)
#    rowArray2 =   [0,0,1,1,0,1,0,0, 0,0,1,0,0,1,0,0] # Row 2 - third row (counting down from the top) 
#    rowArray3 =   [1,0,1,0,1,1,0,0, 0,0,0,0,1,1,0,0] # Row 3 - fourth row (counting down from the top)
#    rowArray4 =   [1,1,0,0,0,1,0,1, 0,0,0,0,0,1,1,0] # Row 4 - fifth row (counting down from the top)
#    rowArray5 =   [1,1,0,1,0,0,1,1, 1,0,0,0,0,0,1,0] # Row 5 - sixth row (counting down from the top)
#    rowArray6 =   [0,0,0,1,0,0,0,1, 1,0,1,0,0,0,0,0] # Row 6 - seventh row (counting down from the top) 
#    rowArray7 =   [0,0,1,1,0,0,0,0, 0,1,1,0,0,0,0,0] # Row 7 - eighth row (counting down from the top)                                              
                                                                                                            
# Second non-equilibrium scale-free set (two side clusters removed); rows 1 - 8 (9 - 16 in graph)
#    rowArray0 =   [1,0,0,0,1,0,0,0, 0,0,1,1,1,1,0,0] # Row 0 - top row 
#    rowArray1 =   [1,0,1,0,0,0,0,0, 0,1,1,1,0,0,0,0] # Row 1 - second row (counting down from the top)
#    rowArray2 =   [0,0,1,1,0,1,0,0, 0,0,1,0,0,1,0,0] # Row 2 - third row (counting down from the top) 
#    rowArray3 =   [0,0,1,0,1,1,0,0, 0,0,0,0,1,1,0,0] # Row 3 - fourth row (counting down from the top)
#    rowArray4 =   [0,0,0,0,0,1,0,1, 0,0,0,0,0,1,1,0] # Row 4 - fifth row (counting down from the top)
#    rowArray5 =   [0,0,0,1,0,0,1,1, 1,0,0,0,0,0,1,0] # Row 5 - sixth row (counting down from the top)
#    rowArray6 =   [1,1,0,1,0,0,0,1, 1,0,1,0,0,0,0,0] # Row 6 - seventh row (counting down from the top) 
#    rowArray7 =   [1,0,1,1,0,0,0,0, 0,1,1,0,0,0,0,0] # Row 7 - eighth row (counting down from the top)                                              
                                                                                                                                                                                                                                            
                                                                                                                                                                                    
# Equilibrium scale-free; rows 8 - 15
    if patternSelect == 1:         
        rowArray8 =   [0,0,0,0,0,1,1,0, 0,0,1,0,1,1,0,1] # Row 8 - ninth row 
        rowArray9 =   [0,1,1,1,0,1,0,1, 1,0,0,0,1,0,1,1] # Row 9 - tenth row (counting down from the top)
        rowArray10 =  [0,0,1,1,1,0,0,1, 1,1,0,0,1,0,1,1] # Row 10 - eleventh row (counting down from the top) 
        rowArray11 =  [0,1,1,1,1,0,0,0, 1,0,1,0,0,0,1,1] # Row 11 - twelfth row (counting down from the top)
        rowArray12 =  [0,1,1,1,0,0,0,1, 0,0,1,1,0,1,0,1] # Row 12 - thirteenth row (counting down from the top)
        rowArray13 =  [1,1,1,0,0,1,0,1, 1,0,1,0,1,1,0,0] # Row 13 - fourteenth row (counting down from the top)
        rowArray14 =  [0,1,0,0,1,1,1,0, 1,1,0,0,0,1,0,1] # Row 14 - fifteenth row (counting down from the top) 
        rowArray15 =  [0,0,1,1,1,1,0,0, 1,0,0,1,0,0,0,1] # Row 15 - sixteenth row (counting down from the top)                                                                   

# Rich Club
    if patternSelect == 2:        
        rowArray8 =   [1,1,1,0,0,0,0,0, 0,0,0,0,0,1,1,1] # Row 8 - ninth row 
        rowArray9 =   [1,1,1,1,0,0,0,0, 0,0,0,0,0,1,1,1] # Row 9 - tenth row (counting down from the top)
        rowArray10 =  [1,1,1,1,1,0,0,0, 0,0,0,0,0,0,1,1] # Row 10 - eleventh row (counting down from the top) 
        rowArray11 =  [1,1,1,1,1,0,0,0, 0,0,0,0,0,0,1,1] # Row 11 - twelfth row (counting down from the top)
        rowArray12 =  [1,1,1,1,1,1,1,0, 0,0,0,0,0,0,1,1] # Row 12 - thirteenth row (counting down from the top)
        rowArray13 =  [1,1,1,1,1,1,1,0, 0,0,0,0,0,1,1,1] # Row 13 - fourteenth row (counting down from the top)
        rowArray14 =  [1,1,1,1,1,1,1,0, 0,0,0,0,0,0,1,1] # Row 14 - fifteenth row (counting down from the top) 
        rowArray15 =  [1,1,1,1,1,0,0,0, 0,0,0,0,1,1,1,1] # Row 15 - sixteenth row (counting down from the top) 

# Rich Club with x1<< 0.5
    if patternSelect == 3:        
        rowArray8 =   [0,0,0,0,1,1,1,1, 1,1,1,0,0,0,0,0] # Row 8 - ninth row 
        rowArray9 =   [0,0,0,0,0,1,1,1, 1,1,0,0,0,0,0,0] # Row 9 - tenth row (counting down from the top)
        rowArray10 =  [0,0,0,0,0,0,1,1, 1,0,0,0,0,0,0,0] # Row 10 - eleventh row (counting down from the top) 
        rowArray11 =  [0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0] # Row 11 - twelfth row (counting down from the top)
        rowArray12 =  [0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0] # Row 12 - thirteenth row (counting down from the top)
        rowArray13 =  [0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0] # Row 13 - fourteenth row (counting down from the top)
        rowArray14 =  [0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0] # Row 14 - fifteenth row (counting down from the top) 
        rowArray15 =  [0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0] # Row 15 - sixteenth row (counting down from the top) 

# Stretched-out Rich Club with x1<< 0.5
    if patternSelect == 4:        
        rowArray8 =   [0,0,0,0,1,1,1,1, 1,1,1,0,0,0,0,0] # Row 8 - ninth row 
        rowArray9 =   [0,0,0,1,1,1,0,0, 0,1,1,0,0,0,0,0] # Row 9 - tenth row (counting down from the top)
        rowArray10 =  [0,0,0,1,1,0,0,0, 0,0,1,0,0,0,0,0] # Row 10 - eleventh row (counting down from the top) 
        rowArray11 =  [0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0] # Row 11 - twelfth row (counting down from the top)
        rowArray12 =  [0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0] # Row 12 - thirteenth row (counting down from the top)
        rowArray13 =  [0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0] # Row 13 - fourteenth row (counting down from the top)
        rowArray14 =  [0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0] # Row 14 - fifteenth row (counting down from the top) 
        rowArray15 =  [0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0] # Row 15 - sixteenth row (counting down from the top) 

# Rich Club
    if patternSelect == 5:        
        rowArray8 =   [1,1,0,0,0,0,0,0, 0,0,0,0,0,1,1,1] # Row 8 - ninth row 
        rowArray9 =   [1,1,1,0,0,0,0,0, 0,0,0,0,0,1,1,1] # Row 9 - tenth row (counting down from the top)
        rowArray10 =  [1,1,1,1,0,0,0,0, 0,0,0,0,0,0,1,1] # Row 10 - eleventh row (counting down from the top) 
        rowArray11 =  [1,1,1,1,0,0,0,0, 0,0,0,0,0,0,1,1] # Row 11 - twelfth row (counting down from the top)
        rowArray12 =  [1,1,1,1,1,1,0,0, 0,0,0,0,0,0,1,1] # Row 12 - thirteenth row (counting down from the top)
        rowArray13 =  [1,1,1,1,1,1,0,0, 0,0,0,0,0,1,1,1] # Row 13 - fourteenth row (counting down from the top)
        rowArray14 =  [1,1,1,1,1,1,0,0, 0,0,0,0,0,0,1,1] # Row 14 - fifteenth row (counting down from the top) 
        rowArray15 =  [1,1,1,1,0,0,0,0, 0,0,0,0,1,1,1,1] # Row 15 - sixteenth row (counting down from the top) 


# Non-equilibrium scale-free; rows 9 - 16 (1 - 8 in graph) -- total number of A units reduced by 17 (out of 256)  
#    rowArray8 =   [0,0,0,0,0,1,0,0, 0,0,1,0,1,1,0,1] # Row 0 - ninth row 
#    rowArray9 =   [0,0,0,0,1,1,0,0, 0,0,0,0,1,0,1,1] # Row 1 - tenth row (counting down from the top)
#    rowArray10 =  [0,1,0,0,0,1,0,1, 1,0,0,0,1,0,1,1] # Row 2 - eleventh row (counting down from the top) 
#    rowArray11 =  [0,1,1,0,0,0,1,1, 1,0,1,0,0,0,1,1] # Row 3 - twelfth row (counting down from the top)
#    rowArray12 =  [0,0,1,1,0,0,0,0, 1,0,1,1,0,1,0,1] # Row 4 - thirteenth row (counting down from the top)
#    rowArray13 =  [0,0,1,0,0,1,0,0, 0,0,1,0,1,1,0,0] # Row 5 - fourteenth row (counting down from the top)
#    rowArray14 =  [0,0,0,0,1,1,1,0, 0,0,0,0,0,1,0,1] # Row 6 - fifteenth row (counting down from the top) 
#    rowArray15 =  [0,0,1,1,1,1,0,0, 1,0,0,1,0,0,0,1] # Row 7 - sixteenth row (counting down from the top)                                                                                         

# # Second non-equilibrium scale-free set (two side clusters removed) 
#    rowArray8 =   [0,0,0,0,0,1,0,0, 0,0,1,0,1,1,0,0] # Row 0 - ninth row 
#    rowArray9 =   [0,0,0,0,1,1,0,0, 0,0,0,0,1,0,0,0] # Row 1 - tenth row (counting down from the top)
#    rowArray10 =  [0,1,0,0,0,1,0,1, 1,0,0,0,1,0,0,0] # Row 2 - eleventh row (counting down from the top) 
#    rowArray11 =  [0,1,1,0,0,0,1,1, 1,0,1,0,0,0,0,0] # Row 3 - twelfth row (counting down from the top)
#    rowArray12 =  [0,0,1,1,0,0,0,0, 1,0,1,1,0,1,0,0] # Row 4 - thirteenth row (counting down from the top)
#    rowArray13 =  [0,0,1,0,0,1,0,0, 0,0,1,0,1,1,0,0] # Row 5 - fourteenth row (counting down from the top)
#    rowArray14 =  [0,0,0,0,1,1,1,0, 0,0,0,0,0,1,0,1] # Row 6 - fifteenth row (counting down from the top) 
#    rowArray15 =  [0,0,1,1,1,1,0,0, 1,0,0,1,0,0,0,1] # Row 7 - sixteenth row (counting down from the top)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
                                                                                                                                    
    if rowNum == 0: rowArray = rowArray0
    if rowNum == 1: rowArray = rowArray1 
    if rowNum == 2: rowArray = rowArray2
    if rowNum == 3: rowArray = rowArray3
    if rowNum == 4: rowArray = rowArray4
    if rowNum == 5: rowArray = rowArray5 
    if rowNum == 6: rowArray = rowArray6
    if rowNum == 7: rowArray = rowArray7
    if rowNum == 8: rowArray = rowArray8
    if rowNum == 9: rowArray = rowArray9 
    if rowNum == 10: rowArray = rowArray10
    if rowNum == 11: rowArray = rowArray11
    if rowNum == 12: rowArray = rowArray12
    if rowNum == 13: rowArray = rowArray13 
    if rowNum == 14: rowArray = rowArray14
    if rowNum == 15: rowArray = rowArray15
                            
    return (rowArray)

####################################################################################################
#
# Procedure to print out the 2-D CVM grid size specifications
#
####################################################################################################
    
def print_debug_status (debug_print_off):

    if not debug_print_off:
        print()
        print( 'Debug printing is on')  #debug_print_off false
    else: 
        print( 'Debug printing is off') #debug_print_off true   
    print()
    print ('------------------------------------------------------------------------------')
    print()
    return

####################################################################################################
#
# Procedure to print out the 2-D CVM grid size specifications
#
####################################################################################################

def print_grid_size_specs ():
    
    print()
    print( 'Grid size specifications: ')
    print( '  This 2-D CVM process works with a matrix of M x L units, where:' )
    print( '    M (the array_length is):', array_length )
    print( '    L (the layers is):      ', array_layers )  
    print()  
    print ('------------------------------------------------------------------------------')
    print() 
    return

####################################################################################################
#
# Procedure to print the parameters for the run
#
####################################################################################################
    
def print_run_parameters (h0, h_incr, h_range):
 
    start_range = h0 + h_incr
    end_range = h0 + h_incr*h_range 
    print()
    print( ' The parameters for this run are: ')
    print( '   Initial interaction enthalpy parameter h is:', h0)
    print( ' The run will begin with an interaction enthalpy of h at: %.3f' % start_range)
    print( '   up through a final value of:  %.3f' % end_range)
    print( ' The increment in h will be: ', h_incr)
    print()
    print()
    return

####################################################################################################
#
# Procedure to print the parameters for the run
#
####################################################################################################

def print_pattern_selection(pattern_select):

    print ()
    if pattern_select == 0:
        print ( 'The selected pattern for this run is randomly-generated.')
    if pattern_select == 1:
        print ( 'The selected pattern for this run is a scale-free-like pattern.')        
    if pattern_select == 2:
        print ( 'The selected pattern for this run is a rich-club-like pattern.') 
    if pattern_select == 3:
        print ( 'The selected pattern for this run is a very small rich-club-like pattern.')        
        
    print()
    print ('------------------------------------------------------------------------------')
    print ()
    return()


####################################################################################################
#
# Procedure to print out the 2-D CVM grid
#
####################################################################################################

def print_grid (unit_array):
       
    for i in range (0, pairs):
        if i<5: # This puts the single-decimal rows a little to the left
            actualEvenRowNum = 2*i
            print( 'Row ', actualEvenRowNum, ':  ', end =" ") 
            for j in range(0, array_length):
                print( unit_array[actualEvenRowNum,j], end =" ") 
            print ()
            actualOddRowNum = 2*i+1
            print( 'Row  ', actualOddRowNum, ':  ', end =" ") 
            for j in range(0,array_length):
                print( unit_array[actualOddRowNum,j], end =" ") 
            print  ()
        else:
            actualEvenRowNum = 2*i
            print( 'Row', actualEvenRowNum, ':  ', end =" ") 
            for j in range(0,array_length):
                print( unit_array[actualEvenRowNum,j], end =" ") 
            print ()
            actualOddRowNum = 2*i+1
            print( 'Row ', actualOddRowNum, ':  ', end =" ") 
            for j in range(0,array_length):
                print( unit_array[actualOddRowNum,j], end =" ") 
            print  ()           
    print ()    
    return


####################################################################################################
#
# Procedure to print out the x1 and x2 variables
#
#    Inputs:    array_size_list: a list of two integers; array_length and layers
#               h: the interaction enthalpy parameter
#
####################################################################################################

def print_x_result (x1, x2, x1_total, x2_total, x1_target, max_x_dif):

# Print the locally-computed values for x1 and x2; these are not passed back to _main__.     
 
    print() 
    print( 'The distribution among states A and B (x1 and x2) units is:' )    
    print( "  ( A ) x1_total =", x1_total, "( B ) x2_total =", x2_total, ' for a total of ', x1_total + x2_total, ' units.' )   
    print() 
    print( ' The fractional values for x are:  x1 = %.4f'  % x1, ' and x2 = %.4f'  % x2)
    print()  
    print( ' The actual value for x1 is %.4f'  % x1, ' and the desired value for x1 is ', x1_target) 
    delta_x = x1 - x1_target
    print( '   ==>> The difference (x1 - x1_target) is %.4f'  % delta_x)
    print( ' The acceptable difference between the two values is ', max_x_dif)
    if abs(delta_x) > max_x_dif:
        if delta_x > 0:
            print( '   so we see that x1 is too large, and we want to decrease x1.')
        if delta_x < 0:
            print( '   so we see that x1 is too small, and we want to increase x1.')
    else:
        print(' The difference between the actual and the target is within accepted bounds.')            
    print()
    print ('------------------------------------------------------------------------------')
    print ()
    return


####################################################################################################
#
# Procedure to print out the fractional values for the configuration variables
#
#    Inputs:    config_vars_frac_list: a list of the configuration variables as fraction values
#
####################################################################################################
    
def print_config_vars_fraction_vals (config_vars_frac_list, pattern):
                       
    x1 = config_vars_frac_list[0]
    x2 = config_vars_frac_list[1]
    
    y1 = config_vars_frac_list[2]
    y2 = config_vars_frac_list[3]
    y3 = config_vars_frac_list[4] 
    
    sumY = y1 + 2.*y2 + y3

    w1 = config_vars_frac_list[5]
    w2 = config_vars_frac_list[6]
    w3 = config_vars_frac_list[7]       

    sumW = w1 + 2.*w2 + w3

    z1 = config_vars_frac_list[8]
    z2 = config_vars_frac_list[9]
    z3 = config_vars_frac_list[10] 
    z4 = config_vars_frac_list[11]
    z5 = config_vars_frac_list[12]
    z6 = config_vars_frac_list[13] 
              
    sumZ = z1 + 2.*z2 + z3 + z4 + 2.*z5 + z6    
    

    
    print() 
    print( ' For pattern number ', pattern, ', the configuration variables have the following values:')        
    print()
    print( '        x1 = %.4f'  % (x1), '    x2 = %.4f'  % (x2))
    print( ' y1 = %.4f'  % (y1), ' y2Total = %.4f'  % (2.*y2), ' y3 = %.4f'  % (y3) ) #, ' sumY = %.4f'  % (sumY)                                                 
    print( ' w1 = %.4f'  % (w1), ' w2Total = %.4f'  % (2.0*w2), ' w3 = %.4f'  % (w3) ) #, ' sumW = %.4f'  % (sumW) 

    print()
    print( ' z1 = %.4f'  % (z1), ' z2Total = %.4f'  % (2.0*z2), ' z3 = %.4f'  % (z3) )
    print( ' z6 = %.4f'  % (z6), ' z5Total = %.4f'  % (2.0*z5), ' z4 = %.4f'  % (z4) ) #, 'sumZ = %.4f'  % (sumW) 
    return


####################################################################################################
#
# Procedure to print out the full set of configuration variables
#
#    Inputs:    config_vars_list
#
####################################################################################################

def print_config_vars_comparison (config_vars_frac_list, zero_activation_analytic_config_vars_list, h):

    total_units = float(array_length*array_layers)
#    printConfigVars (config_vars_list)      
                        
    x1 = config_vars_frac_list[0]
    x2 = config_vars_frac_list[1]  
    y1 = config_vars_frac_list[2]
    y2 = config_vars_frac_list[3]
    y3 = config_vars_frac_list[4] 
    w1 = config_vars_frac_list[5]
    w2 = config_vars_frac_list[6]
    w3 = config_vars_frac_list[7]       
    z1 = config_vars_frac_list[8]
    z2 = config_vars_frac_list[9]
    z3 = config_vars_frac_list[10]
    z4 = config_vars_frac_list[11]
    z5 = config_vars_frac_list[12]
    z6 = config_vars_frac_list[13] 
    
    x1_analyt = zero_activation_analytic_config_vars_list[0]
    x2_analyt = zero_activation_analytic_config_vars_list[1]  
    y1_analyt = zero_activation_analytic_config_vars_list[2]
    y2_analyt = zero_activation_analytic_config_vars_list[3]*2.0
    y3_analyt = zero_activation_analytic_config_vars_list[4] 
    w1_analyt = zero_activation_analytic_config_vars_list[5]
    w2_analyt = zero_activation_analytic_config_vars_list[6]
    w3_analyt = zero_activation_analytic_config_vars_list[7]       
    z1_analyt = zero_activation_analytic_config_vars_list[8]
    z2_analyt = zero_activation_analytic_config_vars_list[9]
    z3_analyt = zero_activation_analytic_config_vars_list[10]
    z4_analyt = zero_activation_analytic_config_vars_list[11]
    z5_analyt = zero_activation_analytic_config_vars_list[12]
    z6_analyt = zero_activation_analytic_config_vars_list[13] 
    h = zero_activation_analytic_config_vars_list[14]
        

    print()
    print(' ----------------------------------------------------------------------')    
    print()
    print(' Comparing the actual configuration variables for the resultant grid ')
    print('   versus those that would result if the analytic case, where eps0 = 0.0, were true.')
    print(' Note that h does not have to be 1.0; i.e., eps1 does not have to be 0.')
    print(' However, the further that h is from 1.0, the more likely it is that the ')
    print('   analytic values of the configuration variables will differ from the actual.')
    
    sumX = x1+x2  
    print()
    print( ' For h =  %.4f' % h)
    print( '             Grid values             Analytic Values')
    print() 
    print( '                  x1 = %.4f'  % x1,  '   x1_analyt = %.4f' % x1_analyt )
    print( '                  x2 = %.4f'  % x2,  '   x2_analyt = %.4f' % x2_analyt )
    print( '   Sum of the x(i) = %.4f'  % sumX )
    print()    

    sumY = y1+2*y2+y3    
    print()
    print( "For the y(i) variables:" )
    print( '            (A-A) y1 = %.4f'  % y1,  '   y1_analyt = %.4f' % y1_analyt )
    print( '    (A-B & B-A) 2*y2 = %.4f'  % y2,  '   y2_analyt = %.4f' % y2_analyt , '(total y2)' )  
    print( '            (B-B) y3 = %.4f'  % y3,  '   y3_analyt = %.4f' % y3_analyt )
    print( '     Sum of the y(i) = %.4f'  % sumY )
    print() 
    
    sumW = w1+2*w2+w3 
    print()
    print( "Totals for the w(i) variables:" )
    print( '          (A---A) w1 = %.4f'  % w1,  '   w1_analyt = %.4f' % y1_analyt )
    print( '  (A---B & B---A) w2 = %.4f'  % w2,  '   w2_analyt = %.4f' % y2_analyt , '(total w2)' )    
    print( '          (B---B) w3 = %.4f'  % w3,  '   w3_analyt = %.4f' % y3_analyt )
    print( '     Sum of the w(i) = %.4f'  % sumW )
    print()         

    sumZ = z1+2*z2+z3 +z4+2*z5+z6
    print()
    print( "  Totals for the z(i) variables:" )
    print( '          (A-A-A) z1 = %.4f'  % z1,  '   z1_analyt = %.4f' % z1_analyt )
    print( '  (A-A-B & B-A-A) z2 = %.4f'  % z2,  '   z2_analyt = %.4f' % z2_analyt , '(total z2)' )   
    print( '          (A-B-A) z3 = %.4f'  % z3,  '   z1_analyt = %.4f' % z3_analyt )
    print( '          (B-A-B) z4 = %.4f'  % z4,  '   z4_analyt = %.4f' % z4_analyt )
    print( '  (A-B-B & B-B-A) z5 = %.4f'  % z5,  '   z5_analyt = %.4f' % z5_analyt , '(total z5)' ) 
    print( '          (B-B-B) z6 = %.4f'  % z6,  '   z6_analyt = %.4f' % z6_analyt )        
    print( '     Sum of the z(i) = %.4f'  % sumZ )
    print()      
    print() 
    return

####################################################################################################
####################################################################################################
#
# Procedure to print the thermodyanmic values for a given the interaction enthalpy paramter h.
#
#    Inputs: 
#           h: the interaction enthalpy parameter
#           sys_vals_list: contains h and the thermodynamic quantities  
#
####################################################################################################
####################################################################################################
    
def print_thermodynamic_values (h, eps0, sys_vals_list):

# Note that: sys_vals_list = (negS, enthalpy0, enthalpy1, free_energy)    
    negS        = sys_vals_list[0]
    enthalpy1   = sys_vals_list[1]    
    enthalpy1   = sys_vals_list[2]
    free_energy = sys_vals_list[3] 
    epsilon1    = log(h)/2.0
    
    print()  
    print() 
    print(' The computed thermodynamic quantities:' ) 
    print('   The activation enthalpy parameter epsilon0 =  %.4f' % (eps0))          
    print( '    h = %.4f' % (h), '    epsilon1 = %.4f' % (epsilon1) )
    print( '    negative entropy s = %.4f' % (negS) )
    print( '    enthalpy = %.4f' % (enthalpy1) )     
    print( '    free energy = %.4f' % (free_energy) )
    print() 
    print() 
    return


####################################################################################################
####################################################################################################
#
# Procedure to print the thermodyanmic values for a given the interaction enthalpy paramter h.
#
#    Inputs: 
#           h: the interaction enthalpy parameter
#           sys_vals_list: contains h and the thermodynamic quantities  
#
####################################################################################################
####################################################################################################
    
def print_one_set_of_thermodynamics (pattern_select, eps0, h_array, h_range, 
                                       neg_S_array, f_eps0_array, f_eps1_array, f_energy_array):   
    

# Note that: sys_vals_list = (negS, enthalpy0, enthalpy1, free_energy)    

    print()  
    print() 
    print(' The computed thermodynamic quantities for the single case where epsilon0 is', eps0, ':' ) 
    print('   Epsilon1 is found as log(h)/2.0, which is from Kikuchi & Brush (1967)')          
    print()    
    print( '  The negEntropy (-S) is constant and is given as: %.4f' % (neg_S_array[0]) )
    print( )

    for i in range (0, h_range, 1):        
        h = h_array[i]    
        epsilon1    = log(h_array[i])/2.0 
        enthalpy = f_eps1_array[i]
        free_energy = f_energy_array[i]
        print( '    h = %.4f' % (h), '    epsilon1 = %.4f' % (epsilon1) , '    enthalpy = %.4f' % (enthalpy), '    free energy = %.4f' % (free_energy) )
    print() 
    return







####################################################################################################
####################################################################################################
#
# Procedure to print the thermodyanmic values for a given the interaction enthalpy paramter h.
#
#    Inputs: 
#           h: the interaction enthalpy parameter
#           sys_vals_list: contains h and the thermodynamic quantities  
#
####################################################################################################
####################################################################################################
    
def print_set_of_three_thermodynamics (pattern_select, eps0a, eps0b, eps0c, h_array, h_range, 
                                       neg_S_array1, f_eps0_array1, f_eps1_array1, f_energy_array1, 
                                       neg_S_array2, f_eps0_array2, f_eps1_array2, f_energy_array2,
                                       neg_S_array3, f_eps0_array3, f_eps1_array3, f_energy_array3):   
    

# Note that: sys_vals_list = (negS, enthalpy0, enthalpy1, free_energy)    

    print()  
    print( '------------------------------------------------------------------------------')
    print()
    print(' (From print_set_of_three_thermodynamics:)')
    print()    
    print(' The computed thermodynamic quantities:' )           
    print()
    print( '    negative entropy s = %.4f' % (neg_S_array1[0]) )
    print()
    print('  For the case where eps0 = ', eps0a)
    print( '   h         eps1    activ-enth   interact-enth   negS      FE')
    for i in range (0, h_range, 1):        
        h = h_array[i]    
        epsilon1    = log(h_array[i])/2.0 
        enthalpy0 = f_eps0_array1[i]
        enthalpy1 = f_eps1_array1[i]
        negS = neg_S_array1[i]
        free_energy = f_energy_array1[i]
        print( '  %.4f' % (h), '   %.4f' % (epsilon1) , '    %.4f' % (enthalpy0), '       %.4f' % (enthalpy1),  '    %.4f' % (negS), '  %.4f' % (free_energy) )
    print() 
    
    print('  For the case where eps0 = ', eps0b)
    print( '   h         eps1    activ-enth   interact-enth   negS      FE')
    for i in range (0, h_range, 1):        
        h = h_array[i]    
        epsilon1    = log(h_array[i])/2.0 
        enthalpy0 = f_eps0_array2[i]
        enthalpy1 = f_eps1_array2[i]
        negS = neg_S_array2[i]        
        free_energy = f_energy_array2[i]
        print( '  %.4f' % (h), '   %.4f' % (epsilon1) , '    %.4f' % (enthalpy0), '       %.4f' % (enthalpy1),  '    %.4f' % (negS), '  %.4f' % (free_energy) )
    print()     

    print()
    print('  For the case where eps0 = ', eps0c)
    print( '   h         eps1    activ-enth   interact-enth   negS      FE')
    for i in range (0, h_range, 1):        
        h = h_array[i]    
        epsilon1    = log(h_array[i])/2.0 
        enthalpy0 = f_eps0_array3[i]        
        enthalpy1 = f_eps1_array3[i]
        negS = neg_S_array3[i]        
        free_energy = f_energy_array3[i]
        print( '  %.4f' % (h), '   %.4f' % (epsilon1) , '    %.4f' % (enthalpy0), '       %.4f' % (enthalpy1),  '    %.4f' % (negS), '  %.4f' % (free_energy) )
    print() 


    return


####################################################################################################
####################################################################################################
#
# Procedure to print a statement about the thermodynamic plots.  
#
####################################################################################################
####################################################################################################


def print_thermodynamic_plot_statement(pattern_select):
    
    print()
    print( '------------------------------------------------------------------------------')
    print()
    print( ' (From print_thermodynamic_plot_statement:)')
    print() 
    print( ' The thermodynamic quantities plot vs. h' )
    print( ' - The negEntropy is in blue, ' )
    print( ' - The per-unit enthalpy is in cyan (and is zero if epsilon0 = 0), ' )
    print( ' - The interaction enthalpy (eps1*(2*y2 - y1 - y3)) is in maroon,' ) 
    print( ' - The free energy is in red.'  ) 
    print()
    
    if not explanation_thermodynamic_plot_off:
        if pattern_select == 0:
            print( ' If a randomly-generated grid has been used, then the distribution of nodes will ')
            print( '    be approximately random, and thus the configuration variables will take on the')          
            print( '    expected probabilistic distributions.')
            print( ' If the interaction enthalpy epsilon1 = 0 (h = 1) then the interaction enthalpy')
            print( '    will be approximately zero:' )
            print( '    2*epsilon1*(2*y2 - y1 - y3) approx== 0.')
            print()
            print( ' In this situation, the interaction enthalpy, shown in maroon, is approx zero.' )
            print( '   Any non-zero values for the interaction enthalpy (most noticeable at higher')
            print( '     values for h) are the result of having x1 being not precisely at x1 = 0.5;')
            print( '     this is due to random pattern generation. As a result, the other configuration values')
            print( '     take on values so that 2*y2 - y1 - y3 is not precisely zero,')
            print( '     so an interaction enthalpy appears.' )
            print( ' The free energy is shown in red.')
            print( ' The entropy is shown in blue.')
            print( '   For the case where the activation enthalpy epsilon0 > 0, there is an '), 
            print( '     activation enthalpy term, and the free energy is not equal to the entropy.' )
            print( '   For the case where the activation enthalpy epsilon0 = 0, there is no ')
            print( '     activation enthalpy, and if epsilon1 also = 0, there is no interaction enthalpy.')
            print( '   In this case, since the free energy curve is printed after the entropy curve, and overlies it, ')
            print('      the entropy curve may not be visible.')
            print('    Trust me; it is there.')
        
        if pattern_select > 0:
            print( ' We have selected a specific pattern for the 2-D CVM grid.')
            print( ' This grid is a designed, rather than random, pattern. ')
            print( ' Because we are not altering the nodes (we are not seeking free energy minimization), ')
            print( '   we will not see any change in the entropy - because the entropy depends strictly on ')
            print( '   the node configurations, and we are not changing those. ')
            print( ' Note that the entropy (blue line) is a flat line.')
            print()
            print( ' We are seeing the interaction enthalpy term (maroon) as a linear expression:')
            print( '   As we vary epsilon1, the interaction enthalpy will vary (linearly).')
            print()
            print( ' Further, the free energy will vary linearly with the interaction enthalpy,')
            print( '   as the interaction enthalpy is the only thing that is changing.')
            print()
            print( ' There is no special meaning to be derived from this, as we are not changing the node')
            print( '   compositions to achieve a free energy minimum.')
            print( ' This is strictly for illustration.')
            print()

 
       
    return


####################################################################################################
####################################################################################################
#
# Procedure to plot the thermodyanmic values versus the interaction enthalpy paramter h.
#
#    Inputs:    
#           h_array: the range of h values, where h =exp(2*epsilon1)
#           neg_S_array: array of the negative entropy values
#           f_eps0_array: array of the activation enthalpies; will be zeros if epsilon0 = 0
#           f_eps1_array: array of the interaction enthalpies; depends on h (actually, on epsllon1)
#           f_energy_array: array of the free energy values, computed directly using the 
#               Kikuchi-Brush equation    
#
####################################################################################################
####################################################################################################


def plot_thermodynamic_vals (h_array, neg_S_array, f_eps0_array, f_eps1_array, f_energy_array):
                                                                                                                                     
    pylab.figure(1)
    pylab.plot (h_array,neg_S_array) 
    pylab.plot (h_array,f_eps0_array,'c')
    pylab.plot (h_array,f_eps1_array,'m')
    pylab.plot (h_array,f_energy_array,'r')
    pylab.show()  
    return

####################################################################################################
####################################################################################################
#
# Procedure to plot the analytic values versus the interaction enthalpy paramter h.
#
#    Inputs:    
#           h_array_analytic: the range of h values, where h =exp(2*epsilon1)
#           z1_array
#           z3_array
#           y2_array
#
####################################################################################################
####################################################################################################


def plot_analytic_vals (h_array, z1_array, z3_array, y2_array): 

    pylab.figure(2)
    pylab.plot (h_array,z1_array, 'g')    
    pylab.plot (h_array,z3_array,'r')
    pylab.plot (h_array,y2_array,'m')
    pylab.show()  
    return



####################################################################################################
####################################################################################################
#
# Function to randomly-generate an array, and then permute it to achieve the desired 
#    z1 & z3 values.
#
#    Inputs:    array_size_list: a list of two integers; array_length and layers
#               h: the interaction enthalpy parameter
#    Return: the matrix unit_array, a matrix of 0's and 1's.
#
####################################################################################################
####################################################################################################

# This function is currently not in use; it's also not completed ... 
def initialize_generated_matrix (array_size_list, h):

    array_length = array_size_list[0]
    array_layers = array_size_list[1]

    hSquared = h*h
    hFourth  = hSquared*hSquared
    denom    = 8.0*(hFourth - 6.0*hSquared + 1.0)  
    z3Analytic = (hSquared - 3.0)*(hSquared + 1.0)/denom      #  App. B, Eqn. 29
    z1Analytic = (1.0 - 3.0*hSquared)*(hSquared + 1.0)/denom  #  App. B, Eqn. 30
    y1Analytic = z1Analytic + 0.5*(0.5 - z1Analytic - z3Analytic)
    y2Analytic = z3Analytic + 0.5*(0.5 - z1Analytic - z3Analytic)

# Create the matrix 'unit_array' so that it has a random population of 0's and 1's.
    unit_array = np.random.choice([0, 1],size=(array_layers, array_length)) # Create an array filled with random values
# Note: this function can be used to create proportional distributions: np.random.choice([0, 1], size=(10,), p=[1./3, 2./3])


    print()
    print ('------------------------------------------------------------------------------')
    print ('------------------------------------------------------------------------------')
    print()
    print( 'With h = %.3f , we begin with a randomly-generated array:' % (h) )
    print()
    
# Bring the array closer to the desired configuration variable values    
    
    return unit_array

####################################################################################################
####################################################################################################
#
# Function to initialize the matrix with EITHER a pre-stored pattern of values
#    OR randomly-generate an array, and then permute it to achieve the desired 
#    z1 & z3 values.
#
#    Inputs:    array_size_list: a list of two integers; array_length and layers
#               pattern_select: an integer indicating whether to randomly-generate
#                   and then permute an array (0), or select a pattern (1 .. N)
#               h: the interaction enthalpy parameter
#    Return: the matrix unit_array, a matrix of 0's and 1's.
#
####################################################################################################
####################################################################################################

def initialize_matrix (array_size_list, pattern_select, h):

           
    array_length = array_size_list[0]
    array_layers = array_size_list[1]


# Note: The passed value patternProb is used to determine if we are returning a stored pattern, or
#       are probabilistically-generating our data.
#       If patternProb = 0: probabilistic generation, dependent on h 
#       If patternProb > 1: select one of the N stored patterns (1 ... N)
        
    if pattern_select == 0:
        unit_array = initialize_generated_matrix (array_size_list,h)


# Create the initial matrix, 'unit_array,' and populate it with zeros
    else:     
        unit_array = np.zeros((array_layers,array_length), dtype=np.int)      
                  
# Read the stored grid into unit_array
        x1_total = x2_total = 0   
        for i in range(0,array_layers):
            dataArray = obtainGridRow (i, pattern_select)
#            rownum = i+1
            for j in range(0, array_length):
                unit_array[i,j]=dataArray[j]
                if unit_array[i,j]==1: 
                    x1_total = x1_total + 1
                else: x2_total = x2_total + 1 
                
    print_grid (unit_array)
    
    
    return unit_array

####################################################################################################
####################################################################################################
#
# Procedure to  compute configuration variables x'i and return as elements of list configXVarsList
# (Yes, the x'i were computed during array creation and randomization. They are being recomputed 
#    as part of computing a list of ALL the configuration variables.)
#
####################################################################################################
####################################################################################################


def compute_config_X_variables (array_size_list, unit_array):


####################################################################################################
# This section unpacks the input variable array_size_list
####################################################################################################

    array_length = array_size_list [0]
    array_layers = array_size_list [1]
    unit_array = unit_array

# Debug print statements
    if not debug_print_off:
        print()  
        print( "Just entered compute_config_X_variables")
                                   
# Initialize the y'i variables
    x1_total = x2_total = 0 

    
    for i in range (0,array_layers):
  
        x1_partial = x2_partial = 0  

    # Compute the x'i values for each sub-row of the zigzag, just to see 
    #   the distribution 
    # Start counting through the array elements, L->R.
        for j in range(0, array_length):
            # If the initial unit is A:
            if unit_array[i,j]>0.1: 
                # The unit is "A," add it to x1 
                x1_partial = x1_partial + 1
            else: # The initial unit is B:
                x2_partial = x2_partial + 1
# debug prints
#        print "In row", i
#        print "x1_partial is", x1_partial, "x2_partial is", x2_partial
        x1_total = x1_total + x1_partial
        x2_total = x2_total + x2_partial
#        print "x1_total (so far) is", x1_total, "x2_total (so far) is", x2_total               
         
    x1 = x1_total
    x2 = x2_total       
    configVarsXList = (x1, x2)
        
#    print  "Leaving compute_config_X_variables for calling procedure"
#    print     
    return (configVarsXList)



####################################################################################################
####################################################################################################
#
# Procedure to compute the set of configuration variables y'i working across a single zigzag chain
# Procedure returns a list configvar containing the three y configuration variables:
#    y1 & y2 & y3
 
#
####################################################################################################
####################################################################################################



def computeConfigYEvenRowZigzagVariables (array_size_list, unit_array, topRow):


####################################################################################################
# This section unpacks the input variable array_size_list
####################################################################################################

    array_length = array_size_list [0]
    array_layers = array_size_list [1]
    unit_array = unit_array
                                 

###################################################################################################
#
# Compute the nearest-neighbor values y(i)
#
###################################################################################################


# y_1 is A-A
# y_3 is B-B
# left_y_2 is A-B
# right_y_2 is B-A
#
# The total number of y'i's is the same as the total number of x'i's.


# Initialize the y'i variables
    y1_total = left_y2_total = right_y2_total = y3_total = 0   

###################################################################################################
#
# Compute the nearest-neighbor values y(i) for the case of 
# downward-right-pointing diagonals, from top to next layer
# going L->R across the zigzag array
#
###################################################################################################

# Start counting through the layers; since we will work with a pair of 
# overlapping layers (for diagonal nearest-neighbors), we use a count of
# layers - 1. 

# commenting out for debug   
#for i in range(0,array_layers-1):
 #       top_row = i
  #      next_row = i+1
    top_row = topRow
    next_row = topRow + 1
  

# Start counting through the array elements, L->R.
    for j in range(0, array_length):
        # If the initial unit is A:
        if unit_array[top_row,j]>0.1: 
            # Compare with the same (jth) unit in the overlapping row 
            # comprising the zigzag chain
            # If the nearest-neighbor unit is also A:
            if unit_array[next_row,j] > 0.1:
                # h_increment the y_1; the count of A-A nearest-neighbor pairs:
                y1_total = y1_total + 1
            else: # The nearest-neighbor unit is B:
                left_y2_total = left_y2_total + 1
        else: # The initial unit is B:
            if unit_array[next_row,j] > 0.1:  # If the nearest-neighbor unit is A:
                right_y2_total = right_y2_total + 1            
            else: # The nearest-neighbor unit is also B:
                y3_total = y3_total + 1 
                
# Debug section: Print totals for right-downwards-pointing diagonals
#    print "Subtotals so far (downward-right-pointing-diagonals):"
#    print "(A-A) y1_total =", y1_total, "(A-B) left_y2_total =", left_y2_total   
#    print "(B-B) y3_total =", y3_total, "(B-A) right_y2_total =", right_y2_total 


###########################################################
#
# Compute the nearest-neighbor values y(i) for the case of 
# upward-right-pointing diagonals, from next-to-top layer up to  
# the top layer, going L->R across the zigzag array
#
###########################################################

# Recall that we are carrying forward previously-computed partial totals
# for the y'i values. 

# Start counting through the layers again, however, the computations will start
# with the lower layer and look in an upward-right-diagonal to the layer above. 

# commenting out for debug   
#for i in range(0,array_layers-1):
 #       top_row = i
  #      next_row = i+1
    

# Start counting through the array elements, L->R.
# Since we are comparing the unit in the lower row to the one shifted diagonally
# above and over to the right, we only step through to the array_length - 1 unit.
# A final step (after this) will be to compute the wrap-around. 
    for j in range(0, array_length-1):
        # If the initial unit is A:
        if unit_array[next_row,j]>0.1: 
            # Compare with the NEXT (j+1) unit in the overlapping top row 
            # comprising the zigzag chain
            # If the nearest-neighbor unit is also A:
            if unit_array[top_row,j+1] > 0.1:
                # h_increment the y_1; the count of A-A nearest-neighbor pairs:
                y1_total = y1_total + 1
            else: # The nearest-neighbor unit is B:
                left_y2_total = left_y2_total + 1
        else: # The initial unit is B:
            if unit_array[top_row,j+1] > 0.1:  # If the nearest-neighbor unit is A:
                right_y2_total = right_y2_total + 1            
            else: # The nearest-neighbor unit is also B:
                y3_total = y3_total + 1     

# Debug section: Print totals for right-upwards-pointing diagonals
#    print "Subtotals so far (downward + upward-right-pointing-diagonals):"
#    print "(A-A) y1_total =", y1_total, "(A-B) left_y2_total =", left_y2_total   
#    print "(B-B) y3_total =", y3_total, "(B-A) right_y2_total =", right_y2_total 

                
                                                
# Only one step remains.
# We need to compute the wrap-around for the zigzag chain (to get the total number
# of y'i's to be the same as the total number of x'i's. 
# We compute the nearest-neighbor pair similarity between the last unit on the 
# lower row with the first unit on the upper row.                                        

    if unit_array[next_row,array_length-1]>0.1: 
        # Compare with the FIRST unit in the overlapping top row 
        # comprising the zigzag chain
        # If the nearest-neighbor unit is also A:
        if unit_array[top_row,0] > 0.1:
        # Increment the y_1; the count of A-A nearest-neighbor pairs:
            y1_total = y1_total + 1
        else: # The nearest-neighbor unit is B:
            left_y2_total = left_y2_total + 1
    else: # The initial unit is B:
        if unit_array[top_row,0] > 0.1:  # If the nearest-neighbor unit is A:
            right_y2_total = right_y2_total + 1            
        else: # The nearest-neighbor unit is also B:
            y3_total = y3_total + 1 

#Debug section: Print message,"Computing last of the y'i values - wraparound"
#    print "Computing last of the y'i values - wraparound"

# This concludes computation of the y'i totals

    
################################################################
        
    if not debug_print_off:
        print()
        print( "Totals for the y'i variables:") 
        print( "(A-A) y1_total =", y1_total, "(A-B) left_y2_total =", left_y2_total )   
        print( "(B-B) y3_total =", y3_total, "(B-A) right_y2_total =", right_y2_total)  
        print()  

################################################################


###################################################################################################
#
# Assign the computed configuration variables to elements of the config_vars_list, 
# which will be passed back to the calling procedure
#
###################################################################################################
    
    y1 = y1_total
    y2 = left_y2_total + right_y2_total
    y3 = y3_total       
    configVarsYList = (y1, y2, y3)
      
    return (configVarsYList)



####################################################################################################
####################################################################################################
#
# Procedure to compute the set of configuration variables y'i
# Procedure returns a list configvar containing the three y configuration variables:
#    y1 & y2 & y3
 
#
####################################################################################################
####################################################################################################



def computeConfigYOddRowZigzagVariables (array_size_list, unit_array, topRow):


####################################################################################################
# This section unpacks the input variable array_size_list
####################################################################################################

    array_length = array_size_list [0]
    array_layers = array_size_list [1]
    unit_array = unit_array

# Initialize the y'i variables
    y1_total = left_y2_total = right_y2_total = y3_total = 0   

###################################################################################################
#
# Compute the nearest-neighbor values y(i) for the case of 
# downward-right-pointing diagonals, from top to next layer
# going L->R across the zigzag array
#
###################################################################################################

# Start counting through the layers; since we will work with a pair of 
# overlapping layers (for diagonal nearest-neighbors), we use a count of
# layers - 1. 

# commenting out for debug   
#for i in range(0,array_layers-1):
 #       top_row = i
  #      next_row = i+1
    top_row = topRow
    next_row = topRow + 1
    if top_row == array_layers-1: next_row = 0  
  

# Start counting through the array elements, L->R.
    for j in range(0, array_length-1): # Same logic as in the Even Row y(i) computation
           # but we go for one (TWO???) less down the array length 
        # If the initial unit is A:
        if unit_array[top_row,j]>0.1: 
            # Compare with the same (jth) unit in the overlapping row 
            # comprising the zigzag chain
            # If the nearest-neighbor unit is also A:
            if unit_array[next_row,j+1] > 0.1:
                # Increment the y_1; the count of A-A nearest-neighbor pairs:
                y1_total = y1_total + 1
            else: # The nearest-neighbor unit is B:
                left_y2_total = left_y2_total + 1
        else: # The initial unit is B:
            if unit_array[next_row,j+1] > 0.1:  # If the nearest-neighbor unit is A:
                right_y2_total = right_y2_total + 1            
            else: # The nearest-neighbor unit is also B:
                y3_total = y3_total + 1 
                
# Debug section: Print totals for right-downwards-pointing diagonals
#    print "Subtotals so far (downward-right-pointing-diagonals):"
#    print "(A-A) y1_total =", y1_total, "(A-B) left_y2_total =", left_y2_total   
#    print "(B-B) y3_total =", y3_total, "(B-A) right_y2_total =", right_y2_total 


###########################################################
#
# Compute the nearest-neighbor values y(i) for the case of 
# upward-right-pointing diagonals, from next-to-top layer up to  
# the top layer, going L->R across the zigzag array
#
###########################################################

# Recall that we are carrying forward previously-computed partial totals
# for the y'i values. 

# Start counting through the layers again, however, the computations will start
# with the lower layer and look in an upward-right-diagonal to the layer above. 

# commenting out for debug   
#for i in range(0,array_layers-1):
 #       top_row = i
  #      next_row = i+1
    

# Start counting through the array elements, L->R.
# Since we are comparing the unit in the lower row to the one shifted diagonally
# above and over to the right, we only step through to the array_length - 1 unit.
# A final step (after this) will be to compute the wrap-around. 
    for j in range(0, array_length): # Same logic as in the Even Row y(i) computation
            # But we can include the full array length (the other was truncated at array_length - 1)
        # If the initial unit is A:
        if unit_array[next_row,j]>0.1: 
            # Compare with the NEXT (j+1) unit in the overlapping top row 
            # comprising the zigzag chain
            # If the nearest-neighbor unit is also A:
            if unit_array[top_row,j] > 0.1:
                # Increment the y_1; the count of A-A nearest-neighbor pairs:
                y1_total = y1_total + 1
            else: # The nearest-neighbor unit is B:
                left_y2_total = left_y2_total + 1
        else: # The initial unit is B:
            if unit_array[top_row,j] > 0.1:  # If the nearest-neighbor unit is A:
                right_y2_total = right_y2_total + 1            
            else: # The nearest-neighbor unit is also B:
                y3_total = y3_total + 1     

# Debug section: Print totals for right-upwards-pointing diagonals
#    print "Subtotals so far (downward + upward-right-pointing-diagonals):"
#    print "(A-A) y1_total =", y1_total, "(A-B) left_y2_total =", left_y2_total   
#    print "(B-B) y3_total =", y3_total, "(B-A) right_y2_total =", right_y2_total 

                
                                                
# Only one step remains.
# We need to compute the wrap-around for the zigzag chain (to get the total number
# of y'i's to be the same as the total number of x'i's. 
# We compute the nearest-neighbor pair similarity between the last unit on the 
# lower row with the first unit on the upper row.                                        

    if unit_array[top_row,array_length-1]>0.1: 
        # Compare with the FIRST unit in the overlapping top row 
        # comprising the zigzag chain
        # If the nearest-neighbor unit is also A:
        if unit_array[next_row,0] > 0.1:
        # Increment the y_1; the count of A-A nearest-neighbor pairs:
            y1_total = y1_total + 1
        else: # The nearest-neighbor unit is B:
            left_y2_total = left_y2_total + 1
    else: # The initial unit is B:
        if unit_array[next_row,0] > 0.1:  # If the nearest-neighbor unit is A:
            right_y2_total = right_y2_total + 1            
        else: # The nearest-neighbor unit is also B:
            y3_total = y3_total + 1 

#Debug section: Print message,"Computing last of the y'i values - wraparound"
#    print "Computing last of the y'i values - wraparound"

# This concludes computation of the y'i totals


    
################################################################  

    if not debug_print_off:
        print()
        print( "Totals for the y'i variables:") 
        print( "(A-A) y1_total =", y1_total, "(A-B) left_y2_total =", left_y2_total )   
        print( "(B-B) y3_total =", y3_total, "(B-A) right_y2_total =", right_y2_total)  
        print()  

################################################################


###################################################################################################
#
# Assign the computed configuration variables to elements of the config_vars_list, 
# which will be passed back to the calling procedure
#
###################################################################################################
    
    y1 = y1_total
    y2 = left_y2_total + right_y2_total
    y3 = y3_total       
    configVarsYList = (y1, y2, y3)
      
    return (configVarsYList)



####################################################################################################
####################################################################################################

# This function runs both the even-to-odd and odd-to-even y(i) nearest neighbors; it combines the two in building
#   another row on top of the basic 1-D zigzag chain

####################################################################################################

def compute_config_Y_variables (array_size_list, unit_array):

# Initialize the y'i variables

    y1 = y2 = y3 = 0

    if not debug_print_off:
        print()
        print( '  Starting to compute Y variables')
        print( '  Total number of pairs of zigzag chains is: ', pairs)
        print() 
    for i in range (0, pairs): 
        topRow = 2*i
        if not debug_print_off:
            print( '  Row: ', topRow)
        # Obtain the y(i) values from the first even-to-odd zigzag chain (0 to 1, running top-to-bottom)
        configVarsYList = computeConfigYEvenRowZigzagVariables (array_size_list, unit_array, topRow) 
        # Assign the returned results to the local sum for each of the z(i) triplets
        y1 = y1+configVarsYList[0]
        y2 = y2+configVarsYList[1]
        y3 = y3+configVarsYList[2]
        topRow = 2*i+1
        if not debug_print_off:
            print()
            print( '  Row: ', topRow)
        configVarsYList = computeConfigYOddRowZigzagVariables (array_size_list, unit_array, topRow) 
        # Assign the returned results to the local sum for each of the z(i) triplets
        y1 = y1+configVarsYList[0]
        y2 = y2+configVarsYList[1]
        y3 = y3+configVarsYList[2]



# Debug section: Print totals for right-downwards-then-upwards triplets
        if not debug_print_off:
            print()
            print( ' -----------')
            print()
            print( "Totals for all y(i), after completing Row: ", topRow)
            print( "           (A-A) y1_total =", y1)
            print( "(A-B) plus (B-A) y2_total =", y2)
            print( "           (B-B) y3_total =", y3 )  
            print()
            print( ' -----------')
            print()
                

    configVarsYList = (y1, y2, y3)                                                                                                                                                                        
    return (configVarsYList)





####################################################################################################
####################################################################################################
#
# Procedure to compute horizontal configuration variables w'i and return as elements of list configXVarsList
#
####################################################################################################
####################################################################################################

def computeConfigWHorizontalRowVariables (array_size_list, unit_array):


####################################################################################################
# This section unpacks the input variable array_size_list
####################################################################################################

    array_length = array_size_list [0]
    array_layers = array_size_list [1]
    unit_array = unit_array

# Debug print statements
    if not debug_print_off:
         print( "Just entered computeConfigWVariables")
                                   
# Initialize the w'i variables
    w1_total = w2_total = w3_total = 0 
    w1_partial = w2_partial = w3_partial = 0  
    
    for i in range (0,array_layers):
        w1_partial = w2_partial = w3_partial = 0  
# Compute the w'i values for each sub-row of the zigzag, just to see 
#   the distribution 
# Start counting through the array elements, L->R.
        for j in range(0, array_length):             
            nextNearestNeighbor = j+1
            rowLimit = array_length-1
            if j == rowLimit: nextNearestNeighbor = 0
            # If the initial unit is A:
            if unit_array[i,j]>0.1: 
                # The unit is "A," see if the next unit is "A" or "B" 
                if unit_array[i,nextNearestNeighbor]>0.1: 
                # Compare with the NEXT (j+1) unit in the SAME row 
                #   comprising a partial row of the the zigzag chain
                #   If this next-nearest-neighbor unit is also "A":                
                    w1_partial = w1_partial + 1
                else: # The next-nearest-neighbor is in "B"
                    w2_partial = w2_partial + 1                    
                                    
            else: # The initial unit is B:
                # The unit is "B," see if the next unit is "A" or "B" 
                if unit_array[i,nextNearestNeighbor]>0.1: 
                # Compare with the NEXT (j+1) unit in the SAME row 
                #   comprising a partial row of the the zigzag chain
                #   If this next-nearest-neighbor unit is also "A":                
                    w2_partial = w2_partial + 1
                else: # The next-nearest-neighbor is in "B"
                    w3_partial = w3_partial + 1                    

        detailed_debug_print_offW = True
        if not detailed_debug_print_offW:
            print()
            print( "In row ", i)
            print( "w1_partial = ", w1_partial)
            print( "w2_partial = ", w2_partial)                 
            print( "w3_partial = ", w3_partial) 
                                                        
        # Check the wrap-around value between the last unit in the row
        #   and the first item of this same row
#        if unit_array[i,array_length-1]>0.1: 
#            # The unit is "A," see if the wraparound unit is "A" or "B" 
#            if unit_array[i,0] > 0.1: #This unit is "A"
#                w1_partial = w1_partial + 1
#            else: w2_partial = w2_partial + 1                    
#        else: 
#            if unit_array[i,0] > 0.1: #This unit is "A"
#                w2_partial = w2_partial + 1
#            else: w3_partial = w3_partial + 1                    
                
                    
#        print "In row", i, "after wrap-around - still testing for A"
#        print "w1_partial = ", w1_partial
#        print "w2_partial = ", w2_partial

        w1_total = w1_total + w1_partial
        w2_total = w2_total + w2_partial 
        w3_total = w3_total + w3_partial                                     
                                                                                                            
    w1 = w1_total
    w2 = w2_total 
    w3 = w3_total      
    configVarsWList = (w1, w2, w3)

################################################################
    if not debug_print_off:
        print()
        print( "Totals for all horizontal w(i)")
        print( "            (A--A) w1_total =", w1)
        print( "(A--B) plus (B--A) w2_total =", w2)
        print( "            (B--B) w3_total =", w3)   
        print()          

################################################################
                                   
                                                                      
                                                                                                                                            
    return (configVarsWList)


####################################################################################################
####################################################################################################
#
# Procedure to compute vertical configuration variables w'i and return as elements of list configWVarsList
#
####################################################################################################
####################################################################################################

def computeConfigWVerticalColVariables (array_size_list, unit_array):


####################################################################################################
# This section unpacks the input variable array_size_list
####################################################################################################

    array_length = array_size_list [0]
    array_layers = array_size_list [1]
    iLimit = array_layers - 2
    
# Debug print statements
    if not debug_print_off:
        print( "Just entered computeConfigWVerticalVariables")
                                   
# Initialize the w'i variables
    w1_total = w2_total = w3_total = 0 
    w1_partial = w2_partial = w3_partial = 0  

#    

    for i in range (0, array_layers): # run through the rows, look at those two rows apart
        vertPair = i+2        
        if i == iLimit: vertPair = 0
        if i == iLimit+1: vertPair = 1
        for j in range(0, array_length):             
            # If the initial unit is A:
            if unit_array[i,j]>0.1: 
                # The unit is "A," see if the next unit is "A" or "B" 
                if unit_array[vertPair,j]>0.1: 
                # Compare with the NEXT (j+1) unit in the SAME row 
                #   comprising a partial row of the the zigzag chain
                #   If this next-nearest-neighbor unit is also "A":                
                    w1_partial = w1_partial + 1
                else: # The next-nearest-neighbor is in "B"
                    w2_partial = w2_partial + 1                    
                                    
            else: # The initial unit is B:
                # The unit is "B," see if the next unit is "A" or "B" 
                if unit_array[vertPair,j]>0.1: 
                # Compare with the NEXT (j+1) unit in the SAME row 
                #   comprising a partial row of the the zigzag chain
                #   If this next-nearest-neighbor unit is also "A":                
                    w2_partial = w2_partial + 1
                else: # The next-nearest-neighbor is in "B"
                    w3_partial = w3_partial + 1     

           
#                        

    w1_total = w1_total + w1_partial
    w2_total = w2_total + w2_partial 
    w3_total = w3_total + w3_partial                                     
                                                                                                            
    w1 = w1_total
    w2 = w2_total 
    w3 = w3_total      
    configVarsWList = (w1, w2, w3)

      

    
################################################################
    if not debug_print_off:
        print()
        print( "Totals for all horizontal w(i)")
        print( "            (A--A) w1_total =", w1)
        print( "(A--B) plus (B--A) w2_total =", w2)
        print( "            (B--B) w3_total =", w3)   
        print()           

################################################################
            
                        
    return (configVarsWList)






####################################################################################################
####################################################################################################

# This function runs both the even-to-odd and odd-to-even zigzags; it is the first step in building
#   another row on top of the basic 1-D zigzag chain

####################################################################################################

def computeConfigWVariables (array_size_list, unit_array):

# Initialize the y'i variables

    w1 = w2 = w3 = 0

    if not debug_print_off:
        print()
        print( '  Starting to compute W variables')
        print( '  Total number of zigzag chains is: ', array_layers)
        print()  
    
    configVarsWList = computeConfigWHorizontalRowVariables (array_size_list, unit_array) 
    # Assign the returned results to the local sum for each of the z(i) triplets
    w1 = w1+configVarsWList[0]
    w2 = w2+configVarsWList[1]
    w3 = w3+configVarsWList[2]


# Debug section: Print totals for right-downwards-then-upwards triplets
    if not debug_print_off:
        print()
        print('  Row: ', array_layers)
        print()
        print( ' -----------')
        print( ' ')
        print( "Totals for all horizontal w(i), after completing Row: ", array_layers)
        print( "            (A--A) w1_total =", w1)
        print( "(A--B) plus (B--A) w2_total =", w2)
        print( "            (B--B) w3_total =", w3)  
        print()
        print( ' -----------')
        print()

  
# NOTE: Still need to write the computation for an extra odd row in grid, if it exists


    configVarsWList = computeConfigWVerticalColVariables (array_size_list, unit_array)
    w1 = w1+configVarsWList[0]
    w2 = w2+configVarsWList[1]
    w3 = w3+configVarsWList[2]       

# Debug section: Print totals for right-downwards-then-upwards triplets
    if not debug_print_off: 
        print()
        print( ' -----------')
        print( ' ')
        print( "Totals for all horizontal and vertical w(i)")
        print( "            (A--A) w1_total =", w1)
        print( "(A--B) plus (B--A) w2_total =", w2)
        print( "            (B--B) w3_total =", w3)  
        print()
        print( ' -----------')
        print()

    configVarsWList = (w1, w2, w3)                                                                                                                                                                        
    return (configVarsWList)





####################################################################################################
####################################################################################################
#
# Procedure to compute the the precise value of a triplet z'i variable given
#   Input: integer values for unit (U), nearest-neighbor (NN), next-nearest-neighbor (= NNN)
#   Output: locally-h_incremented values for z1 & z2 & z3 & z4 & z5 & z6 
#
####################################################################################################
####################################################################################################



def computeSpecificTripletZVariable (U, NN, NNN):

# Debug print statements
    if not detailed_debug_print_off:
        print()
        print( "Just entered computeSpecificTripletZVariable")
    
    localTripletValueList = list()
    
    z1 = 0
    left_z2 = 0
    right_z2 = 0
    z3 = 0
    z4 = 0
    left_z5 = 0
    right_z5 = 0
    z6 = 0    

    if U > 0.1: 
        # Compare with the same (jth) unit in the overlapping row 
        # comprising the zigzag chain
        # If the nearest-neighbor unit is also A:
        if NN > 0.1:
            # We have the first portion of A-A-X triplet:
            if NNN > 0.1: 
                # We have an A-A-A triplet
                z1 = 1
            else: 
                # We have an A-A-B triplet
                left_z2 = 1
        else: # The nearest-neighbor unit is B, we have an A-B-X triplet:
            if NNN > 0.1: 
                # We have an A-B-A triplet
                z3 = 1
            else: 
                # We have an A-B-B triplet
                right_z5 = 1
    else: # The initial unit is B:
        if NN  > 0.1:  # If the nearest-neighbor unit is A:
            # We have the first portion of B-A-X triplet:                    
            if NNN > 0.1: 
                # We have an B-A-A triplet
                right_z2 = 1
            else: 
                # We have an B-A-B triplet                        
                z4 = 1          
        else: # The nearest-neighbor unit is also B:
            # We have the first portion of B-B-X triplet:                     
            if NNN >0.1: 
                # We have an B-B-A triplet
                left_z5 = 1
            else: 
                # We have an B-B-B triplet                        
                z6 = 1 
           
                                 
    localTripletValueList = (z1, left_z2, right_z2, z3, z4, left_z5, right_z5, z6)
      
    return (localTripletValueList)

###################################################################################################
#
# Procedure to debug print the newly computed triplet values
#
###################################################################################################

def debugPrintZ (z1_h_incr, left_z2_h_incr, right_z2_h_incr, z3_h_incr, z4_h_incr, left_z5_h_incr, right_z5_h_incr, z6_h_incr):
        
    print()
    print( '(A-A-A) z1_h_incr =', z1_h_incr, '(A-A-B) left_z2_h_incr =', left_z2_h_incr )
    print( '(A-B-A) z3_h_incr =', z3_h_incr, '(B-A-A) right_z2_h_incr =', right_z2_h_incr )
    print( '(B-A-B) z4_h_incr =', z4_h_incr, '(B-B-A) left_z5_h_incr =', left_z5_h_incr )
    print( '(B-B-B) z6_h_incr =', z6_h_incr, '(A-B-B) right_z5_h_incr =', right_z5_h_incr ) 
    print()    
    return
    

####################################################################################################
####################################################################################################
#
# Function to compute the set of configuration variables z'i going upper-to-lower across two rows
#    starting with an EVEN row (0 to 1, 2 to 3, etc.) 
# Function returns a list configvar containing the six z configuration variables:
#    z1 & z2 & z3 & z4 & z5 & z6 
#
####################################################################################################
####################################################################################################



def computeConfigZEvenUpperToLower (array_size_list, unit_array, top_row):

    array_length = array_size_list [0]
    array_layers = array_size_list [1]
    unit_array = unit_array
    
# Create the array to hold the partial (the h_increments in the) z'i's, and populate it with zeros
    zPartialArray = np.zeros((array_length), dtype=np.int)

  
    z1_partial = left_z2_partial = right_z2_partial = z3_partial = 0
    z4_partial = left_z5_partial = right_z5_partial = z6_partial = 0


    next_row = top_row + 1

# Start counting through the array elements, L->R.
    for j in range(0, array_length-1):
        U = unit_array[top_row,j]
        NN = unit_array[next_row,j]
        NNN = unit_array[top_row, j+1]

        TripletValueList = computeSpecificTripletZVariable (U, NN, NNN)

    

# Debug print statements
#        if not detailed_debug_print_off:
#            print ' '
#            print 'Debug printing: computeConfigZEvenUpperToLower'  #debug_print_off false        
#            print "Returning from computeSpecificTripletZVariable" 
#            print "Unpacking the specific triplet value found"
    
        z1_h_incr = TripletValueList[0]
        left_z2_h_incr = TripletValueList[1]
        right_z2_h_incr = TripletValueList[2]
        z3_h_incr = TripletValueList[3]
        z4_h_incr = TripletValueList[4] 
        left_z5_h_incr = TripletValueList[5]
        right_z5_h_incr = TripletValueList[6] 
        z6_h_incr = TripletValueList[7]
            
        # Debug print statements
        if not detailed_debug_print_off:
            print()
            print( 'Debug printing: computeConfigZEvenUpperToLower, h_incrementing the z(i) in for loop:')  #debug_print_off false
            debugPrintZ (z1_h_incr, left_z2_h_incr, right_z2_h_incr, z3_h_incr, z4_h_incr, left_z5_h_incr, right_z5_h_incr, z6_h_incr)
    
                                    
        z1_partial = z1_partial + z1_h_incr 
        left_z2_partial = left_z2_partial + left_z2_h_incr
        right_z2_partial = right_z2_partial + right_z2_h_incr 
        z3_partial = z3_partial + z3_h_incr 
        z4_partial = z4_partial + z4_h_incr 
        left_z5_partial = left_z5_partial + left_z5_h_incr 
        right_z5_partial = right_z5_partial + right_z5_h_incr 
        z6_partial = z6_partial + z6_h_incr 
      
        # Completed for loop; have gone through entire two-row zigzag with triplets that are upper-to-lower-then-upper; 
        #   no wrap-arounds
        
    zPartialArray[0] = z1_partial
    zPartialArray[1] = left_z2_partial
    zPartialArray[2] = right_z2_partial
    zPartialArray[3] = z3_partial        
    zPartialArray[4] = z4_partial  
    zPartialArray[5] = left_z5_partial
    zPartialArray[6] = right_z5_partial    
    zPartialArray[7] = z6_partial    
    return (zPartialArray)    
            

####################################################################################################
####################################################################################################
#
# Function to compute the set of configuration variables z'i going lower-to-upper across two rows
#    as part of calculations for SET beginning with an EVEN row (0 to 1, 2 to 3, etc.) 
# Function returns a list configvar containing the six z configuration variables:
#    z1 & z2 & z3 & z4 & z5 & z6 
#
####################################################################################################
####################################################################################################



def computeConfigZEvenLowerToUpper (array_size_list, unit_array, top_row):

    array_length = array_size_list [0]
    array_layers = array_size_list [1]
    unit_array = unit_array
    
# Create the array to hold the partial (the h_increments in the) z'i's, and populate it with zeros
    zPartialArray = np.zeros((array_length), dtype=np.int)

  
    z1_partial = left_z2_partial = right_z2_partial = z3_partial = 0
    z4_partial = left_z5_partial = right_z5_partial = z6_partial = 0

    next_row = top_row + 1


# NOTE: We are computing the SECOND row of triplets in a zigzag chain,
#   going from the lower row to the top
# Start counting through the array elements, L->R.
    for j in range(0, array_length-1):
        U = unit_array[next_row,j]
        NN = unit_array[top_row,j+1]
        NNN = unit_array[next_row, j+1]

        TripletValueList = computeSpecificTripletZVariable (U, NN, NNN)

    
    # Debug print statements
        if not detailed_debug_print_off:
            print()
            print( "Computing the SECOND row of a zigzag chain for leading unit ", j)
            print( "Returning from computeSpecificTripletZVariable") 
            print( "Unpacking the specific triplet value found")
    
        z1_h_incr = TripletValueList[0]
        left_z2_h_incr = TripletValueList[1]
        right_z2_h_incr = TripletValueList[2]
        z3_h_incr = TripletValueList[3]
        z4_h_incr = TripletValueList[4] 
        left_z5_h_incr = TripletValueList[5]
        right_z5_h_incr = TripletValueList[6] 
        z6_h_incr = TripletValueList[7]

        # Debug print statements
        if not detailed_debug_print_off:
            print()
            print( 'Debug printing: computeConfigZEvenLowerToUpper, h_incrementing the z(i) in for loop:')  #debug_print_off false
            debugPrintZ (z1_h_incr, left_z2_h_incr, right_z2_h_incr, z3_h_incr, z4_h_incr, left_z5_h_incr, right_z5_h_incr, z6_h_incr)
    
                                    
        z1_partial = z1_partial + z1_h_incr 
        left_z2_partial = left_z2_partial + left_z2_h_incr
        right_z2_partial = right_z2_partial + right_z2_h_incr 
        z3_partial = z3_partial + z3_h_incr 
        z4_partial = z4_partial + z4_h_incr 
        left_z5_partial = left_z5_partial + left_z5_h_incr 
        right_z5_partial = right_z5_partial + right_z5_h_incr 
        z6_partial = z6_partial + z6_h_incr 

        # Completed for loop; have gone through entire two-row zigzag with triplets that are lower-to-upper-to-lower; 
        #   no wrap-arounds
        

    zPartialArray[0] = z1_partial
    zPartialArray[1] = left_z2_partial
    zPartialArray[2] = right_z2_partial
    zPartialArray[3] = z3_partial        
    zPartialArray[4] = z4_partial  
    zPartialArray[5] = left_z5_partial
    zPartialArray[6] = right_z5_partial    
    zPartialArray[7] = z6_partial    
    return (zPartialArray)  



# ************************************************************************************************ #
#
# The following are a collection of print functions
#
# ************************************************************************************************ #

############################################
#
# print function: two rows; EVEN-to-ODD
# 
#-------------------------------------------

def printEvenToOddRows (top_row, unit_array):

    next_row = top_row + 1
    print( ' *************************')
    print()     
    print( 'top_row = ', top_row, ' next_row = ', next_row)
    print( 'Row', top_row, ':', blnkspc, )
    for j in range(0,array_length):
        print( unit_array[top_row,j], blnkspc,)
    print() 

    print( 'Row ', next_row, ':', blnkspc,)
    print (blnkspc,)
    for j in range(0,array_length):
        print( unit_array[next_row,j], blnkspc,)
    print( )
    print( ' *************************')
    return

####################################################################################################
####################################################################################################
#
# Function to compute the set of configuration variables z(i), going from EVEN-to-ODD rows
# Function returns a list configvar containing the six z configuration variables:
#    z1 & z2 & z3 & z4 & z5 & z6 
#
####################################################################################################
####################################################################################################



def compute_config_Z_variablesEvenToOdd (array_size_list, unit_array, topRow):


####################################################################################################
# This section unpacks the input variable array_size_list
####################################################################################################

    array_length = array_size_list [0]
    array_layers = array_size_list [1]
    unit_array = unit_array

# Debug print statements
    if not detailed_debug_print_off:
        print()
        print( "Just entered compute_config_Z_variables: Even-to-Odd")
                                  
                                                               

###################################################################################################
#
# Compute the triplet values z(i)
#
###################################################################################################


# z_1 is A-A-A
# z_3 is A-B-A 
# left_z_2 is A-A-B
# right_z_2 is B-A-A

# z_4 is B-A-B
# z_6 is B-B-B 
# left_z_5 is B-B-A
# right_z_5 is A-B-B

#
# The total number of z'i's is 6.


# Initialize the z'i variables

    z1 = z2 = z3 = z4 = z5 = z6 = 0

    z1_total = left_z2_total = right_z2_total = z3_total = 0   
    z4_total = left_z5_total = right_z5_total = z6_total = 0   
    y1_total = left_y2_total = right_y2_total = y3_total = 0          

###################################################################################################
#
# For an EVEN-to-ODD row combination (0 & 1, 2 & 3, etc): 
# Compute the triplet values z(i) for the case of 
# downward-right-then-upwards-right-pointing triplets, from top to next layer down
# going L->R across the zigzag array
# This step does not include any wrap-arounds
#
###################################################################################################

# Start counting through the layers; since we will work with a pair of 
# overlapping layers (for diagonal nearest-neighbors), we use a count of
# layers - 1. 

# Debug print statements
    if not detailed_debug_print_off:
        print()
        print( "Calling computeSpecificTripletZVariable")

    U  = NN = NNN = 0
    
    TripletValueList = list()  
 
# commenting out for debug   
#for i in range(0,array_layers-1):
 #       top_row = i
 #      next_row = i+1
    top_row = topRow
    next_row = topRow + 1  

    if not z_debug_print_off: 
        printEvenToOddRows (top_row, unit_array)       
      
    z1_partial = left_z2_partial = right_z2_partial = z3_partial = 0
    z4_partial = left_z5_partial = right_z5_partial = z6_partial = 0

# Compute the first contribution to z'i's from the top-to-bottom-to-top row
    zPartialArray = computeConfigZEvenUpperToLower (array_size_list, unit_array, top_row) 
             
# Unpack the new z(i) contributions into the partial values for z(i)                                     
    z1_partial= zPartialArray[0]
    left_z2_partial = zPartialArray[1]
    right_z2_partial = zPartialArray[2]  
    z3_partial = zPartialArray[3]        
    z4_partial = zPartialArray[4]   
    left_z5_partial = zPartialArray[5] 
    right_z5_partial = zPartialArray[6]     
    z6_partial = zPartialArray[7]                            
                                                    
# Update the total z'i values:                                                                                                                                                                 
    z1_total = z1_total + z1_partial
    left_z2_total = left_z2_total + left_z2_partial 
    right_z2_total = right_z2_total + right_z2_partial      
    z3_total = z3_total + z3_partial
    z4_total = z4_total + z4_partial
    left_z5_total = left_z5_total + left_z5_partial
    right_z5_total = right_z5_total + right_z5_partial  
    z6_total = z6_total + z6_partial
    z5_total = left_z5_total + right_z5_total
    z2_total = left_z2_total + right_z2_total
    
# Debug section: Print totals for right-downwards-then-upwards triplets
    if not z_debug_print_off:
        print()
        print( "Subtotals so far (downward-right-then_upwards-right-pointing triplets)")
        print( "Before any wrap-arounds:")
        print( "(A-A-A) z1_total =", z1_total, "(A-A-B) left_z2_total =", left_z2_total )
        print( "(A-B-A) z3_total =", z3_total, "(B-A-A) right_z2_total =", right_z2_total ) 
        print( "(B-A-B) z4_total =", z4_total, "(B-B-A) left_z5_total =", left_z5_total )
        print( "(B-B-B) z6_total =", z6_total, "(A-B-B) right_z5_total =", right_z5_total)



###################################################################################################
#
# Compute the triplet values z(i) for the first wrap-around
# Start at last unit on top row, use downward-pointing-diagonal for last unit on next row
# Then wrap-around to pick up first unit on top row.
#
###################################################################################################

  

# Debug print statements
    if not detailed_debug_print_off: 
        print()
        print( "Computing first wrap-around triplet")
        print( "Calling computeSpecificTripletZVariable")

    U  = NN = NNN = 0

    U = unit_array[top_row,array_length-1]
    NN = unit_array[next_row,array_length-1]
    NNN = unit_array[top_row, 0]

    TripletValueList = computeSpecificTripletZVariable (U, NN, NNN)

    
    # Debug print statements
    if not detailed_debug_print_off:  
        print()
        print( "Returning from computeSpecificTripletZVariable" )
        print( "Unpacking the specific triplet value found")
    
    z1_h_incr = TripletValueList[0]
    left_z2_h_incr = TripletValueList[1]
    right_z2_h_incr = TripletValueList[2]
    z3_h_incr = TripletValueList[3]
    z4_h_incr = TripletValueList[4] 
    left_z5_h_incr = TripletValueList[5]
    right_z5_h_incr = TripletValueList[6] 
    z6_h_incr = TripletValueList[7]

    # Debug print statements

    if not debug_print_off: 
        print()         
        print( "(A-A-A) z1_h_incr =", z1_h_incr, "(A-A-B) left_z2_h_incr =", left_z2_h_incr )
        print( "(A-B-A) z3_h_incr =", z3_h_incr, "(B-A-A) right_z2_h_incr =", right_z2_h_incr ) 
        print( "(B-A-B) z4_h_incr =", z4_h_incr, "(B-B-A) left_z5_h_incr =", left_z5_h_incr )
        print( "(B-B-B) z6_h_incr =", z6_h_incr, "(A-B-B) right_z5_h_incr =", right_z5_h_incr ) 
        print( ' End of h_incremental z(i) for first wrap-around, even-to-odd, top row = ', top_row)
        print()                                   
                                                
# Update the total z'i values:                                                                                                                                                                 
    z1_total = z1_total + z1_h_incr
    left_z2_total = left_z2_total + left_z2_h_incr 
    right_z2_total = right_z2_total + right_z2_h_incr      
    z3_total = z3_total + z3_h_incr
    z4_total = z4_total + z4_h_incr
    left_z5_total = left_z5_total + left_z5_h_incr
    right_z5_total = right_z5_total + right_z5_h_incr  
    z6_total = z6_total + z6_h_incr
    z5_total = left_z5_total + right_z5_total
    z2_total = left_z2_total + right_z2_total
    
# Debug section: Print totals for right-downwards-then-upwards triplets
    if not z_debug_print_off:
        print()
        print( "Subtotals so far (downward-right-then_upwards-right-pointing triplets)")
        print(  "After adding in the top-layer wrap-around:")
        print(  "(A-A-A) z1_total =", z1_total, "(A-A-B) left_z2_total =", left_z2_total )
        print(  "(A-B-A) z3_total =", z3_total, "(B-A-A) right_z2_total =", right_z2_total  )
        print(  "(B-A-B) z4_total =", z4_total, "(B-B-A) left_z5_total =", left_z5_total )
        print(  "(B-B-B) z6_total =", z6_total, "(A-B-B) right_z5_total =", right_z5_total)


###########################################################
#
# Compute the triplet values z(i) for the case of 
# upward-right-pointing diagonals, from next-to-top layer up to  
# the top layer, then going down again, 
# going L->R across the zigzag array
#
###########################################################

# Recall that we are carrying forward previously-computed totals
# for the z'i values. 
# However, we will re-initialize the partial totals
# so that we get a partial total across each layer of the zigzag

# Start counting through the layers again, however, the computations will start
# with the lower layer and look in an upward-right-diagonal to the layer above. 

# commenting out for debug   
#for i in range(0,array_layers-1):
 #       top_row = i
  #      next_row = i+1
     
#    for i in range (0,1):
#        top_row = 0
#        next_row = 1  
  
    zPartialArray = computeConfigZEvenLowerToUpper (array_size_list, unit_array, top_row)
             
# Unpack the new z(i) contributions into the partial values for z(i)                                     
    z1_partial= zPartialArray[0]
    left_z2_partial = zPartialArray[1]
    right_z2_partial = zPartialArray[2]  
    z3_partial = zPartialArray[3]        
    z4_partial = zPartialArray[4]   
    left_z5_partial = zPartialArray[5] 
    right_z5_partial = zPartialArray[6]     
    z6_partial = zPartialArray[7]                            
   
             
# Update the total z'i values:                                                                                                                                                                 
    z1_total = z1_total + z1_partial
    left_z2_total = left_z2_total + left_z2_partial 
    right_z2_total = right_z2_total + right_z2_partial      
    z3_total = z3_total + z3_partial
    z4_total = z4_total + z4_partial
    left_z5_total = left_z5_total + left_z5_partial
    right_z5_total = right_z5_total + right_z5_partial  
    z6_total = z6_total + z6_partial
    z5_total = left_z5_total + right_z5_total
    z2_total = left_z2_total + right_z2_total
    
# Debug section: Print totals for right-downwards-then-upwards triplets
    if not z_debug_print_off:
        print()
        print( "Subtotals so far (adding in upwards-right-then_downwards-right-pointing triplets)")
        print( "Before the last wrap-around:")
        print( "(A-A-A) z1_total =", z1_total, "(A-A-B) left_z2_total =", left_z2_total )
        print( "(A-B-A) z3_total =", z3_total, "(B-A-A) right_z2_total =", right_z2_total ) 
        print( "(B-A-B) z4_total =", z4_total, "(B-B-A) left_z5_total =", left_z5_total )
        print( "(B-B-B) z6_total =", z6_total, "(A-B-B) right_z5_total =", right_z5_total)


                 
                                                
# Only one step remains.
# We need to compute the second wrap-around for the zigzag chain (to get the total number
# of z'i's . 
                                       


###################################################################################################
#
# Compute the triplet values z(i) for the second wrap-around
# Start at last unit on bottom row, use upward-pointing-diagonal for first unit on upper row 
# Then wrap-around to pick up first unit on bottom row.
#
###################################################################################################

  

# Debug print statements
    if not z_debug_print_off:
        print( "Computing second wrap-around triplet")
        print( "Calling computeSpecificTripletZVariable")
    U  = NN = NNN = 0

    U = unit_array[next_row,array_length-1]
    NN = unit_array[top_row,0]
    NNN = unit_array[next_row, 0]

    TripletValueList = computeSpecificTripletZVariable (U, NN, NNN)

    
    # Debug print statements
#    print "Returning from computeSpecificTripletZVariable" 
#    print "Unpacking the specific triplet value found"
    
    z1_h_incr = TripletValueList[0]
    left_z2_h_incr = TripletValueList[1]
    right_z2_h_incr = TripletValueList[2]
    z3_h_incr = TripletValueList[3]
    z4_h_incr = TripletValueList[4] 
    left_z5_h_incr = TripletValueList[5]
    right_z5_h_incr = TripletValueList[6] 
    z6_h_incr = TripletValueList[7]

    # Debug print statements
    if not z_debug_print_off:
        print()          
        print( "(A-A-A) z1_h_incr =", z1_h_incr, "(A-A-B) left_z2_h_incr =", left_z2_h_incr )
        print( "(A-B-A) z3_h_incr =", z3_h_incr, "(B-A-A) right_z2_h_incr =", right_z2_h_incr ) 
        print( "(B-A-B) z4_h_incr =", z4_h_incr, "(B-B-A) left_z5_h_incr =", left_z5_h_incr )
        print( "(B-B-B) z6_h_incr =", z6_h_incr, "(A-B-B) right_z5_h_incr =", right_z5_h_incr ) 
        print() 
                                    
                                                
# Update the total z'i values:                                                                                                                                                                 
    z1_total = z1_total + z1_h_incr
    left_z2_total = left_z2_total + left_z2_h_incr 
    right_z2_total = right_z2_total + right_z2_h_incr      
    z3_total = z3_total + z3_h_incr
    z4_total = z4_total + z4_h_incr
    left_z5_total = left_z5_total + left_z5_h_incr
    right_z5_total = right_z5_total + right_z5_h_incr  
    z6_total = z6_total + z6_h_incr
    z5_total = left_z5_total + right_z5_total
    z2_total = left_z2_total + right_z2_total
    
# Debug section: Print totals for right-downwards-then-upwards triplets
    if not z_debug_print_off:
        print( "Totals for all triplets, after adding in the second wrap-around")
        print( "(A-A-A) z1_total =", z1_total, "(A-A-B) left_z2_total =", left_z2_total )
        print( "(A-B-A) z3_total =", z3_total, "(B-A-A) right_z2_total =", right_z2_total ) 
        print( "(B-A-B) z4_total =", z4_total, "(B-B-A) left_z5_total =", left_z5_total )
        print( "(B-B-B) z6_total =", z6_total, "(A-B-B) right_z5_total =", right_z5_total )



# This concludes computation of the z'i totals for a complete pass through a zigzag chain


###################################################################################################
#
# Assign the computed configuration variables to elements of the config_vars_list, 
# which will be passed back to the calling procedure
#
###################################################################################################
    
    z1 = z1_total
    z2 = left_z2_total + right_z2_total
    z3 = z3_total 
    z4 = z4_total
    z5 = left_z5_total + right_z5_total
    z6 = z6_total           
    configVarsZList = (z1, z2, z3, z4, z5, z6)
      
    return (configVarsZList)



#**************************************************************************************************

####################################################################################################
####################################################################################################
#
# FUNCTION AREA to compute the set of configuration variables z(i) for ODD-to-EVEN rows
# Function returns a list configvar containing the six z configuration variables:
#    z1 & z2 & z3 & z4 & z5 & z6 
#
####################################################################################################
####################################################################################################


####################################################################################################
####################################################################################################
#
# Function to compute the set of configuration variables z'i going upper-to-lower across two rows
#    starting with an ODD row (0 to 1, 2 to 3, etc.) 
# Function returns a list configvar containing the six z configuration variables:
#    z1 & z2 & z3 & z4 & z5 & z6 
#
####################################################################################################
####################################################################################################



def computeConfigZOddUpperToLower (array_size_list, unit_array, top_row):

    array_length = array_size_list [0]
    array_layers = array_size_list [1]
    unit_array = unit_array
    
# Create the array to hold the partial (the h_increments in the) z'i's, and populate it with zeros
    zPartialArray = np.zeros((array_length), dtype=np.int)

  
    z1_partial = left_z2_partial = right_z2_partial = z3_partial = 0
    z4_partial = left_z5_partial = right_z5_partial = z6_partial = 0


    next_row = top_row + 1
    bottomRow = array_layers
    if next_row == bottomRow: next_row = 0


# Start counting through the array elements, L->R.
    for j in range(0, array_length-1):
        U = unit_array[top_row,j]
        NN = unit_array[next_row,j+1]
        NNN = unit_array[top_row, j+1]

        TripletValueList = computeSpecificTripletZVariable (U, NN, NNN)

    
    # Debug print statements
        # Debug print statements
        if not detailed_debug_print_off:
            print()
            print( 'Debug printing: computeConfigZEvenUpperToLower')  #debug_print_off false        
            print( "Returning from computeSpecificTripletZVariable" )
            print( "Unpacking the specific triplet value found")
    
        z1_h_incr = TripletValueList[0]
        left_z2_h_incr = TripletValueList[1]
        right_z2_h_incr = TripletValueList[2]
        z3_h_incr = TripletValueList[3]
        z4_h_incr = TripletValueList[4] 
        left_z5_h_incr = TripletValueList[5]
        right_z5_h_incr = TripletValueList[6] 
        z6_h_incr = TripletValueList[7]
            
        # Debug print statements
        if not detailed_debug_print_off:
            if j == 0:
                print()
                print( 'Debug printing: computeConfigZOddUpperToLower, h_incrementing the z(i) in for loop:' ) #debug_print_off false
                print( ' j = ', j, '   U = ', U, '   NN = ', NN, '   NNN = ', NNN)
                debugPrintZ (z1_h_incr, left_z2_h_incr, right_z2_h_incr, z3_h_incr, z4_h_incr, left_z5_h_incr, right_z5_h_incr, z6_h_incr)
    
                                    
        z1_partial = z1_partial + z1_h_incr 
        left_z2_partial = left_z2_partial + left_z2_h_incr
        right_z2_partial = right_z2_partial + right_z2_h_incr 
        z3_partial = z3_partial + z3_h_incr 
        z4_partial = z4_partial + z4_h_incr 
        left_z5_partial = left_z5_partial + left_z5_h_incr 
        right_z5_partial = right_z5_partial + right_z5_h_incr 
        z6_partial = z6_partial + z6_h_incr 
      
        # Completed for loop; have gone through entire two-row zigzag with triplets that are upper-to-lower-then-upper; 
        #   no wrap-arounds
        
    zPartialArray[0] = z1_partial
    zPartialArray[1] = left_z2_partial
    zPartialArray[2] = right_z2_partial
    zPartialArray[3] = z3_partial        
    zPartialArray[4] = z4_partial  
    zPartialArray[5] = left_z5_partial
    zPartialArray[6] = right_z5_partial    
    zPartialArray[7] = z6_partial    
    return (zPartialArray)   



####################################################################################################
####################################################################################################
#
# Function to compute the set of configuration variables z'i going upper-to-lower across two rows
#    starting with an EVEN row (0 to 1, 2 to 3, etc.), as the SECOND STEP in doing the ODD-to-EVEN
# Function returns a list configvar containing the six z configuration variables:
#    z1 & z2 & z3 & z4 & z5 & z6 
#
####################################################################################################
####################################################################################################



def computeConfigZOddLowerToUpper (array_size_list, unit_array, top_row):

    array_length = array_size_list [0]
    array_layers = array_size_list [1]
    unit_array = unit_array
    
# Create the array to hold the partial (the h_increments in the) z'i's, and populate it with zeros
    zPartialArray = np.zeros((array_length), dtype=np.int)

  
    z1_partial = left_z2_partial = right_z2_partial = z3_partial = 0
    z4_partial = left_z5_partial = right_z5_partial = z6_partial = 0


    next_row = top_row + 1
    bottomRow = array_layers
    if next_row == bottomRow: next_row = 0

    if not z_debug_print_off:
        if top_row == 0: 
            print() 
            print( ' In OddToEven')
            print( ' top_row = ', top_row)
            print( ' bottomRow = ', bottomRow)
            print( ' next_row = ', next_row)
            print()    
        
# NOTE: We are computing the SECOND row of triplets in a zigzag chain,
#   going from the lower row to the top
# Start counting through the array elements, L->R.
    for j in range(0, array_length-1):
        U = unit_array[next_row,j]
        NN = unit_array[top_row,j]
        NNN = unit_array[next_row, j+1]
        if not detailed_debug_print_off:
            if top_row == 15:
                print( ' For top_row = ', top_row, ' and j = ', j, ', then U = ', U, ' and next_row = ', next_row, ' and j+1 is ', j+1, ' and NN = ', NN )
        TripletValueList = computeSpecificTripletZVariable (U, NN, NNN)

    
  


    
    # Debug print statements
        if not detailed_debug_print_off:
            print()
            print( "Computing the SECOND row of a zigzag chain for leading unit ", j)
            print( "Returning from computeSpecificTripletZVariable")
            print( "Unpacking the specific triplet value found")
    
        z1_h_incr = TripletValueList[0]
        left_z2_h_incr = TripletValueList[1]
        right_z2_h_incr = TripletValueList[2]
        z3_h_incr = TripletValueList[3]
        z4_h_incr = TripletValueList[4] 
        left_z5_h_incr = TripletValueList[5]
        right_z5_h_incr = TripletValueList[6] 
        z6_h_incr = TripletValueList[7]

        # Debug print statements
        if not detailed_debug_print_off:
            print()
            print( 'Debug printing: computeConfigZEvenLowerToUpper, h_incrementing the z(i) in for loop:')  #debug_print_off false
            debugPrintZ (z1_h_incr, left_z2_h_incr, right_z2_h_incr, z3_h_incr, z4_h_incr, left_z5_h_incr, right_z5_h_incr, z6_h_incr)
    
                                    
        z1_partial = z1_partial + z1_h_incr 
        left_z2_partial = left_z2_partial + left_z2_h_incr
        right_z2_partial = right_z2_partial + right_z2_h_incr 
        z3_partial = z3_partial + z3_h_incr 
        z4_partial = z4_partial + z4_h_incr 
        left_z5_partial = left_z5_partial + left_z5_h_incr 
        right_z5_partial = right_z5_partial + right_z5_h_incr 
        z6_partial = z6_partial + z6_h_incr 

        # Completed for loop; have gone through entire two-row zigzag with triplets that are lower-to-upper-to-lower; 
        #   no wrap-arounds
        

    zPartialArray[0] = z1_partial
    zPartialArray[1] = left_z2_partial
    zPartialArray[2] = right_z2_partial
    zPartialArray[3] = z3_partial        
    zPartialArray[4] = z4_partial  
    zPartialArray[5] = left_z5_partial
    zPartialArray[6] = right_z5_partial    
    zPartialArray[7] = z6_partial    
    return (zPartialArray)  


############################################
#
# print function: two rows
# 
#-------------------------------------------

def printOddToEvenRows (top_row, unit_array):

    next_row = top_row + 1
    bottom_row = array_layers - 1
    if top_row == bottom_row: next_row = 0
    print( ' *************************')
    print()     
    print( 'top_row = ', top_row, ' next_row = ', next_row)
    print( 'Row', top_row, ':', blnkspc, )
    print (blnkspc),
    for j in range(0,array_length):
        print( unit_array[top_row,j], blnkspc,)
    print ()

    print( 'Row ', next_row, ':', blnkspc,)
    for j in range(0,array_length):
        print( unit_array[next_row,j], blnkspc,)
    print() 
    print( ' *************************')
    return

####################################################################################################
####################################################################################################
#
# Function to compute the set of configuration variables z'i starting with an ODD row (0 to 1, 2 to 3, etc.) 
# Function returns a list configvar containing the six z configuration variables:
#    z1 & z2 & z3 & z4 & z5 & z6 
#
####################################################################################################
####################################################################################################



def compute_config_z_variables_odd_to_even (array_size_list, unit_array, topRow):


####################################################################################################
# This section unpacks the input variable array_size_list
####################################################################################################

    array_length = array_size_list [0]
    array_layers = array_size_list [1]
    unit_array = unit_array

    top_row = topRow
    next_row = topRow + 1
    if top_row == array_layers-1: next_row = 0

# Debug print statements
    if not detailed_debug_print_off:
        print()
        print( "Just entered compute_config_Z_variables: Odd-to-Even")

    if not z_debug_print_off:
        printOddToEvenRows (top_row, unit_array)

###################################################################################################
#
# Compute the triplet values z(i)
#
###################################################################################################

# Initialize the z'i variables

    z1 = z2 = z3 = z4 = z5 = z6 = 0

    z1_total = left_z2_total = right_z2_total = z3_total = 0   
    z4_total = left_z5_total = right_z5_total = z6_total = 0   
    y1_total = left_y2_total = right_y2_total = y3_total = 0          

###################################################################################################
#
# For an ODD-to-EVEN row combination (1 & 2, 3 & 4, etc): 
# Compute the triplet values z(i) for the case of 
# downward-right-then-upwards-right-pointing triplets, from top to next layer down
# going L->R across the zigzag array
# This step does not include any wrap-arounds
#
# This is identical with the corresponding step in Even-to-Odd; the variance is in 
#   (1) The nearest-neighbor on the L-R going down is at j+1, not j, and
#   (2) how the wrap-arounds are computed
#
###################################################################################################

# Start counting through the layers; since we will work with a pair of 
# overlapping layers (for diagonal nearest-neighbors), we use a count of
# layers - 1. 

# Debug print statements
    if not detailed_debug_print_off:
        print()
        print( "Calling computeSpecificTripletZVariable")

    U  = NN = NNN = 0
    
    TripletValueList = list()  
 
# commenting out for debug   
#for i in range(0,array_layers-1):
 #       top_row = i
 #      next_row = i+1

# multiply defined with earlier in function
#    top_row = topRow
#    next_row = topRow + 1  
  
    z1_partial = left_z2_partial = right_z2_partial = z3_partial = 0
    z4_partial = left_z5_partial = right_z5_partial = z6_partial = 0

# Compute the first contribution to z'i's from the top-to-bottom-to-top row
# Even though we're computing the zigzags for an ODD-to-EVEN zigzag chain, 
#   we can start with the EvenUpperToLower computation (same as with EVEN-to-ODD)
#   and the endpoint for the chain is the same. 
# Note that in the next step (after this), the wrap-around triplet will be different. 
    zPartialArray = computeConfigZOddUpperToLower (array_size_list, unit_array, top_row) 
             
# Unpack the new z(i) contributions into the partial values for z(i)                                     
    z1_partial= zPartialArray[0]
    left_z2_partial = zPartialArray[1]
    right_z2_partial = zPartialArray[2]  
    z3_partial = zPartialArray[3]        
    z4_partial = zPartialArray[4]   
    left_z5_partial = zPartialArray[5] 
    right_z5_partial = zPartialArray[6]     
    z6_partial = zPartialArray[7]                            
                                                    
# Update the total z'i values:                                                                                                                                                                 
    z1_total = z1_total + z1_partial
    left_z2_total = left_z2_total + left_z2_partial 
    right_z2_total = right_z2_total + right_z2_partial      
    z3_total = z3_total + z3_partial
    z4_total = z4_total + z4_partial
    left_z5_total = left_z5_total + left_z5_partial
    right_z5_total = right_z5_total + right_z5_partial  
    z6_total = z6_total + z6_partial
    z5_total = left_z5_total + right_z5_total
    z2_total = left_z2_total + right_z2_total
    
# Debug section: Print totals for right-downwards-then-upwards triplets
    if not debug_print_off:
        print()
        print( "Subtotals so far (downward-right-then_upwards-right-pointing triplets)")
        print( "Before any wrap-arounds:")
        print( "(A-A-A) z1_total =", z1_total, "(A-A-B) left_z2_total =", left_z2_total )
        print( "(A-B-A) z3_total =", z3_total, "(B-A-A) right_z2_total =", right_z2_total ) 
        print( "(B-A-B) z4_total =", z4_total, "(B-B-A) left_z5_total =", left_z5_total )
        print( "(B-B-B) z6_total =", z6_total, "(A-B-B) right_z5_total =", right_z5_total)



###################################################################################################
#
# Compute the triplet values z(i) for the first wrap-around on an ODD-to-EVEN zigzag
# Start at last unit on top row, use downward-pointing-diagonal for last unit on next row
# Then wrap-around to pick up first unit on top row.
#
###################################################################################################

  

# Debug print statements
    if not detailed_debug_print_off:
        print()
        print( "Computing first wrap-around triplet")
        print( "Calling computeSpecificTripletZVariable")

    U  = NN = NNN = 0

# Note that when we are working an ODD-to-EVEN zigzag chain, the
#   first wrap-around triplet works on different units than when 
#   we are working an EVEN-to-ODD chain. 

    if not detailed_debug_print_off:
        print()
        print( '  debug in first wrap-around triplet in an ODD-to-EVEN zigzag computation')



    U = unit_array[top_row,array_length-1]
    NN = unit_array[next_row,0]
    NNN = unit_array[top_row, 0]

    TripletValueList = computeSpecificTripletZVariable (U, NN, NNN)

    
    # Debug print statements
    if not detailed_debug_print_off:
        print()
        print( "Returning from computeSpecificTripletZVariable" )
        print( "Unpacking the specific triplet value found")
    
    z1_h_incr = TripletValueList[0]
    left_z2_h_incr = TripletValueList[1]
    right_z2_h_incr = TripletValueList[2]
    z3_h_incr = TripletValueList[3]
    z4_h_incr = TripletValueList[4] 
    left_z5_h_incr = TripletValueList[5]
    right_z5_h_incr = TripletValueList[6] 
    z6_h_incr = TripletValueList[7]

    # Debug print statements
    if not debug_print_off:
        print()         
        print( "(A-A-A) z1_h_incr =", z1_h_incr, "(A-A-B) left_z2_h_incr =", left_z2_h_incr )
        print( "(A-B-A) z3_h_incr =", z3_h_incr, "(B-A-A) right_z2_h_incr =", right_z2_h_incr ) 
        print( "(B-A-B) z4_h_incr =", z4_h_incr, "(B-B-A) left_z5_h_incr =", left_z5_h_incr )
        print( "(B-B-B) z6_h_incr =", z6_h_incr, "(A-B-B) right_z5_h_incr =", right_z5_h_incr ) 
        print()
                                    
                                                
# Update the total z'i values:                                                                                                                                                                 
    z1_total = z1_total + z1_h_incr
    left_z2_total = left_z2_total + left_z2_h_incr 
    right_z2_total = right_z2_total + right_z2_h_incr      
    z3_total = z3_total + z3_h_incr
    z4_total = z4_total + z4_h_incr
    left_z5_total = left_z5_total + left_z5_h_incr
    right_z5_total = right_z5_total + right_z5_h_incr  
    z6_total = z6_total + z6_h_incr
    z5_total = left_z5_total + right_z5_total
    z2_total = left_z2_total + right_z2_total
    
# Debug section: Print totals for right-downwards-then-upwards triplets
    if not debug_print_off:
        print() 
        print( "Subtotals so far (downward-right-then_upwards-right-pointing triplets)")
        print( "After adding in the top-layer wrap-around:")
        print( "(A-A-A) z1_total =", z1_total, "(A-A-B) left_z2_total =", left_z2_total )
        print( "(A-B-A) z3_total =", z3_total, "(B-A-A) right_z2_total =", right_z2_total ) 
        print( "(B-A-B) z4_total =", z4_total, "(B-B-A) left_z5_total =", left_z5_total )
        print( "(B-B-B) z6_total =", z6_total, "(A-B-B) right_z5_total =", right_z5_total )



###########################################################
#
# Compute the triplet values z(i) for the case of an on an ODD-to-EVEN zigzag pass
# with upward-right-pointing diagonals, from next-to-top layer up to  
# the top layer, then going down again, going L->R across the zigzag array
#
###########################################################

# Recall that we are carrying forward previously-computed totals
# for the z(i) values. 
# However, we will re-initialize the partial totals
# so that we get a partial total across each layer of the zigzag

# WThe same count of units will work for the zigzag chain as what we did when we
#   had an EVEN-to-ODD zigzag. 
    
    zPartialArray = computeConfigZOddLowerToUpper (array_size_list, unit_array, top_row)
             
# Unpack the new z(i) contributions into the partial values for z(i)                                     
    z1_partial= zPartialArray[0]
    left_z2_partial = zPartialArray[1]
    right_z2_partial = zPartialArray[2]  
    z3_partial = zPartialArray[3]        
    z4_partial = zPartialArray[4]   
    left_z5_partial = zPartialArray[5] 
    right_z5_partial = zPartialArray[6]     
    z6_partial = zPartialArray[7]                            
   
             
# Update the total z'i values:                                                                                                                                                                 
    z1_total = z1_total + z1_partial
    left_z2_total = left_z2_total + left_z2_partial 
    right_z2_total = right_z2_total + right_z2_partial      
    z3_total = z3_total + z3_partial
    z4_total = z4_total + z4_partial
    left_z5_total = left_z5_total + left_z5_partial
    right_z5_total = right_z5_total + right_z5_partial  
    z6_total = z6_total + z6_partial
    z5_total = left_z5_total + right_z5_total
    z2_total = left_z2_total + right_z2_total
    
# Debug section: Print totals for right-downwards-then-upwards triplets
    if not z_debug_print_off:
        print( "Subtotals so far (adding in upwards-right-then_downwards-right-pointing triplets)")
        print( "Before the last wrap-around:")
        print( "(A-A-A) z1_total =", z1_total, "(A-A-B) left_z2_total =", left_z2_total )
        print( "(A-B-A) z3_total =", z3_total, "(B-A-A) right_z2_total =", right_z2_total ) 
        print( "(B-A-B) z4_total =", z4_total, "(B-B-A) left_z5_total =", left_z5_total )
        print( "(B-B-B) z6_total =", z6_total, "(A-B-B) right_z5_total =", right_z5_total )


                 
                                                
# Only one step remains.
# We need to compute the second wrap-around for the zigzag chain (to get the total number
# of z'i's . 
                                       


###################################################################################################
#
# Compute the triplet values z(i) for the second wrap-around on an ODD-to-EVEN zigzag pass
# Start at last unit on bottom row, use upward-pointing-diagonal for first unit on upper row 
# Then wrap-around to pick up first unit on bottom row.
#
###################################################################################################

  

# Debug print statements
    if not z_debug_print_off:
        print( "Computing second wrap-around triplet" )
        print( "Calling computeSpecificTripletZVariable" )
    U  = NN = NNN = 0

# Note that the triplet wraparound for the second row of an ODD-to-EVEN zigzag chain
#   uses different units as compared with this same triplet on an EVEN-to-ODD zigzag
    U = unit_array[next_row,array_length-1]
    NN = unit_array[top_row,array_length-1]
    NNN = unit_array[next_row, 0]

    TripletValueList = computeSpecificTripletZVariable (U, NN, NNN)

    
    # Debug print statements
    if not detailed_debug_print_off:
        print( "Returning from computeSpecificTripletZVariable" )
        print( "Unpacking the specific triplet value found" )
    
    z1_h_incr = TripletValueList[0]
    left_z2_h_incr = TripletValueList[1]
    right_z2_h_incr = TripletValueList[2]
    z3_h_incr = TripletValueList[3]
    z4_h_incr = TripletValueList[4] 
    left_z5_h_incr = TripletValueList[5]
    right_z5_h_incr = TripletValueList[6] 
    z6_h_incr = TripletValueList[7]

    # Debug print statements
    if not z_debug_print_off:
        print()          
        print( "(A-A-A) z1_h_incr =", z1_h_incr, "(A-A-B) left_z2_h_incr =", left_z2_h_incr) 
        print( "(A-B-A) z3_h_incr =", z3_h_incr, "(B-A-A) right_z2_h_incr =", right_z2_h_incr ) 
        print( "(B-A-B) z4_h_incr =", z4_h_incr, "(B-B-A) left_z5_h_incr =", left_z5_h_incr )
        print( "(B-B-B) z6_h_incr =", z6_h_incr, "(A-B-B) right_z5_h_incr =", right_z5_h_incr ) 
        print()
                                    
                                                
# Update the total z'i values:                                                                                                                                                                 
    z1_total = z1_total + z1_h_incr
    left_z2_total = left_z2_total + left_z2_h_incr 
    right_z2_total = right_z2_total + right_z2_h_incr      
    z3_total = z3_total + z3_h_incr
    z4_total = z4_total + z4_h_incr
    left_z5_total = left_z5_total + left_z5_h_incr
    right_z5_total = right_z5_total + right_z5_h_incr  
    z6_total = z6_total + z6_h_incr
    z5_total = left_z5_total + right_z5_total
    z2_total = left_z2_total + right_z2_total
    
# Debug section: Print totals for right-downwards-then-upwards triplets
    if not z_debug_print_off:
        print( "Totals for all triplets, after adding in the second wrap-around" )
        print( "(A-A-A) z1_total =", z1_total, "(A-A-B) left_z2_total =", left_z2_total )
        print( "(A-B-A) z3_total =", z3_total, "(B-A-A) right_z2_total =", right_z2_total ) 
        print( "(B-A-B) z4_total =", z4_total, "(B-B-A) left_z5_total =", left_z5_total )
        print( "(B-B-B) z6_total =", z6_total, "(A-B-B) right_z5_total =", right_z5_total )



# This concludes computation of the z'i totals for a complete pass through a zigzag chain



# Lots of steps need to happen next

####################################################################################################
#
# Assign the computed configuration variables to elements of the config_vars_list, 
# which will be passed back to the calling procedure
#
###################################################################################################
    
    z1 = z1_total
    z2 = left_z2_total + right_z2_total
    z3 = z3_total 
    z4 = z4_total
    z5 = left_z5_total + right_z5_total
    z6 = z6_total           
    configVarsZList = (z1, z2, z3, z4, z5, z6)

    return (configVarsZList)



####################################################################################################
####################################################################################################

# This function runs both the even-to-odd and odd-to-even zigzags; it is the first step in building
#   another row on top of the basic 1-D zigzag chain

####################################################################################################

def compute_config_Z_variables (array_size_list, unit_array):

# Initialize the z'i variables

    z1 = z2 = z3 = z4 = z5 = z6 = 0

#    z1_total = left_z2_total = right_z2_total = z3_total = 0   
#    z4_total = left_z5_total = right_z5_total = z6_total = 0   
#    y1_total = left_y2_total = right_y2_total = y3_total = 0  

    if not debug_print_off:
        print()
        print( '  Starting to compute Z variables')
        print( '  Total number of pairs of zigzag chains is: ', pairs)
        print()  
    for i in range (0, pairs):
        topRow = 2*i
        if not z_debug_print_off:
            print( '  Row: ', topRow )
        # Obtain the z(i) values from the first even-to-odd zigzag chain (0 to 1, running top-to-bottom)
        config_vars_z_list_even_to_odd = compute_config_Z_variablesEvenToOdd (array_size_list, unit_array, topRow)
        # Assign the returned results to the local sum for each of the z(i) triplets
        z1 = z1+config_vars_z_list_even_to_odd[0]
        z2 = z2+config_vars_z_list_even_to_odd[1]
        z3 = z3+config_vars_z_list_even_to_odd[2]
        z4 = z4+config_vars_z_list_even_to_odd[3]
        z5 = z5+config_vars_z_list_even_to_odd[4]
        z6 = z6+config_vars_z_list_even_to_odd[5]


# Debug section: Print totals for right-downwards-then-upwards triplets
        if not z_debug_print_off:
            print()
            print( 'Starting for loop with i = ', i)
            print( ' -----------')
            print()
            print( "Totals for all triplets, after completing Row: ", topRow, "downwards-to-upwards")
            print( "             (A-A-A) z1_total =", z1)
            print( "(A-A-B) plus (B-A-A) z2_total =", z2)
            print( "             (A-B-A) z3_total =", z3 )  
            print( "             (B-A-B) z4_total =", z4)
            print( "(B-B-A) plus (A-B-B) z5_total =", z5) 
            print( "             (B-B-B) z6_total =", z6)
            print()
            print( ' -----------')
            print()
                
        # Start working on the next zigzag chain        
        topRow = 2*i+1
        # Obtain the z(i) values from the first odd-to-even-to zigzag chain (1 to 2, running top-to-bottom)

        config_vars_z_list_odd_to_even = compute_config_z_variables_odd_to_even (array_size_list, unit_array, topRow)                    
                                        
        # Add the returned results to the local sum for each of the z(i) triplets
        z1 = z1+config_vars_z_list_odd_to_even[0]
        z2 = z2+config_vars_z_list_odd_to_even[1]
        z3 = z3+config_vars_z_list_odd_to_even[2]
        z4 = z4+config_vars_z_list_odd_to_even[3]
        z5 = z5+config_vars_z_list_odd_to_even[4]
        z6 = z6+config_vars_z_list_odd_to_even[5]


# Debug section: Print totals for right-upwards-then-downwards triplets
        if not z_debug_print_off:        
            print( ' -----------')
            print()
            print( "Totals for all triplets, after completing Row: ", topRow, "upwards-to-downwards")
            print( "             (A-A-A) z1_total =", z1)
            print( "(A-A-B) plus (B-A-A) z2_total =", z2)
            print( "             (A-B-A) z3_total =", z3)   
            print( "             (B-A-B) z4_total =", z4)
            print( "(B-B-A) plus (A-B-B) z5_total =", z5) 
            print( "             (B-B-B) z6_total =", z6)
            print()
            print( ' -----------')
            print()
            print( 'Closing a pass through for loop with i = ', i )    
            print()    
# NOTE: Still need to write the computation for an extra odd row in grid, if it exists

    configVarsZList = (z1, z2, z3, z4, z5, z6)                                                                                                                                                                        
    return (configVarsZList)

####################################################################################################
####################################################################################################
#
# Function to compute the set of configuration variables as fractions, given the input of the 
#   whole-number (total count) of configuration variables
# Function returns a list configvar containing the fourteen configuration variables:
#   x1 & x2
#   y1 & y2 & y3
#   w1 & w2 & w3
#   z1 & z2 & z3 & z4 & z5 & z6
#
####################################################################################################
####################################################################################################

def compute_config_vars_fractions (config_vars_list, config_vars_frac_list):
    
    total_units = float(array_length*array_layers)
    total_units_times_two = total_units*2.0
                       
    x1 = float(config_vars_list[0])/total_units
    x2 = float(config_vars_list[1])/total_units 
    
    y1 = float(config_vars_list[2])/total_units_times_two
    y2 = float(config_vars_list[3])/(2.0*total_units_times_two)
    y3 = float(config_vars_list[4])/total_units_times_two 
    
    sumY = y1 + 2.*y2 + y3

    w1 = float(config_vars_list[5])/total_units_times_two
    w2 = float(config_vars_list[6])/(2*total_units_times_two)
    w3 = float(config_vars_list[7])/total_units_times_two        

    sumW = w1 + 2.*w2 + w3

    z1 = float(config_vars_list[8])/total_units_times_two
    z2 = float(config_vars_list[9])/(2.0*total_units_times_two)
    z3 = float(config_vars_list[10])/total_units_times_two 
    z4 = float(config_vars_list[11])/total_units_times_two
    z5 = float(config_vars_list[12])/(2.0*total_units_times_two)
    z6 = float(config_vars_list[13])/total_units_times_two 
              
    sumZ = z1 + 2.*z2 + z3 + z4 + 2.*z5 + z6

    
# Create the master configuration variables fraction List; config_vars_frac_list, 
#   assign configuration variable values, and return the list to the calling procedure 
# Note: Unlike the primary (whole number) config_vars_list, this list DOES NOT contain unit_array    
               
    config_vars_frac_list = (x1, x2, y1, y2, y3, w1, w2, w3, z1, z2, z3, z4, z5, z6)      
    
    return(config_vars_frac_list)
    

####################################################################################################
####################################################################################################
#
# Function to compute the set of configuration variables x'i, y'i, w'i, and z'i
# Note that the x'i were also computed during matrix initialization (cross-check)
# Function returns a list configvar containing the fourteen configuration variables:
#   X1 & X2
#   Y1 & Y2 & Y3
#   W1 & W2 & W3
#   Z1 & Z2 & Z3 & Z4 & Z5 & Z6
#
####################################################################################################
####################################################################################################


def compute_config_variables (array_size_list, unit_array):
    
# Define all the configuration variables (x, y, w, and z) as elements of their respective lists,
#   and assign them their equilibrium values when all enthalpy parameters are set to 0.
#   Then, obtain the actual values for the configuration variables by calling procedures
#   specific to each configuration variable type (x, y, w, and z). 
#   Assign the returned configuration variables to list elements in the main
#   config_vars_list, which is returned to teh calling procedure. 

    if not debug_print_off:
        print ("In compute_config_variables")

#   Initialize the empty list for the full set of configuration variables
    config_vars_list = list() # empty list
    
    configXVarsList = list() # empty list
    configXVarsList = compute_config_X_variables (array_size_list, unit_array)
    X1 = configXVarsList[0]
    X2 = configXVarsList[1]    
    
    configYVarsList = list() # empty list

    configYVarsList = compute_config_Y_variables (array_size_list, unit_array)
    Y1 = configYVarsList[0]
    Y2 = configYVarsList[1]
    Y3 = configYVarsList[2]     

    configWVarsList = list() # empty list    

    configWVarsList = computeConfigWVariables (array_size_list, unit_array)
    W1 = configWVarsList[0]
    W2 = configWVarsList[1]
    W3 = configWVarsList[2]     

    configZVarsList = list() # empty list 
    configZVarsList = compute_config_Z_variables (array_size_list, unit_array)    
    Z1 = configZVarsList[0]
    Z2 = configZVarsList[1]
    Z3 = configZVarsList[2]  
    Z4 = configZVarsList[3]
    Z5 = configZVarsList[4]
    Z6 = configZVarsList[5]     

# Create the master Configuration Variables List; config_vars_list, assign configuration variable
#   values, and return the list to the calling procedure 
               
    config_vars_list = (X1, X2, Y1, Y2, Y3, W1, W2, W3, Z1, Z2, Z3, Z4, Z5, Z6, unit_array)   
         
    if not debug_print_off:
        print ("In compute_config_variables, about to return to **main**")                
    return (config_vars_list)    






####################################################################################################
####################################################################################################

def compute_zero_activation_analytic_config_variables (zero_activation_analytic_config_vars_list, h):    # empty listh):
    
# Define all the configuration variables (x, y, w, and z) as elements of their respective lists,
#   and assign them their equilibrium values when all enthalpy parameters are set to 0.
#   Then, obtain the actual values for the configuration variables by calling procedures
#   specific to each configuration variable type (x, y, w, and z). 
#   Assign the returned configuration variables to list elements in the main
#   config_vars_list, which is returned to teh calling procedure. 


    x1 = 0.5
    x2 = 0.5
    h_squared = h*h
    den = -h_squared + 6.0*h - 1.0
    z3 = (-h+3.0)*(h+1.0)/(8.0*den)
    s = (3.0*h-1.0)/(-h+3.0)
    z1 = s*z3
    z2 = (0.5-z1-z3)/2.0
    z4 = z3
    z5 = z2
    z6 = z1
    y1 = z1 + z2
    y3 = z5 + z6
    y2 = (1.0 - 2.0*z1 + 2.0*z3)/4.0
    w1 = z1 + z3
    w3 = z6 + z4
    w2 = (1.0 - w1 - w3)/2.0
    
# Create the master Configuration Variables List; config_vars_list, assign configuration variable
#   values, and return the list to the calling procedure 
               
    zero_activation_analytic_config_vars_list = (x1, x2, y1, y2, y3, w1, w2, w3, z1, z2, z3, z4, z5, z6, h)  
                
    return (zero_activation_analytic_config_vars_list)    


    


####################################################################################################
####################################################################################################


def compute_config_variables_analytic_delta (h, config_vars_frac_list):
     

    hSquared = h*h
    z3Deltaorig = (hSquared-3.0)*(hSquared+1.0)/(8.0*(hSquared*hSquared-6*hSquared+1.0))
    s = (1.0-3.0*hSquared)/(hSquared-3.0)
    z3Delta = z3Deltaorig - delta
    z4Delta = z3Deltaorig + delta
    z1Delta = s*z3Deltaorig - delta*0.5
    z6Delta = s*z3Deltaorig + delta*0.5
    z2Delta = (0.5-z1Delta-z3Delta)/2.0 - delta*0.8
    z5Delta = (0.5-z6Delta-z4Delta)/2.0 + delta*0.8
    y1Delta = z1Delta + z2Delta
    y3Delta = z5Delta + z6Delta
    y2Delta = (1.0 - z1Delta + z3Delta - z6Delta + z4Delta)/4.0
    w1Delta = z1Delta + z3Delta
    w3Delta = z6Delta + z4Delta
    w2Delta = (1.0 - w1Delta - w3Delta)/2.0
    x1Delta = y1Delta + y2Delta
    x2Delta = y3Delta + y2Delta       

# Create the master Configuration Variables List; config_vars_list, assign configuration variable
#   values, and return the list to the calling procedure 
               
    compute_config_vars_list_delta = (x1Delta, x2Delta, y1Delta, y2Delta, y3Delta, 
    w1Delta, w2Delta, w3Delta, z1Delta, z2Delta, z3Delta, z4Delta, z5Delta, z6Delta)   
             
    return (compute_config_vars_list_delta)    


####################################################################################################
####################################################################################################
#
# Procedure to print out the final set of configuration variables x'i, y'i, w'i, and z'i
#
####################################################################################################
####################################################################################################


def print_config_vars_comparison_delta (h, config_vars_list, compute_config_vars_list_delta):
  
            
    x1 = config_vars_list[0]
    x2 = config_vars_list[1]
    y1 = config_vars_list[2]
    y2 = config_vars_list[3]
    y3 = config_vars_list[4] 
    w1 = config_vars_list[5]
    w2 = config_vars_list[6]
    w3 = config_vars_list[7]        
    z1 = config_vars_list[8]
    z2 = config_vars_list[9]
    z3 = config_vars_list[10] 
    z4 = config_vars_list[11]
    z5 = config_vars_list[12]
    z6 = config_vars_list[13] 

    x1Delta = compute_config_vars_list_delta[0]
    x2Delta = compute_config_vars_list_delta[1]
    y1Delta = compute_config_vars_list_delta[2]
    y2Delta = compute_config_vars_list_delta[3]
    y3Delta = compute_config_vars_list_delta[4] 
    w1Delta = compute_config_vars_list_delta[5]
    w2Delta = compute_config_vars_list_delta[6]
    w3Delta = compute_config_vars_list_delta[7]        
    z1Delta = compute_config_vars_list_delta[8]
    z2Delta = compute_config_vars_list_delta[9]
    z3Delta = compute_config_vars_list_delta[10] 
    z4Delta = compute_config_vars_list_delta[11]
    z5Delta = compute_config_vars_list_delta[12]
    z6Delta = compute_config_vars_list_delta[13] 
        
    sumX = x1+x2 
    sumXDelta =  x1Delta+x2Delta 
    print()
    print( '  For h =', h, 'The configuration variables are: ')
    print()  
    print( '                x1 = %.4f'  % x1,   '   x1Delta  = %.4f' % x1Delta )
    print( '                x2 = %.4f'  % x2,   '   x2Delta  = %.4f' % x2Delta )
    print( '   Sum of the x(i) = %.4f'  % sumX, ' sumXDelta  = %.4f' % sumXDelta )        
                        
    sumZ = z1+2.0*z2+z3 +z4+2.0*z5+z6
    sumZDelta = z1Delta+2.0*z2Delta+z3Delta +z4Delta+2.0*z5Delta+z6Delta    
    print()
    print( "Totals for the z(i) variables:" )
    print( '        (A-A-A) z1 = %.4f'  % z1,   '   z1Delta  = %.4f' % z1Delta ) 
    print( '(A-A-B & B-A-A) z2 = %.4f'  % z2,   '   z2Delta  = %.4f' % z2Delta )  
    print( '        (A-B-A) z3 = %.4f'  % z3,   '   z3Delta  = %.4f' % z3Delta ) 
    print( '        (B-A-B) z4 = %.4f'  % z4,   '   z4Delta  = %.4f' % z4Delta )
    print( '(A-B-B & B-B-A) z5 = %.4f'  % z5,   '   z5Delta  = %.4f' % z5Delta )
    print( '        (B-B-B) z6 = %.4f'  % z6,   '   z6Delta  = %.4f' % z6Delta )        
    print( '   Sum of the z(i) = %.4f'  % sumZ, ' sumZDelta  = %.4f' % sumZDelta  )       
                        
             
                                     
    sumY = y1+2.0*y2+y3
    sumYDelta = y1Delta+2.0*y2Delta+y3Delta        
    print()
    print( 'Totals for the y(i) variables:' )
    print( '          (A-A) y1 = %.4f'  % y1,   '   y1Delta  = %.4f' % y1Delta )  
    print( '    (A-B & B-A) y2 = %.4f'  % y2,   '   y2Delta  = %.4f' % y2Delta )   
    print( '          (B-B) y3 = %.4f'  % y3,   '   y3Delta  = %.4f' % y3Delta ) 
    print( ' Multiplying by degeneracy factors:' )
    print( '   Sum of the y(i) = %.4f'  % sumY, ' sumYDelta  = %.4f' % sumYDelta )

    
    sumW = w1+2.0*w2+w3
    sumWDelta = w1Delta+2.0*w2Delta+w3Delta      
    print()
    print( 'Totals for the w(i) variables:' )
    print( '        (A---A) w1 = %.4f'  % w1,   '   w1Delta  = %.4f' % w1Delta )
    print( '(A---B & B---A) w2 = %.4f'  % w2,   '   w2Delta  = %.4f' % w2Delta )  
    print( '        (B---B) w3 = %.4f'  % w3,   '   w3Delta  = %.4f' % w3Delta )
    print( ' Multiplying by degeneracy factors:' )
    print( '   Sum of the w(i) = %.4f'  % sumW, ' sumWDelta  = %.4f' % sumWDelta )
    print()         

    return

####################################################################################################
#
# Function to increase the x1 value in the array
#
####################################################################################################

def adjust_matrix_x1_up (array_size_list, unit_array, config_vars_list, h):
    
    total_units = float(array_length*array_layers)
    total_units_times_two = total_units*2.0
            
    x1 = float(config_vars_list[0])/total_units
    x2 = float(config_vars_list[1])/total_units 


# Randomly select a unit; if it is 1, change to 0
    unit_row = randrange(0, array_length)       
    unit_col = randrange(0, array_layers)        

    if not show_progress_adjust_matrix_off: 
        print( '   In adjust_matrix_x1_up;', end =" " )  
        if unit_array[unit_row, unit_col] == 0: 
            unit_array[unit_row, unit_col] = 1
            print ( ' successfully changed unit in [', unit_row, ',', unit_col, '] to 1')
        else: print( ' selected unit is already 1' )  
        print()
            
    return unit_array


####################################################################################################
#
# Function to decrease the x1 value in the array
#
####################################################################################################

def adjust_matrix_x1_down (array_size_list, unit_array, config_vars_list, h):
    
    total_units = float(array_length*array_layers)
    total_units_times_two = total_units*2.0
                       
    x1 = float(config_vars_list[0])/total_units
    x2 = float(config_vars_list[1])/total_units 
            
# Randomly select a unit; if it is 1, change to 0
    unit_row = randrange(0, array_length)       
    unit_col = randrange(0, array_layers) 

    if not show_progress_adjust_matrix_off:       
        print( '   In adjust_matrix_x1_down;', end =" " )   
        if unit_array[unit_row, unit_col] == 1: 
            unit_array[unit_row, unit_col] = 0
            print( ' successfully changed unit in [', unit_row, ',', unit_col, '] to 0' )
        else: print( ' selected unit is already 0'  )      
        print()                        
                                                                        
    return unit_array


####################################################################################################
#
# Function to adjust the array and bring x1 and x2 closer togetehr
#
####################################################################################################
    
def adjust_matrix (array_size_list, unit_array, h, jrange, max_x_dif, x1_target):

    config_vars_list = compute_config_variables (array_size_list, unit_array, h)
        
    total_units = float(array_length*array_layers)
    total_units_times_two = total_units*2.0

    x1_total = config_vars_list[0]
    x2_total = config_vars_list[1]  
    x1 = float(config_vars_list[0])/total_units # Obtain the initial value for x1
    x2 = float(config_vars_list[1])/total_units # Obtain the initial value for x2

    print_x_result (x1, x2, x1_total, x2_total, x1_target, max_x_dif)                 

    steps_taken = 0
    break_out = False
    pos_delta = False
    neg_delta = False
    for j in range (0, jrange, 1):
        x1 = float(config_vars_list[0])/total_units
        x1DeltaPos = x1 - x1_target
        if x1DeltaPos > 0: 
            pos_delta = True
            if x1DeltaPos > max_x_dif: # x1 is too large           
                unit_array = adjust_matrix_x1_down (array_size_list, unit_array, config_vars_list, h)
                config_vars_list = compute_config_variables (array_size_list, unit_array, h)
                steps_taken = steps_taken + 1
                if not show_progress_adjust_matrix_off: 
                    print( ' The actual x1 value is too large; improving the unit array for step j = ', j)
            else: 
                break_out = True
        x1DeltaNeg = x1_target - x1
        if x1DeltaNeg > 0:
            if x1DeltaNeg > max_x_dif: # x1 is too small    
                unit_array = adjust_matrix_x1_up (array_size_list, unit_array, config_vars_list, h)
                config_vars_list = compute_config_variables (array_size_list, unit_array, h)
                steps_taken = steps_taken + 1
                if not show_progress_adjust_matrix_off:                                 
                    print( ' The actual x1 value is too small; improving the unit array for j = ', j )
            else:
                break_out = True
        if break_out: break  

    x1 = float(config_vars_list[0])/total_units # Obtain the final rsulting value for x1

    print()
    if steps_taken == 0:
        print ( ' No changes were made to the initial random array.')
    else:    
        print( ' The resultant value for x1 after ', steps_taken, 'steps is %.4f'  % x1)
        if j < jrange-1:
            print( ' This was done within the allowable range of ', jrange, ' steps')
        else:
            print ( ' More steps than the allowed range of ', jrange, ' steps would be needed to bring x1 into desired tolerance.')

    return unit_array               

                
####################################################################################################
#
# Function to compute the entropy
#
####################################################################################################

def LfFunc(x):
    Lfval = x*log(x)-x
    return (Lfval)

    
####################################################################################################
####################################################################################################
#
# Function to compute the entropy, enthalpy, and free energy for a smooth calculation through 
#    a stepped range of h values
#
####################################################################################################
####################################################################################################


def compute_thermodynamic_vars(eps0, h, config_vars_frac_list):
    
    x1 = config_vars_frac_list[0]
    x2 = config_vars_frac_list[1]  
    
    y1 = config_vars_frac_list[2] 
    y2 = config_vars_frac_list[3]  
    y3 = config_vars_frac_list[4]     
    
    w1 = config_vars_frac_list[5] 
    w2 = config_vars_frac_list[6]  
    w3 = config_vars_frac_list[7]  
     
    z1 = config_vars_frac_list[8] 
    z2 = config_vars_frac_list[9]  
    z3 = config_vars_frac_list[10]      
    z4 = config_vars_frac_list[11] 
    z5 = config_vars_frac_list[12]  
    z6 = config_vars_frac_list[13]  
        
                                                                                                                                                    
    Lfx = LfFunc(x1)+LfFunc(x2)
#    print ' lnFunc (x) = %.4f'  % (Lfx)

    Lfy = LfFunc(y1) + 2.0*LfFunc(y2) + LfFunc(y3)
#    print ' lnFunc (y) = %.4f'  % (Lfy)  
      
    Lfw = LfFunc(w1) + 2.*LfFunc(w2) + LfFunc(w3)
#    print ' lnFunc (w) = %.4f'  % (Lfw)

    Lfz = LfFunc(z1) + 2.*LfFunc(z2) + LfFunc(z3) + LfFunc(z4) + 2.*LfFunc(z5) + LfFunc(z6)
#    print ' lnFunc (z) = %.4f'  % (Lfz)
        
    negS = -(2*Lfy+Lfw-Lfx-2*Lfz)
   
    epsilon1 =  log(h)/2.0  # Epsilon1 is an enthalpy parameter, defined in terms of the actual interaction enthalpy       
                            #  See Kikuchi and Brush (1067) for details    
    enthalpy0 = eps0*x1
    enthalpy1 = 2.*epsilon1*(2.*y2 - y1 - y3)    
    
    free_energy = enthalpy0 + enthalpy1 + negS          
        
    sys_vals_list = (negS, enthalpy0, enthalpy1, free_energy)
                                  
    return (sys_vals_list)   
         


####################################################################################################
#
# Function to fill a new unit matrix so it is identical to the starting matrix
#
####################################################################################################
    
def create_identical_unit_array (array_size_list, unit_array):
    
    local_array_length = array_size_list[0]
    local_array_layers = array_size_list[1]
    new_unit_array = np.random.choice([0, 1],size=(local_array_layers,local_array_length))

#  How to copy a list: b = [x for x in a] 
#  Another way: b = list(a)   
         
# We have an  L x M array of units, where M (across) =', localArrayLength, 'and L (layers) =', localArrayLayers
             
    for row in range (0, local_array_layers, 1):
        for column in range (0, local_array_length, 1): 
            new_unit_array[row, column] =  unit_array[row, column]              
                        
    return new_unit_array

####################################################################################################
#
# Function to find an x1 node, which will be a candidate for flipping to lower the total FE value
#
####################################################################################################

def find_candidate_x1_node (array_size_list, unit_array, max_range, titillate_FEMinimum_vals_details_Bool_off):
    
    array_length = array_size_list[0]
    array_layers = array_size_list[1]
    
    candidate_X1_row_col_list = (0,0,0)    
                       
    successBool = 0
    for k in range (0, max_range, 1):
        unit_row = randrange(0, array_length)       
        unit_col = randrange(0, array_layers)

        if unit_array[unit_row, unit_col] == 1:
            candidate_X1_row_col_list = [unit_row, unit_col, k]
            numCandidates = k
            successBool = 1
            if not titillate_FEMinimum_vals_details_Bool_off:
                print()
                print( '  For k = ', k, 'Candidate x1 found at Row: ', unit_row, 'Col: ', unit_col)
                num_candidates = k
                successBool = 1
            break
        # if not findFEMinimumValsDetailsBoolOff:
    
    if successBool == 0:           
        if not titillate_FEMinimum_vals_details_Bool_off:
            print()
            print( ' In findCandidateX2node' ) 
            print( ' Did not find a suitable x2 candidate after ', max_range, ' attempts' )
            print()  
                                    
    return candidate_X1_row_col_list
    # END findCandidateX1 node


####################################################################################################
#
# Function to find an x2 node, which will be a candidate for flipping to lower the total FE value
#
####################################################################################################

def find_candidate_x2_node (array_size_list, unit_array, max_range, titillate_FEMinimum_vals_details_Bool_off):
    
    array_length = array_size_list[0]
    array_layers = array_size_list[1]
    
    candidate_X2_row_col_list = (0,0,0)    

    successBool = 0                                              
    for k in range (0, max_range, 1):
        unit_row = randrange(0, array_length)       
        unit_col = randrange(0, array_layers)

        if unit_array[unit_row, unit_col] == 0:
            candidate_X2_row_col_list = [unit_row, unit_col, k]
            successBool = 1
            num_candidates = k            
            if not titillate_FEMinimum_vals_details_Bool_off:
                print()
                print( '  For k = ', k, 'Candidate x2 found at Row: ', unit_row, 'Col: ', unit_col)
            break           
        # if not findFEMinimumValsDetailsBoolOff:            
        
    if successBool == 0: 
        if not titillate_FEMinimum_vals_details_Bool_off:
            print() 
            print( ' In findCandidateX2node' )
            print( ' Did not find a suitable x2 candidate after ', max_range, ' attempts')
            print()  
                                                                                        
    return candidate_X2_row_col_list
    # END findCandidateX2 node


####################################################################################################
#
# Function to adjust the array (while keeping x1, x2 const) to change the other config variables
#   and bring the FE value to a minimum
#
####################################################################################################
    
def titillate_to_reach_FEMinimum (array_size_list, unit_array, eps0, h, maxRange):

# Note:     sysValsList = (negS, enthalpy0, enthalpy1, freeEnergy): returned from computeThermValues          
# Note:     configVarsList = (x1, x2, y1, y2, y3, w1, w2, w3, z1, z2, z3, z4, z5, z6, unitArray); returned from computeConfigVars     

# Typical value for max_range is about 30
    max_range = 3
    
    titillate_FEminimum_vals_Bool_off = False
    titillate_FEMinimum_vals_details_Bool_off = False 
    
# typical value for total_trials is about 200 or more
    totalTrials = 5
    step = 1 

# initiate all list that will be used

    config_vars_list_orig       = list() 
    config_vars_frac_list_orig  = list()
    x1_candidate_row_col_list   = list()
    x2_candidate_row_col_list   = list()
    config_vars_list_old        = list() 
    config_vars_frac_list_old   = list()
    config_vars_list_new        = list() 
    config_vars_frac_list_new   = list()   
    sys_vals_list_old           = list()
    sys_vals_list_new           = list()   
    
    x1_array       = np.zeros(totalTrials, dtype=np.float)
    y1_array       = np.zeros(totalTrials, dtype=np.float)
    y2_array       = np.zeros(totalTrials, dtype=np.float)
    y3_array       = np.zeros(totalTrials, dtype=np.float)
    w1_array       = np.zeros(totalTrials, dtype=np.float)
    w2_array       = np.zeros(totalTrials, dtype=np.float)
    w3_array       = np.zeros(totalTrials, dtype=np.float)
    z1_array       = np.zeros(totalTrials, dtype=np.float)
    z3_array       = np.zeros(totalTrials, dtype=np.float)       
    negS_array     = np.zeros(totalTrials, dtype=np.float)
    enthalpy0_array     = np.zeros(totalTrials, dtype=np.float)
    enthalpy1_array     = np.zeros(totalTrials, dtype=np.float)
    free_energy_array   = np.zeros(totalTrials, dtype=np.float)  
    
    new_unit_array = create_identical_unit_array (array_size_list, unit_array)            

    # Obtain fractional values for the configuration variables (previously retrieved total counts)        
    totalUnits = float(array_length*array_layers)
    totalUnitsTimesTwo = totalUnits*2.0   

    # Obtain the starting values for the unitArray configuration variables and associated thermodynamic variables        
    config_vars_list_orig = compute_config_variables (array_size_list, unit_array)
    config_vars_frac_list_orig = compute_config_vars_fractions (config_vars_list_orig, config_vars_frac_list_orig)

 
    x1_orig = config_vars_frac_list_orig[0]
    y1_orig = config_vars_frac_list_orig[2]
    y2_orig = config_vars_frac_list_orig[3]
    y3_orig = config_vars_frac_list_orig[4]
    z1_orig = config_vars_frac_list_orig[8]    
    z3_orig = config_vars_frac_list_orig[10]   

    sys_vals_list_orig = compute_thermodynamic_vars(eps0, h, config_vars_frac_list_orig)


    negS_orig       = sys_vals_list_orig[0]
    enthalpy0_orig  = sys_vals_list_orig[1] 
    enthalpy1_orig  = sys_vals_list_orig[2]  
    FE_orig         = sys_vals_list_orig[3]                  
                                
    successful_flips = 0 # Count the total number of successful flips during the totalTrials

    if not titillate_FEminimum_vals_Bool_off:
        print() 
        print( ' In titillate_to_reach_FEMinimum with h =  %.2f' % (h))
        print( ' Starting to adjust config vars to minimize free energy' )                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
        print( '   Trial parameters: x1(starting) =  %.3f' % (x1_orig), ' h =  %.2f' % (h))
        print()
    
            
    for i in range (0, totalTrials, step): 
        x1_candidate_row_col_list = find_candidate_x1_node (array_size_list, unit_array, max_range, titillate_FEMinimum_vals_details_Bool_off)
        x1_row   = x1_candidate_row_col_list[0]
        x1_col   = x1_candidate_row_col_list[1]
        x1_flips = x1_candidate_row_col_list[2]

        x2_candidate_row_col_list = find_candidate_x2_node (array_size_list, unit_array, max_range, titillate_FEMinimum_vals_details_Bool_off)
        x2_row   = x2_candidate_row_col_list[0]
        x2_col   = x2_candidate_row_col_list[1]
        x2_flips = x2_candidate_row_col_list[2]  
             
        if not titillate_FEMinimum_vals_details_Bool_off:
            print() 
            print( '  For i = ', i, ' the candidate x1 node is at Row: ', x1_row, ' Column: ', x1_col, ' after ', x1_flips, 'flips.')
            print( '  For i = ', i, ' the candidate x2 node is at Row: ', x2_row, ' Column: ', x2_col, ' after ', x2_flips, 'flips.')

        # Swap the unit values in newUnitArray: keeping x1 the same, but changing the other config variable values
        new_unit_array[x1_row,x1_col] = 0
        new_unit_array[x2_row,x2_col] = 1
                                                                                                      
        # Obtain the configuration variable values from EACH of the old and new unitArrays
        config_vars_list_old = compute_config_variables (array_size_list, unit_array)
        config_vars_frac_list_old = compute_config_vars_fractions (config_vars_list_old, config_vars_frac_list_old)
        config_vars_list_new = compute_config_variables (array_size_list, new_unit_array)
        config_vars_frac_list_new = compute_config_vars_fractions (config_vars_list_new, config_vars_frac_list_new)    
        # Obtain the thermodynamic variables corresponding to each of the old and new unitArrrays
        sys_vals_list_old = compute_thermodynamic_vars(eps0, h, config_vars_frac_list_old)
        
        negS_old        = sys_vals_list_old[0]
        enthalpy0_old   = sys_vals_list_old[1]          
        enthalpy1_old   = sys_vals_list_old[2]  
        FE_old          = sys_vals_list_old[3]    
   
        sys_vals_list_new = compute_thermodynamic_vars(eps0, h, config_vars_frac_list_new)

        negS_new        = sys_vals_list_new[0]
        enthalpy0_new   = sys_vals_list_new[1]           
        enthalpy1_new   = sys_vals_list_new[2] 
        FE_new          = sys_vals_list_new[3]      
                
        x1_old = config_vars_frac_list_old[0]
        x1_new = config_vars_frac_list_new[0]
        y2_old = config_vars_frac_list_old[3]
        y2_new = config_vars_frac_list_new[3]        
        z1_old = config_vars_frac_list_old[8]
        z1_new = config_vars_frac_list_new[8] 
        z3_old = config_vars_frac_list_old[10]
        z3_new = config_vars_frac_list_new[10] 


        x1_array[i]   = x1_new
        y2_array[i]   = y2_new
        z1_array[i]   = z1_new
        z1_array[i]   = z3_new      
        negS_array[i]      = negS_new
        enthalpy0_array[i]      = enthalpy1_new
        enthalpy1_array[i]      = enthalpy1_new        
        free_energy_array[i]    = FE_new 
        
               
        
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
        # Print the results of the single flip: Note thermodynamic variables
        if not titillate_FEMinimum_vals_details_Bool_off:                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
            print()
            print() 
            print( ' The results of a single flip in two unitArray positions, at trial number', i, ' of ', totalTrials, ' total trials') 
            print()        
            print( '          x1      y2       z1       z3       negS    enthalpy1   FE ')
            print( ' Old:   %.4f' % (x1_old), '  %.4f' % (y2_old), '  %.4f' % (z1_old) , '  %.4f' % (z3_old) , '  %.4f' % (negS_old),   '  %.4f' % (enthalpy1_old),  '  %.4f' % (FE_old)  )     
            print( ' New:   %.4f' % (x1_new), '  %.4f' % (y2_new), '  %.4f' % (z1_new) , '  %.4f' % (z3_new), '  %.4f' % (negS_new),   '  %.4f' % (enthalpy1_new),  '  %.4f' % (FE_new)   ) 
            print() 

        # If the flip was successful, keep the values (assign to unitArray)
        #  Otherwise, revert newUnitArray back to its earlier state
        successBool = 0
        if FE_new < FE_old: successBool = 1
        if successBool:
#            print ' ' 
#             print ' Original unitArray at Row = ', x1Row, ' Col = ', x1Col, ' is ', unitArray[x1Row,x1Col]
#             print ' Original unitArray at Row = ', x2Row, ' Col = ', x2Col, ' is ', unitArray[x2Row,x2Col]
            unit_array[x1_row, x1_col] = 0
            unit_array[x2_row, x2_col] = 1
#             print ' After flip:'
#             print '      New unitArray at Row = ', x1Row, ' Col = ', x1Col, ' is ', unitArray[x1Row,x1Col]
#             print '      New unitArray at Row = ', x2Row, ' Col = ', x2Col, ' is ', unitArray[x2Row,x2Col]
            successful_flips = successful_flips + 1
            x1_end = x1_new
            y2_end = y2_new
            z1_end = z1_new
            z3_end = z3_new
            negS_end = negS_new
            enthalpy1_end = enthalpy1_new
            FE_end = FE_new                           
            if not titillate_FEMinimum_vals_details_Bool_off:                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
                print()
                print( ' Successful flip: free energy reduced, keeping the change' )        
        else: 
            new_unit_arrayy[x1Row,x1Col] = 1
            new_unit_array[x2Row,x2Col] = 0 
            x1End = x1_old
            y2End = y2_old
            z1End = z1_old
            z3End = z3_old
            negS_end = negS_old
            enthalpy1_end = enthalpy1_old
            FE_end = FE_old 
            if not titillate_FEMinimum_vals_details_Bool_off:                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
                print()
                print( ' Unsuccessful flip: free energy increased, NOT keeping the change')        

#    print ' '
#    print ' In adjustMatrixFEMinimum, completed all computations, for h = ', h
#    print ' '
#    print ' The total number of successful flips out of ', totalTrials, 'is: ', successfulFlips 
#    print ' ' 
#    print ' The total change in the configuration and thermodynamic variables:'
#    print '          x1      y2       z1       z3       negS    enthalpy1   FE '
#    print ' Orig:   %.4f' % (x1Orig), '  %.4f' % (y2Orig), '  %.4f' % (z1Orig) , '  %.4f' % (z3Orig) , '  %.4f' % (negEntropyOrig),   '  %.4f' % (enthalpy1Orig),  '  %.4f' % (FEValueOrig)       
#    print ' End:    %.4f' % (x1End), '  %.4f' % (y2End), '  %.4f' % (z1End) , '  %.4f' % (z3End), '  %.4f' % (negEntropyEnd),   '  %.4f' % (enthalpy1End),  '  %.4f' % (FEValueEnd)    
#    print ' '
#    print ' '                                                                                                          

#    if not titillate_FEMinimum_vals_Bool_off:
#        print_unit_array_modification_results (x1ValsArray, y2ValsArray, 
#            z1ValsArray, z3ValsArray, negSValsArray, enthalpy1Array, freeEnergyArray, h, totalTrials, findFEMinimumValsBoolOff)             
                                     
       
    return unit_array
    # END titillate_to_reach_FEMinimum      
      

        
####################################################################################################
####################################################################################################
#
# Code Documentation - top-down: 
# The MAIN module comprising of calls to:
#  (1) Welcome
#  (2) Obtain array size specifications for an MxN array of 0 or 1 units (currently pre-defined patterns)
#  (3) compute_config_variables: Compute the configuration variables for the entire grid
#    -- NOTE: This is the major work-horse function in this program; see documentation below
#   
#  Documentation for compute_config_variables and its supporting functions:  
#  (3) compute_config_variables: Computes the configuration variables for the grid; x, y, w, and z.
#    (3.1) compute_config_X_variables
#    (3.2) compute_config_Y_variables
#    (3.3) computeConfigWVariables
#    (3.4) compute_config_Z_variables
#
#    NOTE: The configuration variables that are computed here are the TOTALS; 
#        that is, X(i) (2 vars), Y(i) (3 vars), W(i) (3 vars), and Z(i) (6 vars)
#        They can be converted to the fractional values by dividing by the total number of units;
#        this is not needed until computing energy terms.  
#    NOTE: This function calls two primary functions, alternating as it works through the grid rows: 
#      - Compute the X, Y, W, and Z (total) contributions from the initial (or next) odd-to-even rows
#      - Compute the X, Y, W, and Z (total) contributions from the initial (or next) even-to-odd rows
#    NOTE: Because of the way in which the zigzag rows are aligned, computing the wrap-around 
#      configuration variables differs depending on whether working with odd-to-even or even-to-odd.
#
#    (3.1) compute_config_X_variables: Computes the total values for X1 & X2 in the grid, going row-by-row.
#        NOTE: Unlike all the other configuration variable computations, requiring a wrap-around to the
#        units at the beginning of the row(s), this computation is a simple count-up of "A" and "B" units. 
#
#    (3.2) compute_config_Y_variables: Computes the total values for Y1, Y2, and Y3 in the grid.
#        (3.2.1) computeConfigYEvenRowZigzagVariables
#        (3.2.1) computeConfigYOddRowZigzagVariables
#        NOTE: The "Y" nearest-neighbor variables are computed going left-to-right across a row. 
#        There are separate calculations for the even and odd rows due to how the wrap-arounds 
#           are computed. 
#
#    (3.3) computeConfigWVariables: Computes the total values for W1, W2, and W3 in the grid.
#        (3.3.1) computeConfigWHorizontalRowVariables  
#        (3.3.2) computeConfigWVerticalColVariables
#        NOTE: For each unit in the grid, there are four "w" next-nearest-neighbors.
#            If we computed each next-nearest-neighbor for each unit, we'd have to divide by two, 
#            to account for degeneracy. However, in this code, we only compute each NNN once, and 
#            do not need to divide by two.  
#            Because the "w" variables are computed BOTH for the next-nearest-neighbors DIRECTLY ABOVE
#            and DIRECTLY BELOW a given unit, as well as to the LEFT and RIGHT of each unit, there are
#            two kinds of computations ... one working with the vertical, and another with the horizontal.
#            The VERTICAL computations are done once for each pass through a pair of rows; they do not need to 
#            be divided by two (for degeneracy) as each vertical next-nearest-neighbor is only computed once. 
#            Similarly, the horizontal next-nearest-neighbors are also computed going left-to-right across 
#            a row, and are only computed once. 
#
#    (3.4) compute_config_Z_variables: Computes the total values for Z1, Z2, Z3, Z4, Z5, and Z6 in the grid.
#
#
####################################################################################################
####################################################################################################

def main():

####################################################################################################
# Obtain unit array size in terms of array_length (M) and layers (N)
####################################################################################################                

    global array_length
    global array_layers
    global even_layers
    global pairs
    
    global blnkspc  
    
    global debug_print_off
    global detailed_debug_print_off
    global z_debug_print_off
    global show_progress_adjust_matrix_off
    global explanation_thermodynamic_plot_off

    even_layers = True

    blnkspc=' '

# Define values for the global debug variables
    debug_print_off = True
    detailed_debug_print_off = True
    z_debug_print_off = True
    show_progress_adjust_matrix_off = True
    explanation_thermodynamic_plot_off = True
    perform_analytics = True

# This is a local variable; it will be passed to compute_config_variables
#  It will determine whether we print the contents of the x-array at the
#  beginning and end of the adjust-matrix step. 
    beforeAndAfterAdjustedMatrixPrintOff = True
    
    welcome()
    print_debug_status (debug_print_off)

####################################################################################################
# Parmaters needed to define the size of the 2-D CVM grid
#################################################################################################### 
            
    array_size_list   = list() # empty list
    array_size_list   = obtain_array_size_specs () # function call to get the actual dimensions of the 2-D grid
    array_length     = array_size_list[0]
    array_layers     = array_size_list [1]


    if array_layers % 2 == 0: even_layers == True #then an even number of layers

    # Determine the total number of PAIRS of zigzag chains
    pairs_layers = array_layers/2
    pairs = int(pairs_layers + 0.01) 

    print_grid_size_specs ()

####################################################################################################
# Lists needed to hold the configuration variable valuess and the thermodynamic (system) values
#################################################################################################### 

    config_vars_list        = list()    # empty list
    config_vars_frac_list   = list()    # empty list
    sys_vals_list           = list()    # empty list
    zero_activation_analytic_config_vars_list = list()    # empty list
    
####################################################################################################
# Parmaters needed to compute configuration variables
####################################################################################################      
                       
    eps0a = 0.0  # The enthalpy for single unit activation is set to zero for this code
    eps0b = 0.1           
    eps0c = 0.2

####################################################################################################
# User-definable variables influencing the numbers of steps, ranges, etc. to correct the
#  randomly-generated distribution so that it is close to the desired value for x1.  
####################################################################################################                 

    x1_target = 0.45        # Initially defining the target value for x1 as 0.5
    max_x_dif = 0.015       # Maximum difference between x1 from randomly-generated array and the
                            #   equiprobable distribution value of x1 = x1_target
                            #   Typical values are in the range 0.01 - 0.015.
    jrange = 10             # Maximal number of steps allowed to improve the x1 distribution
                            #   Typical values are in the range 100 - 200.             

    max_range = 5           # Used to control range in titillate; typical values > 200 

####################################################################################################
# Compute configuration variables
####################################################################################################                

     

# Determine if we are probabilistically generating a pattern (patternProb = 0)
#   or if we are selecting a pre-stored pattern (patternProb = 1 ... N)

    pattern_select = obtain_pattern_selection()
    print_pattern_selection (pattern_select)
    
 
# Use the initial value of eps0a .. c for eps0
    eps0 = eps0a    
    
# Pick a starting value for h; it should be less than 1; it should be in the realm of 0.7 - 0.8    
    h0 = 1.0
# Pick the number of steps (increasing 0.01) for increasing h
    h_range = 10
# Pick the increment for increasing h
    h_incr = 0.1 
       
    print_run_parameters (h0, h_incr, h_range)             
 
    h_array = np.zeros(h_range, dtype=np.float)

    neg_S_array1 = np.zeros(h_range, dtype=np.float)
    f_eps0_array1 = np.zeros(h_range, dtype=np.float)         
    f_eps1_array1 = np.zeros(h_range, dtype=np.float)  
    f_energy_array1 = np.zeros(h_range, dtype=np.float) 
 
    neg_S_array2 = np.zeros(h_range, dtype=np.float)
    f_eps0_array2 = np.zeros(h_range, dtype=np.float)         
    f_eps1_array2 = np.zeros(h_range, dtype=np.float)  
    f_energy_array2 = np.zeros(h_range, dtype=np.float) 

    neg_S_array3 = np.zeros(h_range, dtype=np.float)
    f_eps0_array3 = np.zeros(h_range, dtype=np.float)         
    f_eps1_array3 = np.zeros(h_range, dtype=np.float)  
    f_energy_array3 = np.zeros(h_range, dtype=np.float) 
        

    h = h0
    if pattern_select == 0:    
        config_vars_list = list ()  # redefine this as an empty list
        unit_array    = initialize_matrix (array_size_list, pattern_select, h)            
    # adjust the unit array so that x1 approximately = x2
        unit_array = adjust_matrix (array_size_list, unit_array, h, jrange, max_x_dif, x1_target)

    # obtain the configuration variables (this step was also done while adjusting the array)       
        config_vars_list = compute_config_variables (array_size_list, unit_array, h)
        config_vars_frac_list = compute_config_vars_fractions (config_vars_list, config_vars_frac_list)
        print_config_vars_fraction_vals (config_vars_frac_list, pattern_select)
        zero_activation_analytic_config_vars_list = compute_zero_activation_analytic_config_variables (zero_activation_analytic_config_vars_list, h)
        print_config_vars_comparison (config_vars_frac_list, zero_activation_analytic_config_vars_list, h)

        h = h - h_incr        
        for i in range (0, h_range, 1):        
            h = h + h_incr

    # Use the initial value of eps0a .. c for eps0
            eps0 = eps0a                                                                     
    # obtain the thermodynamic variables 
            sys_vals_list = compute_thermodynamic_vars (eps0, h, config_vars_frac_list) 
            
# Detailed print statement, not needed after trace            
#            print_thermodynamic_values (h, eps0, sys_vals_list)  
            
    # store the thermodynamic variables to plot later
            h_array[i]=h
            neg_S_array1[i]    = sys_vals_list[0] 
            f_eps0_array1[i]   = sys_vals_list[1]
            f_eps1_array1[i]   = sys_vals_list[2]    
            f_energy_array1[i] = sys_vals_list[3]                                                  

    # Use the second value of eps0a .. c for eps0
            eps0 = eps0b                                                                     
    # obtain the thermodynamic variables 
            sys_vals_list = compute_thermodynamic_vars (eps0, h, config_vars_frac_list)   

# Detailed print statement, not needed after trace
#            print_thermodynamic_values (h, eps0, sys_vals_list)  
            
    # store the thermodynamic variables to plot later
            h_array[i]=h
            neg_S_array2[i]    = sys_vals_list[0] 
            f_eps0_array2[i]   = sys_vals_list[1]
            f_eps1_array2[i]   = sys_vals_list[2]    
            f_energy_array2[i] = sys_vals_list[3]   

    # Use the thirdl value of eps0a .. c for eps0
            eps0 = eps0c                                                                    
    # obtain the thermodynamic variables 
            sys_vals_list = compute_thermodynamic_vars (eps0, h, config_vars_frac_list)  
            
# Detailed print statement, not needed after trace
#            print_thermodynamic_values (h, eps0, sys_vals_list)   
            
    # store the thermodynamic variables to plot later
            h_array[i]=h
            neg_S_array3[i]    = sys_vals_list[0] 
            f_eps0_array3[i]   = sys_vals_list[1]
            f_eps1_array3[i]   = sys_vals_list[2]    
            f_energy_array3[i] = sys_vals_list[3]   



    # end loop of computing configuration variables and thermodyanamic quantities for various values of h


    if pattern_select > 0:   
        h = h0 - h_incr        
        config_vars_list = list ()  # redefine this as an empty list
        unit_array    = initialize_matrix (array_size_list, pattern_select, h)  
        config_vars_list = compute_config_variables (array_size_list, unit_array)
        config_vars_frac_list = compute_config_vars_fractions (config_vars_list, config_vars_frac_list)
        print_config_vars_fraction_vals (config_vars_frac_list, pattern_select)
#        zero_activation_analytic_config_vars_list = compute_zero_activation_analytic_config_variables (h)
#        print_config_vars_comparison (config_vars_frac_list, zero_activation_analytic_config_vars_list) 
        for i in range (0, h_range, 1):        
            h = h + h_incr        
    # obtain the first set of thermodynamic variables 
            eps0 = eps0a
#            print('  For the first epsilon instance, with eps0 = ', eps0a, ':')
            sys_vals_list = compute_thermodynamic_vars (eps0, h, config_vars_frac_list)   
#            zero_activation_analytic_config_vars_list = compute_zero_activation_analytic_config_variables (h)
#            print_config_vars_comparison (config_vars_frac_list, zero_activation_analytic_config_vars_list) 
#            print_thermodynamic_values (h, eps0, sys_vals_list)
    # store the thermodynamic variables to plot later
            h_array[i]=h
            neg_S_array1[i]    = sys_vals_list[0] 
            f_eps0_array1[i]   = sys_vals_list[1]
            f_eps1_array1[i]   = sys_vals_list[2]    
            f_energy_array1[i] = sys_vals_list[3]                                                

    # obtain the second set of thermodynamic variables 
            eps0 = eps0b
#            print('  For the second epsilon instance, with eps0 = ', eps0b, ':')            
            sys_vals_list = compute_thermodynamic_vars (eps0, h, config_vars_frac_list)    

# Detailed print statement, not needed after trace
#            print_thermodynamic_values (h, eps0, sys_vals_list)
    # store the thermodynamic variables to plot later            
            neg_S_array2[i]    = sys_vals_list[0] 
            f_eps0_array2[i]   = sys_vals_list[1]
            f_eps1_array2[i]   = sys_vals_list[2]    
            f_energy_array2[i] = sys_vals_list[3] 

    # obtain the second set of thermodynamic variables 
            eps0 = eps0c
#            print('  For the third epsilon instance, with eps0 = ', eps0c, ':')            
            sys_vals_list = compute_thermodynamic_vars (eps0, h, config_vars_frac_list)   
#            print_thermodynamic_values (h, eps0, sys_vals_list)
    # store the thermodynamic variables to plot later            
            neg_S_array3[i]    = sys_vals_list[0] 
            f_eps0_array3[i]   = sys_vals_list[1]
            f_eps1_array3[i]   = sys_vals_list[2]    
            f_energy_array3[i] = sys_vals_list[3] 

#    print_one_set_of_thermodynamics (pattern_select, eps0, h_array, h_range, 
#                                       neg_S_array1, f_eps0_array1, f_eps1_array1, f_energy_array1)   
    


    print_set_of_three_thermodynamics (pattern_select, eps0a, eps0b, eps0c, h_array, h_range,
                                       neg_S_array1, f_eps0_array1, f_eps1_array1, f_energy_array1, 
                                       neg_S_array2, f_eps0_array2, f_eps1_array2, f_energy_array2,
                                       neg_S_array3, f_eps0_array3, f_eps1_array3, f_energy_array3)

    # end loop of computing configuration variables and thermodyanamic quantities for various values of h
            

    print_thermodynamic_plot_statement(pattern_select)
    
    print()
    print(' For the case where eps0 =  %.4f' % eps0a) 
    print()
    plot_thermodynamic_vals (h_array, neg_S_array1, f_eps0_array1, f_eps1_array1, f_energy_array1)
                                                                                                                              
    print()
    print(' For the case where eps0 =  %.4f' % eps0b) 
    print()
    plot_thermodynamic_vals (h_array, neg_S_array2, f_eps0_array2, f_eps1_array2, f_energy_array2)  
           
    print()
    print(' For the case where eps0 =  %.4f' % eps0c) 
    print()
    plot_thermodynamic_vals (h_array, neg_S_array3, f_eps0_array3, f_eps1_array3, f_energy_array3)


# This ends the primary portion of __main__
# The following code addresses a different task, for the case of using a pre-defined pattern    
    total_units = float(array_length*array_layers)
    x1 = float(config_vars_list[0])/total_units
    y1 = float(config_vars_list[2])/total_units
    y2 = float(config_vars_list[3])/total_units
    y3 = float(config_vars_list[4])/total_units    

    
    eps0_start = eps0a
    eps0_range = 7
    eps0_incr = 0.01
    if pattern_select > 0: 
        eps0 = eps0_start-eps0_incr
        interact_enth_y_delta = 2.0*y2-y1-y3
        print()
        print(' The interactive enthalpy y delta is %.4f' % (interact_enth_y_delta) )
        print('   and is defined as 2*y2 - y1 - y3.')
        print(' The starting value for epsilon0 is  %.4f' % (eps0a) )
        print()
        print(' The enth-delta in the following tables is defined as')
        print('   the activation enthalpy plus the interaction enthalpy.')
        for j in range (0, eps0_range, 1):
            eps0 = eps0 + eps0_incr
            act_enth = eps0*x1            
#            print()
#            print( ' eps0: %.4f' % (eps0), '  act_enth: %.4f' % (act_enth)  )
            h = h0 - h_incr
#            print ('     h       eps1   int_enth   enth-delta')
            for i in range (0, h_range, 1):
                h = h + h_incr
                eps1 = log(h)/2.0 
                int_enth = 2.0*eps1*interact_enth_y_delta
                enthalpy_delta = act_enth + int_enth
#                print( '   %.4f' % (h) , '  %.4f' % (eps1), '   %.4f' % (int_enth) ,'   %.4f' % (enthalpy_delta)  ) 

    h = h0
    eps0 = eps0a
    new_array = titillate_to_reach_FEMinimum (array_size_list, unit_array, eps0, h, max_range)
         
####################################################################################################
# Conclude specification of the MAIN procedure
####################################################################################################                
    
if __name__ == "__main__": main()

####################################################################################################
# End program
####################################################################################################  