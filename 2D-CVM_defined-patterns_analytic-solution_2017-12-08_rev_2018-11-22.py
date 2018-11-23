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
    print( '  and updated 11/22/2018, by A.J. Maren')
    print( 'This version computes the behavior of a perturbed unit_array,')
    print(  ' based on minimizing the free energy both before and after perturbation.')
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

    pattern_select = 0
    return(pattern_select)

####################################################################################################
####################################################################################################
#
# Function to obtain a row of 2-D CVM data - PRESTORED pattern (currently part of a 16x16 grid)
#
####################################################################################################
####################################################################################################


def obtainGridRow (rowNum, patternSelect, h):


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
    if patternSelect == 2:
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
#    rowArray8 =   [0,0,0,0,0,1,1,0, 0,0,1,0,1,1,0,1] # Row 0 - ninth row 
#    rowArray9 =   [0,1,1,1,0,1,0,1, 1,0,0,0,1,0,1,1] # Row 1 - tenth row (counting down from the top)
#    rowArray10 =  [0,0,1,1,1,0,0,1, 1,1,0,0,1,0,1,1] # Row 2 - eleventh row (counting down from the top) 
#    rowArray11 =  [0,1,1,1,1,0,0,0, 1,0,1,0,0,0,1,1] # Row 3 - twelfth row (counting down from the top)
#    rowArray12 =  [0,1,1,1,0,0,0,1, 0,0,1,1,0,1,0,1] # Row 4 - thirteenth row (counting down from the top)
#    rowArray13 =  [1,1,1,0,0,1,0,1, 1,0,1,0,1,1,0,0] # Row 5 - fourteenth row (counting down from the top)
#    rowArray14 =  [0,1,0,0,1,1,1,0, 1,1,0,0,0,1,0,1] # Row 6 - fifteenth row (counting down from the top) 
#    rowArray15 =  [0,0,1,1,1,1,0,0, 1,0,0,1,0,0,0,1] # Row 7 - sixteenth row (counting down from the top)                                                                   

# Rich Club
    rowArray8 =   [1,1,1,0,0,0,0,0, 0,0,0,0,0,1,1,1] # Row 0 - ninth row 
    rowArray9 =   [1,1,1,1,0,0,0,0, 0,0,0,0,0,1,1,1] # Row 1 - tenth row (counting down from the top)
    rowArray10 =  [1,1,1,1,1,0,0,0, 0,0,0,0,0,0,1,1] # Row 2 - eleventh row (counting down from the top) 
    rowArray11 =  [1,1,1,1,1,0,0,0, 0,0,0,0,0,0,1,1] # Row 3 - twelfth row (counting down from the top)
    rowArray12 =  [1,1,1,1,1,1,1,0, 0,0,0,0,0,0,1,1] # Row 4 - thirteenth row (counting down from the top)
    rowArray13 =  [1,1,1,1,1,1,1,0, 0,0,0,0,0,1,1,1] # Row 5 - fourteenth row (counting down from the top)
    rowArray14 =  [1,1,1,1,1,1,1,0, 0,0,0,0,0,0,1,1] # Row 6 - fifteenth row (counting down from the top) 
    rowArray15 =  [1,1,1,1,1,0,0,0, 0,0,0,0,1,1,1,1] # Row 7 - sixteenth row (counting down from the top) 


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

    print()
    print( ' The parameters for this run are: ')
    print( '   Initial interaction enthalpy parameter h is:', h0)
    print( ' The run will begin with an interaction enthalpy of h at:', h0 + h_incr)
    print( '   up through a final value of: ', h0 + h_incr*h_range)
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
# Procedure to print out the full set of configuration variables
#
#    Inputs:    config_vars_list
#
####################################################################################################

def print_config_vars (config_vars_list):


    print () 
    return


####################################################################################################
####################################################################################################
#
# Procedure to randomly-generate an array, and then permute it to achieve the desired 
#    z1 & z3 values.
#
#    Inputs:    array_size_list: a list of two integers; array_length and layers
#               h: the interaction enthalpy parameter
#    Return: the matrix unit_array, a matrix of 0's and 1's.
#
####################################################################################################
####################################################################################################

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
    print( 'With h = ', h, ', we begin with a randomly-generated array:' )
    print()
    
# Bring the array closer to the desired configuration variable values    
    
    return unit_array

####################################################################################################
####################################################################################################
#
# Procedure to initialize the matrix with EITHER a pre-stored pattern of values
#    OR randomly-generate an array, and then permute it to achieve the desired 
#    z1 & z3 values.
#
#    Inputs:    array_size_list: a list of two integers; array_length and layers
#               patternProb: an integer indicating whether to randomly-generate
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

# Determining "pairs" - the total number of pairs of zigzag chains - is done in __main__; "pairs" is a global variable
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


def computeConfigXVariables (array_size_list, unit_array):


####################################################################################################
# This section unpacks the input variable array_size_list
####################################################################################################

    array_length = array_size_list [0]
    array_layers = array_size_list [1]
    unit_array = unit_array

# Debug print statements
    if not debug_print_off:
        print()  
        print( "Just entered computeConfigXVariables")
                                   
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
        
#    print  "Leaving computeConfigXVariables for calling procedure"
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
        # h_increment the y_1; the count of A-A nearest-neighbor pairs:
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
                # h_increment the y_1; the count of A-A nearest-neighbor pairs:
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
                # h_increment the y_1; the count of A-A nearest-neighbor pairs:
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
        # h_increment the y_1; the count of A-A nearest-neighbor pairs:
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

def computeConfigYVariables (array_size_list, unit_array):

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
            print( "           (A-A) z1_total =", y1)
            print( "(A-B) plus (B-A) z2_total =", y2)
            print( "           (B-B) z3_total =", y3 )  
            print()
            print( ' -----------')
            print()
                
        # Start working on the next zigzag chain        
#        topRow = 2*i+1
        # Obtain the z(i) values from the first odd-to-even-to zigzag chain (1 to 2, running top-to-bottom)

#        config_vars_z_list_odd_to_even = compute_config_z_variables_odd_to_even (array_size_list, unit_array, topRow)                    
                                        
        # Add the returned results to the local sum for each of the z(i) triplets
#        z1 = z1+config_vars_z_list_odd_to_even[0]
#        z2 = z2+config_vars_z_list_odd_to_even[1]
#        z3 = z3+config_vars_z_list_odd_to_even[2]
#        z4 = z4+config_vars_z_list_odd_to_even[3]
#        z5 = z5+config_vars_z_list_odd_to_even[4]
#        z6 = z6+config_vars_z_list_odd_to_even[5]


# Debug section: Print totals for right-upwards-then-downwards triplets
#        print ' '
#        print ' -----------'
 #       print ' '
#        print "Totals for all triplets, after completing Row: ", topRow
#        print "             (A-A-A) z1_total =", z1
#        print "(A-A-B) plus (B-A-A) z2_total =", z2
#        print "             (A-B-A) z3_total =", z3   
#        print "             (B-A-B) z4_total =", z4
#        print "(B-B-A) plus (A-B-B) z5_total =", z5 
#        print "             (B-B-B) z6_total =", z6
#        print ' '
#        print ' -----------'
#        print ' '    
# NOTE: Still need to write the computation for an extra odd row in grid, if it exists

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



def computeConfigZVariablesEvenToOdd (array_size_list, unit_array, topRow):


####################################################################################################
# This section unpacks the input variable array_size_list
####################################################################################################

    array_length = array_size_list [0]
    array_layers = array_size_list [1]
    unit_array = unit_array

# Debug print statements
    if not detailed_debug_print_off:
        print()
        print( "Just entered computeConfigZVariables: Even-to-Odd")
                                  
                                                               

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
        print( "Just entered computeConfigZVariables: Odd-to-Even")

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

def computeConfigZVariables (array_size_list, unit_array):

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
        config_vars_z_list_even_to_odd = computeConfigZVariablesEvenToOdd (array_size_list, unit_array, topRow)
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
# Procedure to compute the set of configuration variables x'i, y'i, w'i, and z'i
# Note that the x'i were also computed during matrix initialization (cross-check)
# Procedure returns a list configvar containing the eightfourteen configuration variables:
#   x1 & x2
#   y1 & y2 & y3
#   w1 & w2 & w3
#   z1 & z2 & z3 & z4 & z5 & z6
#
####################################################################################################
####################################################################################################


def compute_config_variables (array_size_list, unit_array, h):
    
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
    configXVarsList = computeConfigXVariables (array_size_list, unit_array)
    x1 = configXVarsList[0]
    x2 = configXVarsList[1]    
    
    configYVarsList = list() # empty list

    configYVarsList = computeConfigYVariables (array_size_list, unit_array)
    y1 = configYVarsList[0]
    y2 = configYVarsList[1]
    y3 = configYVarsList[2]     

    configWVarsList = list() # empty list    

    configWVarsList = computeConfigWVariables (array_size_list, unit_array)
    w1 = configWVarsList[0]
    w2 = configWVarsList[1]
    w3 = configWVarsList[2]     

    configZVarsList = list() # empty list 
    configZVarsList = computeConfigZVariables (array_size_list, unit_array)    
    z1 = configZVarsList[0]
    z2 = configZVarsList[1]
    z3 = configZVarsList[2]  
    z4 = configZVarsList[3]
    z5 = configZVarsList[4]
    z6 = configZVarsList[5]     

# Create the master Configuration Variables List; config_vars_list, assign configuration variable
#   values, and return the list to the calling procedure 
               
    config_vars_list = (x1, x2, y1, y2, y3, w1, w2, w3, z1, z2, z3, z4, z5, z6, unit_array)   
         
    if not debug_print_off:
        print ("In compute_config_variables, about to return to **main**")                
    return (config_vars_list)    

####################################################################################################
####################################################################################################


def compute_config_variables_delta (h):
    
# Define all the configuration variables (x, y, w, and z) as elements of their respective lists,
#   and assign them their equilibrium values when all enthalpy parameters are set to 0.
#   Then, obtain the actual values for the configuration variables by calling procedures
#   specific to each configuration variable type (x, y, w, and z). 
#   Assign the returned configuration variables to list elements in the main
#   config_vars_list, which is returned to teh calling procedure. 


    x1 = 0.5
    x2 = 0.5
    hSquared = h*h
    z3 = (hSquared-3.0)*(hSquared+1.0)/(8.0*(hSquared*hSquared-6*hSquared+1.0))
    s = (1.0-3.0*hSquared)/(hSquared-3.0)
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
               
    config_vars_list = (x1, x2, y1, y2, y3, w1, w2, w3, z1, z2, z3, z4, z5, z6)  
                
    return (config_vars_list)    


    


####################################################################################################
####################################################################################################


def compute_config_variables_analytic (h, config_vars_list):
    
# Define all the configuration variables (x, y, w, and z) as elements of their respective lists,
#   and assign them their equilibrium values when all enthalpy parameters are set to 0.
#   Then, obtain the actual values for the configuration variables by calling procedures
#   specific to each configuration variable type (x, y, w, and z). 
#   Assign the returned configuration variables to list elements in the main
#   config_vars_list, which is returned to teh calling procedure. 

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
    w1Delta, w2Delta, w3Delta, 
    z1Delta, z2Delta, z3Delta, z4Delta, z5Delta, z6Delta)   
             
    return (compute_config_vars_list_delta)    


####################################################################################################
####################################################################################################
#
# Procedure to print out the final set of configuration variables x'i, y'i, w'i, and z'i
#
####################################################################################################
####################################################################################################


def print_config_vars_comparison(h, config_vars_list, compute_config_vars_list_delta):
  
            
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
# Function to h_increase the x1 value in the array
#
####################################################################################################

def adjust_matrix_x1_up (array_size_list, unit_array, config_vars_list, h):
    
    total_units = float(array_length*array_layers)
    total_units_times_two = total_units*2.0
            
    x1 = float(config_vars_list[0])/total_units
    x2 = float(config_vars_list[1])/total_units 

    print( '   In adjust_matrix_x1_up;', end =" " )  
# Randomly select a unit; if it is 1, change to 0
    unit_row = randrange(0, array_length)       
    unit_col = randrange(0, array_layers)        

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
    
    print( '   In adjust_matrix_x1_down;', end =" " )   
        
# Randomly select a unit; if it is 1, change to 0
    unit_row = randrange(0, array_length)       
    unit_col = randrange(0, array_layers)        

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
                print( ' The actual x1 value is too large; improving the unit array for step j = ', j)
                unit_array = adjust_matrix_x1_down (array_size_list, unit_array, config_vars_list, h)
                config_vars_list = compute_config_variables (array_size_list, unit_array, h)
                steps_taken = steps_taken + 1
            else: 
                break_out = True
        x1DeltaNeg = x1_target - x1
        if x1DeltaNeg > 0:
            if x1DeltaNeg > max_x_dif: # x1 is too small    
                print( ' The actual x1 value is too small; improving the unit array for j = ', j )
                unit_array = adjust_matrix_x1_up (array_size_list, unit_array, config_vars_list, h)
                config_vars_list = compute_config_variables (array_size_list, unit_array, h)
                steps_taken = steps_taken + 1
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
# Function to compute the entropy
#
####################################################################################################
####################################################################################################


def computeEntropyV2(config_vars_list):
    print()
    print( '  In computeEntropy ' )
    print()

    total_units = float(array_length*array_layers)
#    printConfigVars (config_vars_list)      
                        
    x1 = float(config_vars_list[0])/total_units
    x2 = float(config_vars_list[1])/total_units

    sumX = x1+x2  
    print() 
    print( '                x1 = %.4f'  % x1 )
    print( '                x2 = %.4f'  % x2 )
    print( '   Sum of the x(i) = %.4f'  % sumX )
    print() 
    
    y1 = float(config_vars_list[2])/(2*total_units)
    y2 = float(config_vars_list[3])/(2*total_units)
    y3 = float(config_vars_list[4])/(2*total_units) 
    w1 = float(config_vars_list[5])/(2*total_units)
    w2 = float(config_vars_list[6])/(2*total_units)
    w3 = float(config_vars_list[7])/(2*total_units)        
    z1 = float(config_vars_list[8])/(2*total_units)
    z2 = float(config_vars_list[9])/(2*total_units)
    z3 = float(config_vars_list[10])/(2*total_units) 
    z4 = float(config_vars_list[11])/(2*total_units)
    z5 = float(config_vars_list[12])/(2*total_units)
    z6 = float(config_vars_list[13])/(2*total_units) 

    sumY = y1+y2+y3    
    print()
    print( "Totals for the y(i) variables:" )
    print( '          (A-A) y1 = %.4f'  % y1 )
    print( '    (A-B & B-A) y2 = %.4f'  % y2 )  
    print( '          (B-B) y3 = %.4f'  % y3 )
    print( '   Sum of the y(i) = %.4f'  % sumY )
    print() 
    
    sumW = w1+w2+w3 
    print()
    print( "Totals for the w(i) variables:" )
    print( '        (A---A) w1 = %.4f'  % w1 )
    print( '(A---B & B---A) w2 = %.4f'  % w2 )  
    print( '        (B---B) w3 = %.4f'  % w3 )
    print( '   Sum of the w(i) = %.4f'  % sumW )
    print()         

    sumZ = z1+z2+z3 +z4+z5+z6
    print()
    print( "Totals for the z(i) variables:" )
    print( '        (A-A-A) z1 = %.4f'  % z1 )
    print( '(A-A-B & B-A-A) z2 = %.4f'  % z2 )  
    print( '        (A-B-A) z3 = %.4f'  % z3 )
    print( '        (B-A-B) z4 = %.4f'  % z4 )
    print( '(A-B-B & B-B-A) z5 = %.4f'  % z5 )
    print( '        (B-B-B) z6 = %.4f'  % z6 )        
    print( '   Sum of the z(i) = %.4f'  % sumZ )
    print()                    
                                                
    Lfx = LfFunc(x1)+LfFunc(x2)
    print( ' lnFunc (x) = %.4f'  % (Lfx) )

    Lfy = LfFunc(y1) + 2*LfFunc(y2) + LfFunc(y3)
    print( ' lnFunc (y) = %.4f'  % (Lfy) ) 
      
    Lfw = LfFunc(w1) + 2*LfFunc(w2) + LfFunc(w3)
    print( ' lnFunc (w) = %.4f'  % (Lfw))

    Lfz = LfFunc(z1) + 2*LfFunc(z2) + LfFunc(z3) + LfFunc(z4) + 2*LfFunc(z5) + LfFunc(z6) 
    print( ' lnFunc (z) = %.4f'  % (Lfz))
        
    s = -(2*Lfy+Lfw-Lfx-2*Lfz)

    
    print()
    print( ' The entropy is: %.4f' % (s) ) 
  
   
     
    return (s)

    
####################################################################################################
####################################################################################################
#
# Function to compute the entropy, enthalpy, and free energy for a smooth calculation through 
#    a stepped range of h values
#
####################################################################################################
####################################################################################################


def compute_thermodynamic_vars(h, config_vars_list):

    total_units = float(array_length*array_layers)
    total_units_times_two = total_units*2.0
                       
    x1 = float(config_vars_list[0])/total_units
    x2 = float(config_vars_list[1])/total_units 
    
    y1 = float(config_vars_list[2])/total_units_times_two
    y2 = float(config_vars_list[3])/total_units_times_two
    y3 = float(config_vars_list[4])/total_units_times_two 
    
    sumY = y1 + y2 + y3

    w1 = float(config_vars_list[5])/total_units_times_two
    w2 = float(config_vars_list[6])/total_units_times_two
    w3 = float(config_vars_list[7])/total_units_times_two        

    sumW = w1 + w2 + w3

    z1 = float(config_vars_list[8])/total_units_times_two
    z2 = float(config_vars_list[9])/total_units_times_two
    z3 = float(config_vars_list[10])/total_units_times_two 
    z4 = float(config_vars_list[11])/total_units_times_two
    z5 = float(config_vars_list[12])/total_units_times_two
    z6 = float(config_vars_list[13])/total_units_times_two 
              

    sumZ = z1 + z2 + z3 + z4 + z5 + z6

    print()         
    print()
    print( ' x1 = %.4f'  % (x1), '      x2 = %.4f'  % (x2))
    print( ' y1 = %.4f'  % (y1), ' y2Total = %.4f'  % (2.*y2), ' y3 = %.4f'  % (y3) ) #, ' sumY = %.4f'  % (sumY)                                                 
    print( ' w1 = %.4f'  % (w1), ' w2Total = %.4f'  % (2.0*w2), ' w3 = %.4f'  % (w3) ) #, ' sumW = %.4f'  % (sumW) 

    print()
    print( ' z1 = %.4f'  % (z1), ' z2Total = %.4f'  % (2.0*z2), ' z3 = %.4f'  % (z3) )
    print( ' z6 = %.4f'  % (z6), ' z5Total = %.4f'  % (2.0*z5), ' z4 = %.4f'  % (z4) ) #, 'sumZ = %.4f'  % (sumW)      
        
                                                                                                                                                    
    Lfx = LfFunc(x1)+LfFunc(x2)
#    print ' lnFunc (x) = %.4f'  % (Lfx)

    Lfy = LfFunc(y1) + 2.0*LfFunc(y2) + LfFunc(y3)
#    print ' lnFunc (y) = %.4f'  % (Lfy)  
      
    Lfw = LfFunc(w1) + 2*LfFunc(w2) + LfFunc(w3)
#    print ' lnFunc (w) = %.4f'  % (Lfw)

    Lfz = LfFunc(z1) + 2*LfFunc(z2) + LfFunc(z3) + LfFunc(z4) + 2*LfFunc(z5) + LfFunc(z6)
#    print ' lnFunc (z) = %.4f'  % (Lfz)
        
    negS = -(2*Lfy+Lfw-Lfx-2*Lfz)

    


#   eps1 = 4.0*log(h)
#    fEps1 = eps1*(2.0*y2)
#    fEps0 = eps0*x2
#    fEnergy = fEps0 + fEps1 + negS 

        
#    print '  For h =', h, 's= %.4f' % (negS), 'fEps0= %.4f' % (fEps0),   'fEps1= %.4f' % (fEps1), 'fEnergy= %.4f' % (fEnergy)

    
    h = ((z1*y2)/(z3*y1))**(0.5)
    
 
    
    epsilon1 =  4*log(h)         
    


    enthalpy0 = 0.0
    enthalpy1 = 2.*epsilon1*y2      
    
    freeEnergy = enthalpy0 + enthalpy1 + negS          

        
    sysValsList = (negS, enthalpy0, enthalpy1, freeEnergy)

    print()  
    print() 
    print(' The computed thermodynamic quantities:' )           
    print( '    h = %.4f' % (h), '    epsilon1 = %.4f' % (epsilon1) )
    print( '    negative entropy s = %.4f' % (negS) )
    print( '    enthalpy = %.4f' % (enthalpy1) )     
    print( '    free energy = %.4f' % (freeEnergy) )
    print() 
    print()  
                                  
    return (sysValsList)   
         
    

        
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
#    (3.1) computeConfigXVariables
#    (3.2) computeConfigYVariables
#    (3.3) computeConfigWVariables
#    (3.4) computeConfigZVariables
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
#    (3.1) computeConfigXVariables: Computes the total values for X1 & X2 in the grid, going row-by-row.
#        NOTE: Unlike all the other configuration variable computations, requiring a wrap-around to the
#        units at the beginning of the row(s), this computation is a simple count-up of "A" and "B" units. 
#
#    (3.2) computeConfigYVariables: Computes the total values for Y1, Y2, and Y3 in the grid.
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
#    (3.4) computeConfigZVariables: Computes the total values for Z1, Z2, Z3, Z4, Z5, and Z6 in the grid.
#
#
####################################################################################################
####################################################################################################

def main():

####################################################################################################
# Obtain unit array size in terms of array_length (M) and layers (N)
####################################################################################################                

    global debug_print_off
    global detailed_debug_print_off
    global z_debug_print_off
    global blnkspc
    global array_length
    global array_layers
    global even_layers
    global pairs
    even_layers = True

    blnkspc=' '

# Define values for the global debug variables
    debug_print_off = True
    detailed_debug_print_off = True
    z_debug_print_off = True

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

    config_vars_list    = list()    # empty list
    sys_vals_list       = list()    # empty list
    
####################################################################################################
# Parmaters needed to compute configuration variables
####################################################################################################      
                       
    eps0 = 0.0  # The enthalpy for single unit activation is set to zero for this code
          


####################################################################################################
# User-definable variables influencing the numbers of steps, ranges, etc. to correct the
#  randomly-generated distribution so that it is close to the desired value for x1.  
####################################################################################################                 

    x1_target = 0.5     # Initially defining the target value for x1 as 0.5
    max_x_dif = 0.015   # Maximum difference between x1 from randomly-generated array and the
                        #   equiprobable distribution value of x1 = x1_target
                        #   Typical values are in the range 0.01 - 0.015.
    jrange = 3          # Maximal number of steps allowed to improve the x1 distribution
                        #   Typical values are in the range 100 - 200.             

####################################################################################################
# Compute configuration variables
####################################################################################################                

     

# Determine if we are probabilistically generating a pattern (patternProb = 0)
#   or if we are selecting a pre-stored pattern (patternProb = 1 ... N)

    pattern_select = obtain_pattern_selection()
    print_pattern_selection (pattern_select)
    
    
# Pick a starting value for h; it should be less than 1; it should be in the realm of 0.7 - 0.8    
    h0 = 0.8
# Pick the number of steps (increasing 0.01) for increasing h
    h_range = 3
# Pick the increment for increasing h
    h_incr = 0.1 
       
    print_run_parameters (h0, h_incr, h_range)             
 
    x_array = np.zeros(h_range, dtype=np.float)
    neg_S_array = np.zeros(h_range, dtype=np.float)
    f_eps0_array = np.zeros(h_range, dtype=np.float)         
    f_eps1_array = np.zeros(h_range, dtype=np.float)  
    f_energy_array = np.zeros(h_range, dtype=np.float) 
         

    
    h = h0
    for i in range (0, h_range, 1):        
        h = h + h_incr
        config_vars_list = list ()  # redefine this as an empty list
        unit_array    = initialize_matrix (array_size_list, pattern_select, h)        

# adjust the unit array so that x1 approximately = x2
        unit_array = adjust_matrix (array_size_list, unit_array, h, jrange, max_x_dif, x1_target)

# obtain the configuration variables (this step was also done while adjusting the array)       
        config_vars_list = compute_config_variables (array_size_list, unit_array, h)
                                                                
# obtain the thermodynamic variables 
        sys_vals_list = compute_thermodynamic_vars (h, config_vars_list)   

# store the thermodynamic variables to plot later
        x_array[i]=h
        neg_S_array[i]    = sys_vals_list[0] 
        f_eps0_array[i]   = sys_vals_list[1]
        f_eps1_array[i]   = sys_vals_list[2]-0.6  
        f_energy_array[i] = sys_vals_list[3]-0.6                                                
                                                                                                                              

    print() 
    print( ' The thermodynamic quantities plot vs. h' )
    print( '   The negEntropy is in blue, ' )
    print( '   The per-unit enthalpy is zero, and is not shown, ' )
    print( '   The interaction enthalpy (eps1*y2) is in maroon, and is shifted by 0.6,' ) 
    print( '     and the free energy is in red, also shifted by 0.6.'  ) 
    pylab.figure(1)
    pylab.plot (x_array,neg_S_array)    
    pylab.plot (x_array,f_eps1_array,'m')
    pylab.plot (x_array,f_energy_array,'r')
    pylab.show()  
 

          
####################################################################################################
# Conclude specification of the MAIN procedure
####################################################################################################                
    
if __name__ == "__main__": main()

####################################################################################################
# End program
####################################################################################################  