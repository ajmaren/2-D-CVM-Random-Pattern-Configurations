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


    print
    print '******************************************************************************'
    print
    print 'Welcome to the 2-D Cluster Variation Method'
    print 'Version 1.1, 12/08/2017, A.J. Maren'
    print 'This version computes thermodynamic quantities for equiprobable 2-D arrays'
    print 'By changing parameters in the main code, the user can select:'
    print '  (O) Randomly generating (and then improving) an array, or'
    print '  (1 .. N) Selecting a pre-stored array'
    print 'For comments, questions, or bug-fixes, contact: alianna.maren@northwestern.edu'
    print 'Alternate email address: alianna@aliannajmaren.com'
    print ' '
    print '  NOTE: In these calculations, x1 = A = 1 (units are at value 1),'
    print '                           and x2 = B = 0 (units are at value 0).'
    print
    print '******************************************************************************'
    print
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

def obtainArraySizeSpecs ():
    
#    x = input('Enter arraylength: ')
#    arraylength = int(x)
#    print 'arraylength is', arraylength  
          
#    x = input('Enter layers: ')
#    layers = int(x)
#    print 'layers is', layers
 

    arraylength = 16
    layers = 16
            
                
    arraySizeList = (arraylength, layers)  
    return (arraySizeList)  

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
# Function to obtain a row of 2-D CVM data - PRESTORED pattern (currently part of a 16x16 grid)
#
####################################################################################################
####################################################################################################


def obtainGridRow (rowNum, patternProb, h):


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
    if patternProb == 2:
        rowArray0 =   [1,0,0,0,1,0,0,1, 0,0,1,1,1,1,0,0] # Row 0 - top row 
        rowArray1 =   [1,0,1,0,0,0,1,1, 0,1,1,1,0,0,1,0] # Row 1 - second row (counting down from the top)
        rowArray2 =   [0,0,1,1,0,1,0,1, 1,0,1,0,0,1,1,1] # Row 2 - third row (counting down from the top) 
        rowArray3 =   [1,0,1,0,1,1,0,0, 1,0,0,0,1,1,1,0] # Row 3 - fourth row (counting down from the top)
        rowArray4 =   [1,1,0,0,0,1,0,1, 0,0,0,1,1,1,1,0] # Row 4 - fifth row (counting down from the top)
        rowArray5 =   [1,1,0,1,0,0,1,1, 1,0,0,1,1,1,0,0] # Row 5 - sixth row (counting down from the top)
        rowArray6 =   [1,1,0,1,0,0,0,1, 1,0,1,0,1,1,1,0] # Row 6 - seventh row (counting down from the top) 
        rowArray7 =   [1,0,1,1,0,1,0,0, 0,1,1,0,0,0,0,0] # Row 7 - eighth row (counting down from the top)  
    
# Rich Club
    if patternProb == 2:
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
####################################################################################################
#
# Procedure to randomly-generate an array, and then permute it to achieve the desired 
#    z1 & z3 values.
#
#    Inputs:    arraySizeList: a list of two integers; arrayLength and layers
#               h: the interaction enthalpy parameter
#    Return: the matrix unit_array, a matrix of 0's and 1's.
#
####################################################################################################
####################################################################################################

def initializeGeneratedMatrix (arraySizeList,h):

    localArrayLength = arraySizeList[0]
    localArrayLayers = arraySizeList[1]

    hSquared = h*h
    hFourth  = hSquared*hSquared
    denom    = 8.0*(hFourth - 6.0*hSquared + 1.0)  
    z3Analytic = (hSquared - 3.0)*(hSquared + 1.0)/denom      #  App. B, Eqn. 29
    z1Analytic = (1.0 - 3.0*hSquared)*(hSquared + 1.0)/denom  #  App. B, Eqn. 30
    y1Analytic = z1Analytic + 0.5*(0.5 - z1Analytic - z3Analytic)
    y2Analytic = z3Analytic + 0.5*(0.5 - z1Analytic - z3Analytic)

# Create the matrix 'unit_array' so that it has a random population of 0's and 1's.
    unit_array = np.random.choice([0, 1],size=(localArrayLayers,localArrayLength)) # Create an array filled with random values
# Note: this function can be used to create proportional distributions: np.random.choice([0, 1], size=(10,), p=[1./3, 2./3])


    print 'A randomly-generated array:'
    print
    
# Bring the array closer to the desired configuration variable values    
    
    return unit_array

####################################################################################################
####################################################################################################
#
# Procedure to initialize the matrix with EITHER a pre-stored pattern of values
#    OR randomly-generate an array, and then permute it to achieve the desired 
#    z1 & z3 values.
#
#    Inputs:    arraySizeList: a list of two integers; arrayLength and layers
#               patternProb: an integer indicating whether to randomly-generate
#                   and then permute an array (0), or select a pattern (1 .. N)
#               h: the interaction enthalpy parameter
#    Return: the matrix unit_array, a matrix of 0's and 1's.
#
####################################################################################################
####################################################################################################

 
def initializeMatrix (arraySizeList,patternProb, h):

           
    localArrayLength = arraySizeList[0]
    localArrayLayers = arraySizeList[1]





# Note: The passed value patternProb is used to determine if we are returning a stored pattern, or
#       are probabilistically-generating our data.
#       If patternProb = 0: probabilistic generation, dependent on h 
#       If patternProb > 1: select one of the N stored patterns (1 ... N)
        
    if patternProb == 0:
        unit_array = initializeGeneratedMatrix (arraySizeList,h)


# Create the initial matrix, 'unit_array,' and populate it with zeros
    else:     
        unit_array = np.zeros((localArrayLayers,localArrayLength), dtype=np.int)
      
                  
# Read the stored grid into unit_array
        x1_total = x2_total = 0   
        for i in range(0,localArrayLayers):
            dataArray = obtainGridRow (i, patternProb)
#            rownum = i+1
            for j in range(0, localArrayLength):
                unit_array[i,j]=dataArray[j]
                if unit_array[i,j]==1: 
                    x1_total = x1_total + 1
                else: x2_total = x2_total + 1 
             
# Print the unit_array so that it shows as a zigzag chain (or in the 2-D case, layers of zigzag chains).

    print 
    print 'The L x M array of units, where M (across) =', localArrayLength, 'and L (layers) =', localArrayLayers
    print 

# Determining if the grid has an even or odd number of layers is done as a global value determination in __main__

# Determine if we have an even or odd number of layers

# Determining "pairs" - the total number of pairs of zigzag chains - is done in __main__; "pairs" is a global variable

        
    for i in range (0, pairs):
        actualEvenRowNum = 2*i
        print 'Row', actualEvenRowNum, ':', blnkspc, 
        for j in range(0,localArrayLength):
            print unit_array[actualEvenRowNum,j], blnkspc,
        print 
        actualOddRowNum = 2*i+1
        print 'Row ', actualOddRowNum, ':', blnkspc,
        print (blnkspc),
        for j in range(0,localArrayLength):
            print unit_array[actualOddRowNum,j], blnkspc,
        print 
    print ' '
    
# Still need to write the print for an odd number of rows in grid


    
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


def computeConfigXVariables (arraySizeList, unitArray):


####################################################################################################
# This section unpacks the input variable arraySizeList
####################################################################################################

    arrayLength = arraySizeList [0]
    arrayLayers = arraySizeList [1]
    unit_array = unitArray

# Debug print statements
    if not debugPrintOff:
        print ' '   
        print "Just entered computeConfigXVariables"
                                   
# Initialize the y'i variables
    x1_total = x2_total = 0 

    
    for i in range (0,arrayLayers):
  
        x1_partial = x2_partial = 0  

    # Compute the x'i values for each sub-row of the zigzag, just to see 
    #   the distribution 
    # Start counting through the array elements, L->R.
        for j in range(0, arrayLength):
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

# Print the locally-computed values for x1 and x2; these are not passed back to Main.     
    print
    print 'The distribution among states A and B (x1 and x2) units is:'    
    print "  ( A ) x1_total =", x1_total, "( B ) x2_total =", x2_total  
    print '  For a total of ', x1_total + x2_total, ' units.'
    print ' ' 

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



def computeConfigYEvenRowZigzagVariables (arraySizeList, unitArray, topRow):


####################################################################################################
# This section unpacks the input variable arraySizeList
####################################################################################################

    arrayLength = arraySizeList [0]
    arrayLayers = arraySizeList [1]
    unit_array = unitArray
                                 

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
#for i in range(0,arrayLayers-1):
 #       top_row = i
  #      next_row = i+1
    top_row = topRow
    next_row = topRow + 1
  

# Start counting through the array elements, L->R.
    for j in range(0, arrayLength):
        # If the initial unit is A:
        if unit_array[top_row,j]>0.1: 
            # Compare with the same (jth) unit in the overlapping row 
            # comprising the zigzag chain
            # If the nearest-neighbor unit is also A:
            if unit_array[next_row,j] > 0.1:
                # Increment the y_1; the count of A-A nearest-neighbor pairs:
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
#for i in range(0,arrayLayers-1):
 #       top_row = i
  #      next_row = i+1
    

# Start counting through the array elements, L->R.
# Since we are comparing the unit in the lower row to the one shifted diagonally
# above and over to the right, we only step through to the arraylength - 1 unit.
# A final step (after this) will be to compute the wrap-around. 
    for j in range(0, arrayLength-1):
        # If the initial unit is A:
        if unit_array[next_row,j]>0.1: 
            # Compare with the NEXT (j+1) unit in the overlapping top row 
            # comprising the zigzag chain
            # If the nearest-neighbor unit is also A:
            if unit_array[top_row,j+1] > 0.1:
                # Increment the y_1; the count of A-A nearest-neighbor pairs:
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

    if unit_array[next_row,arrayLength-1]>0.1: 
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
    if not debugPrintOff:
        print
        print "Totals for the y'i variables:"
        print "(A-A) y1_total =", y1_total, "(A-B) left_y2_total =", left_y2_total   
        print "(B-B) y3_total =", y3_total, "(B-A) right_y2_total =", right_y2_total 
        print        

################################################################


###################################################################################################
#
# Assign the computed configuration variables to elements of the configVarsList, 
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



def computeConfigYOddRowZigzagVariables (arraySizeList, unitArray, topRow):


####################################################################################################
# This section unpacks the input variable arraySizeList
####################################################################################################

    arrayLength = arraySizeList [0]
    arrayLayers = arraySizeList [1]
    unit_array = unitArray

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
#for i in range(0,arrayLayers-1):
 #       top_row = i
  #      next_row = i+1
    top_row = topRow
    next_row = topRow + 1
    if top_row == arrayLayers-1: next_row = 0  
  

# Start counting through the array elements, L->R.
    for j in range(0, arrayLength-1): # Same logic as in the Even Row y(i) computation
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
#for i in range(0,arrayLayers-1):
 #       top_row = i
  #      next_row = i+1
    

# Start counting through the array elements, L->R.
# Since we are comparing the unit in the lower row to the one shifted diagonally
# above and over to the right, we only step through to the arraylength - 1 unit.
# A final step (after this) will be to compute the wrap-around. 
    for j in range(0, arrayLength): # Same logic as in the Even Row y(i) computation
            # But we can include the full array length (the other was truncated at arrayLength - 1)
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

    if unit_array[top_row,arrayLength-1]>0.1: 
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
    if not debugPrintOff:
        print
        print "Totals for the y'i variables:"
        print "(A-A) y1_total =", y1_total, "(A-B) left_y2_total =", left_y2_total   
        print "(B-B) y3_total =", y3_total, "(B-A) right_y2_total =", right_y2_total 
        print        

################################################################


###################################################################################################
#
# Assign the computed configuration variables to elements of the configVarsList, 
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

def computeConfigYVariables (arraySizeList, unitArray):

# Initialize the y'i variables

    y1 = y2 = y3 = 0

    if not debugPrintOff:
        print ' '
        print '  Starting to compute Y variables'
        print '  Total number of pairs of zigzag chains is: ', pairs
        print ' ' 
    for i in range (0, pairs): 
        topRow = 2*i
        if not debugPrintOff:
            print '  Row: ', topRow
        # Obtain the y(i) values from the first even-to-odd zigzag chain (0 to 1, running top-to-bottom)
        configVarsYList = computeConfigYEvenRowZigzagVariables (arraySizeList, unitArray, topRow) 
        # Assign the returned results to the local sum for each of the z(i) triplets
        y1 = y1+configVarsYList[0]
        y2 = y2+configVarsYList[1]
        y3 = y3+configVarsYList[2]
        topRow = 2*i+1
        if not debugPrintOff:
            print ' '
            print '  Row: ', topRow
        configVarsYList = computeConfigYOddRowZigzagVariables (arraySizeList, unitArray, topRow) 
        # Assign the returned results to the local sum for each of the z(i) triplets
        y1 = y1+configVarsYList[0]
        y2 = y2+configVarsYList[1]
        y3 = y3+configVarsYList[2]



# Debug section: Print totals for right-downwards-then-upwards triplets
        if not debugPrintOff:
            print ' '
            print ' -----------'
            print ' '
            print "Totals for all y(i), after completing Row: ", topRow
            print "           (A-A) z1_total =", y1
            print "(A-B) plus (B-A) z2_total =", y2
            print "           (B-B) z3_total =", y3   
            print ' '
            print ' -----------'
            print ' '
                
        # Start working on the next zigzag chain        
#        topRow = 2*i+1
        # Obtain the z(i) values from the first odd-to-even-to zigzag chain (1 to 2, running top-to-bottom)

#        configVarsZListOddToEven = computeConfigZVariablesOddToEven (arraySizeList, unitArray, topRow)                    
                                        
        # Add the returned results to the local sum for each of the z(i) triplets
#        z1 = z1+configVarsZListOddToEven[0]
#        z2 = z2+configVarsZListOddToEven[1]
#        z3 = z3+configVarsZListOddToEven[2]
#        z4 = z4+configVarsZListOddToEven[3]
#        z5 = z5+configVarsZListOddToEven[4]
#        z6 = z6+configVarsZListOddToEven[5]


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

def computeConfigWHorizontalRowVariables (arraySizeList, unitArray):


####################################################################################################
# This section unpacks the input variable arraySizeList
####################################################################################################

    arrayLength = arraySizeList [0]
    arrayLayers = arraySizeList [1]
    unit_array = unitArray

# Debug print statements
    if not debugPrintOff:
        print "Just entered computeConfigWVariables"
                                   
# Initialize the w'i variables
    w1_total = w2_total = w3_total = 0 
    w1_partial = w2_partial = w3_partial = 0  
    
    for i in range (0,arrayLayers):
        w1_partial = w2_partial = w3_partial = 0  
# Compute the w'i values for each sub-row of the zigzag, just to see 
#   the distribution 
# Start counting through the array elements, L->R.
        for j in range(0, arrayLength):             
            nextNearestNeighbor = j+1
            rowLimit = arrayLength-1
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

        detailedDebugPrintOffW = True
        if not detailedDebugPrintOffW:
            print ' '
            print "In row ", i
            print "w1_partial = ", w1_partial
            print "w2_partial = ", w2_partial                 
            print "w3_partial = ", w3_partial 
                                                        
        # Check the wrap-around value between the last unit in the row
        #   and the first item of this same row
#        if unit_array[i,arrayLength-1]>0.1: 
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
    if not debugPrintOff:
        print ' '
        print "Totals for all horizontal w(i)"
        print "            (A--A) w1_total =", w1
        print "(A--B) plus (B--A) w2_total =", w2
        print "            (B--B) w3_total =", w3   
        print ' '          

################################################################
                                   
                                                                      
                                                                                                                                            
    return (configVarsWList)


####################################################################################################
####################################################################################################
#
# Procedure to compute vertical configuration variables w'i and return as elements of list configWVarsList
#
####################################################################################################
####################################################################################################

def computeConfigWVerticalColVariables (arraySizeList, unitArray):


####################################################################################################
# This section unpacks the input variable arraySizeList
####################################################################################################

    arrayLength = arraySizeList [0]
    arrayLayers = arraySizeList [1]
    iLimit = arrayLayers - 2
    
# Debug print statements
    if not debugPrintOff:
        print "Just entered computeConfigWVerticalVariables"
                                   
# Initialize the w'i variables
    w1_total = w2_total = w3_total = 0 
    w1_partial = w2_partial = w3_partial = 0  

#    

    for i in range (0, arrayLayers): # run through the rows, look at those two rows apart
        vertPair = i+2        
        if i == iLimit: vertPair = 0
        if i == iLimit+1: vertPair = 1
        for j in range(0, arrayLength):             
            # If the initial unit is A:
            if unitArray[i,j]>0.1: 
                # The unit is "A," see if the next unit is "A" or "B" 
                if unitArray[vertPair,j]>0.1: 
                # Compare with the NEXT (j+1) unit in the SAME row 
                #   comprising a partial row of the the zigzag chain
                #   If this next-nearest-neighbor unit is also "A":                
                    w1_partial = w1_partial + 1
                else: # The next-nearest-neighbor is in "B"
                    w2_partial = w2_partial + 1                    
                                    
            else: # The initial unit is B:
                # The unit is "B," see if the next unit is "A" or "B" 
                if unitArray[vertPair,j]>0.1: 
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
    if not debugPrintOff:
        print ' '
        print "Totals for all vertical w(i)"
        print "            (A--A) w1_total =", w1
        print "(A--B) plus (B--A) w2_total =", w2
        print "            (B--B) w3_total =", w3   
        print ' '          

################################################################
            
                        
    return (configVarsWList)






####################################################################################################
####################################################################################################

# This function runs both the even-to-odd and odd-to-even zigzags; it is the first step in building
#   another row on top of the basic 1-D zigzag chain

####################################################################################################

def computeConfigWVariables (arraySizeList, unitArray):

# Initialize the y'i variables

    w1 = w2 = w3 = 0

    if not debugPrintOff:
        print ' '
        print '  Starting to compute W variables'
        print '  Total number of zigzag chains is: ', arrayLayers
        print ' ' 
    
    configVarsWList = computeConfigWHorizontalRowVariables (arraySizeList, unitArray) 
    # Assign the returned results to the local sum for each of the z(i) triplets
    w1 = w1+configVarsWList[0]
    w2 = w2+configVarsWList[1]
    w3 = w3+configVarsWList[2]


# Debug section: Print totals for right-downwards-then-upwards triplets
    if not debugPrintOff:
        print ' '
        print '  Row: ', arrayLayers
        print ' '
        print ' -----------'
        print ' '
        print "Totals for all horizontal w(i), after completing Row: ", arrayLayers
        print "            (A--A) w1_total =", w1
        print "(A--B) plus (B--A) w2_total =", w2
        print "            (B--B) w3_total =", w3   
        print ' '
        print ' -----------'
        print ' '

  
# NOTE: Still need to write the computation for an extra odd row in grid, if it exists


    configVarsWList = computeConfigWVerticalColVariables (arraySizeList, unitArray)
    w1 = w1+configVarsWList[0]
    w2 = w2+configVarsWList[1]
    w3 = w3+configVarsWList[2]       

# Debug section: Print totals for right-downwards-then-upwards triplets
    if not debugPrintOff: 
        print ' '
        print ' -----------'
        print ' '
        print "Totals for all horizontal and vertical w(i)"
        print "            (A--A) w1_total =", w1
        print "(A--B) plus (B--A) w2_total =", w2
        print "            (B--B) w3_total =", w3   
        print ' '
        print ' -----------'
        print ' '

    configVarsWList = (w1, w2, w3)                                                                                                                                                                        
    return (configVarsWList)





####################################################################################################
####################################################################################################
#
# Procedure to compute the the precise value of a triplet z'i variable given
#   Input: integer values for unit (U), nearest-neighbor (NN), next-nearest-neighbor (= NNN)
#   Output: locally-incremented values for z1 & z2 & z3 & z4 & z5 & z6 
#
####################################################################################################
####################################################################################################



def computeSpecificTripletZVariable (U, NN, NNN):

# Debug print statements
    if not detailedDebugPrintOff:
        print ' '
        print "Just entered computeSpecificTripletZVariable"
    
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

def debugPrintZ (z1_incr, left_z2_incr, right_z2_incr, z3_incr, z4_incr, left_z5_incr, right_z5_incr, z6_incr):
        
    print ' '
    print '(A-A-A) z1_incr =', z1_incr, '(A-A-B) left_z2_incr =', left_z2_incr 
    print '(A-B-A) z3_incr =', z3_incr, '(B-A-A) right_z2_incr =', right_z2_incr  
    print '(B-A-B) z4_incr =', z4_incr, '(B-B-A) left_z5_incr =', left_z5_incr 
    print '(B-B-B) z6_incr =', z6_incr, '(A-B-B) right_z5_incr =', right_z5_incr  
    print  ' '  
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



def computeConfigZEvenUpperToLower (arraySizeList, unitArray, top_row):

    arrayLength = arraySizeList [0]
    arrayLayers = arraySizeList [1]
    unit_array = unitArray
    
# Create the array to hold the partial (the increments in the) z'i's, and populate it with zeros
    zPartialArray = np.zeros((arrayLength), dtype=np.int)

  
    z1_partial = left_z2_partial = right_z2_partial = z3_partial = 0
    z4_partial = left_z5_partial = right_z5_partial = z6_partial = 0


    next_row = top_row + 1

# Start counting through the array elements, L->R.
    for j in range(0, arrayLength-1):
        U = unit_array[top_row,j]
        NN = unit_array[next_row,j]
        NNN = unit_array[top_row, j+1]

        TripletValueList = computeSpecificTripletZVariable (U, NN, NNN)

    

# Debug print statements
#        if not detailedDebugPrintOff:
#            print ' '
#            print 'Debug printing: computeConfigZEvenUpperToLower'  #debugPrintOff false        
#            print "Returning from computeSpecificTripletZVariable" 
#            print "Unpacking the specific triplet value found"
    
        z1_incr = TripletValueList[0]
        left_z2_incr = TripletValueList[1]
        right_z2_incr = TripletValueList[2]
        z3_incr = TripletValueList[3]
        z4_incr = TripletValueList[4] 
        left_z5_incr = TripletValueList[5]
        right_z5_incr = TripletValueList[6] 
        z6_incr = TripletValueList[7]
            
        # Debug print statements
        if not detailedDebugPrintOff:
            print ' '
            print 'Debug printing: computeConfigZEvenUpperToLower, incrementing the z(i) in for loop:'  #debugPrintOff false
            debugPrintZ (z1_incr, left_z2_incr, right_z2_incr, z3_incr, z4_incr, left_z5_incr, right_z5_incr, z6_incr)
    
                                    
        z1_partial = z1_partial + z1_incr 
        left_z2_partial = left_z2_partial + left_z2_incr
        right_z2_partial = right_z2_partial + right_z2_incr 
        z3_partial = z3_partial + z3_incr 
        z4_partial = z4_partial + z4_incr 
        left_z5_partial = left_z5_partial + left_z5_incr 
        right_z5_partial = right_z5_partial + right_z5_incr 
        z6_partial = z6_partial + z6_incr 
      
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



def computeConfigZEvenLowerToUpper (arraySizeList, unitArray, top_row):

    arrayLength = arraySizeList [0]
    arrayLayers = arraySizeList [1]
    unit_array = unitArray
    
# Create the array to hold the partial (the increments in the) z'i's, and populate it with zeros
    zPartialArray = np.zeros((arrayLength), dtype=np.int)

  
    z1_partial = left_z2_partial = right_z2_partial = z3_partial = 0
    z4_partial = left_z5_partial = right_z5_partial = z6_partial = 0

    next_row = top_row + 1


# NOTE: We are computing the SECOND row of triplets in a zigzag chain,
#   going from the lower row to the top
# Start counting through the array elements, L->R.
    for j in range(0, arrayLength-1):
        U = unit_array[next_row,j]
        NN = unit_array[top_row,j+1]
        NNN = unit_array[next_row, j+1]

        TripletValueList = computeSpecificTripletZVariable (U, NN, NNN)

    
    # Debug print statements
        if not detailedDebugPrintOff:
            print ' '
            print "Computing the SECOND row of a zigzag chain for leading unit ", j
            print "Returning from computeSpecificTripletZVariable" 
            print "Unpacking the specific triplet value found"
    
        z1_incr = TripletValueList[0]
        left_z2_incr = TripletValueList[1]
        right_z2_incr = TripletValueList[2]
        z3_incr = TripletValueList[3]
        z4_incr = TripletValueList[4] 
        left_z5_incr = TripletValueList[5]
        right_z5_incr = TripletValueList[6] 
        z6_incr = TripletValueList[7]

        # Debug print statements
        if not detailedDebugPrintOff:
            print ' '
            print 'Debug printing: computeConfigZEvenLowerToUpper, incrementing the z(i) in for loop:'  #debugPrintOff false
            debugPrintZ (z1_incr, left_z2_incr, right_z2_incr, z3_incr, z4_incr, left_z5_incr, right_z5_incr, z6_incr)
    
                                    
        z1_partial = z1_partial + z1_incr 
        left_z2_partial = left_z2_partial + left_z2_incr
        right_z2_partial = right_z2_partial + right_z2_incr 
        z3_partial = z3_partial + z3_incr 
        z4_partial = z4_partial + z4_incr 
        left_z5_partial = left_z5_partial + left_z5_incr 
        right_z5_partial = right_z5_partial + right_z5_incr 
        z6_partial = z6_partial + z6_incr 

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
    print ' *************************'
    print ' '     
    print 'top_row = ', top_row, ' next_row = ', next_row
    print 'Row', top_row, ':', blnkspc, 
    for j in range(0,arrayLength):
        print unit_array[top_row,j], blnkspc,
    print 

    print 'Row ', next_row, ':', blnkspc,
    print (blnkspc),
    for j in range(0,arrayLength):
        print unit_array[next_row,j], blnkspc,
    print 
    print ' *************************'
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



def computeConfigZVariablesEvenToOdd (arraySizeList, unitArray, topRow):


####################################################################################################
# This section unpacks the input variable arraySizeList
####################################################################################################

    arrayLength = arraySizeList [0]
    arrayLayers = arraySizeList [1]
    unit_array = unitArray

# Debug print statements
    if not detailedDebugPrintOff:
        print ' '
        print "Just entered computeConfigZVariables: Even-to-Odd"
                                  
                                                               

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
    if not detailedDebugPrintOff:
        print ' '
        print "Calling computeSpecificTripletZVariable"

    U  = NN = NNN = 0
    
    TripletValueList = list()  
 
# commenting out for debug   
#for i in range(0,arrayLayers-1):
 #       top_row = i
 #      next_row = i+1
    top_row = topRow
    next_row = topRow + 1  

    if not ZDebugPrintOff: 
        printEvenToOddRows (top_row, unit_array)       
      
    z1_partial = left_z2_partial = right_z2_partial = z3_partial = 0
    z4_partial = left_z5_partial = right_z5_partial = z6_partial = 0

# Compute the first contribution to z'i's from the top-to-bottom-to-top row
    zPartialArray = computeConfigZEvenUpperToLower (arraySizeList, unitArray, top_row) 
             
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
    if not ZDebugPrintOff:
        print ' '
        print "Subtotals so far (downward-right-then_upwards-right-pointing triplets)"
        print "Before any wrap-arounds:"
        print "(A-A-A) z1_total =", z1_total, "(A-A-B) left_z2_total =", left_z2_total 
        print "(A-B-A) z3_total =", z3_total, "(B-A-A) right_z2_total =", right_z2_total  
        print "(B-A-B) z4_total =", z4_total, "(B-B-A) left_z5_total =", left_z5_total 
        print "(B-B-B) z6_total =", z6_total, "(A-B-B) right_z5_total =", right_z5_total



###################################################################################################
#
# Compute the triplet values z(i) for the first wrap-around
# Start at last unit on top row, use downward-pointing-diagonal for last unit on next row
# Then wrap-around to pick up first unit on top row.
#
###################################################################################################

  

# Debug print statements
    if not detailedDebugPrintOff: 
        print ' '
        print "Computing first wrap-around triplet"
        print "Calling computeSpecificTripletZVariable"

    U  = NN = NNN = 0

    U = unit_array[top_row,arrayLength-1]
    NN = unit_array[next_row,arrayLength-1]
    NNN = unit_array[top_row, 0]

    TripletValueList = computeSpecificTripletZVariable (U, NN, NNN)

    
    # Debug print statements
    if not detailedDebugPrintOff:  
        print ' '
        print "Returning from computeSpecificTripletZVariable" 
        print "Unpacking the specific triplet value found"
    
    z1_incr = TripletValueList[0]
    left_z2_incr = TripletValueList[1]
    right_z2_incr = TripletValueList[2]
    z3_incr = TripletValueList[3]
    z4_incr = TripletValueList[4] 
    left_z5_incr = TripletValueList[5]
    right_z5_incr = TripletValueList[6] 
    z6_incr = TripletValueList[7]

    # Debug print statements

    if not debugPrintOff: 
        print ' '         
        print "(A-A-A) z1_incr =", z1_incr, "(A-A-B) left_z2_incr =", left_z2_incr 
        print "(A-B-A) z3_incr =", z3_incr, "(B-A-A) right_z2_incr =", right_z2_incr  
        print "(B-A-B) z4_incr =", z4_incr, "(B-B-A) left_z5_incr =", left_z5_incr 
        print "(B-B-B) z6_incr =", z6_incr, "(A-B-B) right_z5_incr =", right_z5_incr  
        print ' End of incremental z(i) for first wrap-around, even-to-odd, top row = ', top_row
        print ' '                                     
                                                
# Update the total z'i values:                                                                                                                                                                 
    z1_total = z1_total + z1_incr
    left_z2_total = left_z2_total + left_z2_incr 
    right_z2_total = right_z2_total + right_z2_incr      
    z3_total = z3_total + z3_incr
    z4_total = z4_total + z4_incr
    left_z5_total = left_z5_total + left_z5_incr
    right_z5_total = right_z5_total + right_z5_incr  
    z6_total = z6_total + z6_incr
    z5_total = left_z5_total + right_z5_total
    z2_total = left_z2_total + right_z2_total
    
# Debug section: Print totals for right-downwards-then-upwards triplets
    if not ZDebugPrintOff:
        print ' '
        print "Subtotals so far (downward-right-then_upwards-right-pointing triplets)"
        print "After adding in the top-layer wrap-around:"
        print "(A-A-A) z1_total =", z1_total, "(A-A-B) left_z2_total =", left_z2_total 
        print "(A-B-A) z3_total =", z3_total, "(B-A-A) right_z2_total =", right_z2_total  
        print "(B-A-B) z4_total =", z4_total, "(B-B-A) left_z5_total =", left_z5_total 
        print "(B-B-B) z6_total =", z6_total, "(A-B-B) right_z5_total =", right_z5_total


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
#for i in range(0,arrayLayers-1):
 #       top_row = i
  #      next_row = i+1
     
#    for i in range (0,1):
#        top_row = 0
#        next_row = 1  
  
    zPartialArray = computeConfigZEvenLowerToUpper (arraySizeList, unitArray, top_row)
             
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
    if not ZDebugPrintOff:
        print ' '
        print "Subtotals so far (adding in upwards-right-then_downwards-right-pointing triplets)"
        print "Before the last wrap-around:"
        print "(A-A-A) z1_total =", z1_total, "(A-A-B) left_z2_total =", left_z2_total 
        print "(A-B-A) z3_total =", z3_total, "(B-A-A) right_z2_total =", right_z2_total  
        print "(B-A-B) z4_total =", z4_total, "(B-B-A) left_z5_total =", left_z5_total 
        print "(B-B-B) z6_total =", z6_total, "(A-B-B) right_z5_total =", right_z5_total


                 
                                                
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
    if not ZDebugPrintOff:
        print "Computing second wrap-around triplet"
        print "Calling computeSpecificTripletZVariable"
    U  = NN = NNN = 0

    U = unit_array[next_row,arrayLength-1]
    NN = unit_array[top_row,0]
    NNN = unit_array[next_row, 0]

    TripletValueList = computeSpecificTripletZVariable (U, NN, NNN)

    
    # Debug print statements
#    print "Returning from computeSpecificTripletZVariable" 
#    print "Unpacking the specific triplet value found"
    
    z1_incr = TripletValueList[0]
    left_z2_incr = TripletValueList[1]
    right_z2_incr = TripletValueList[2]
    z3_incr = TripletValueList[3]
    z4_incr = TripletValueList[4] 
    left_z5_incr = TripletValueList[5]
    right_z5_incr = TripletValueList[6] 
    z6_incr = TripletValueList[7]

    # Debug print statements
    if not ZDebugPrintOff:
        print          
        print "(A-A-A) z1_incr =", z1_incr, "(A-A-B) left_z2_incr =", left_z2_incr 
        print "(A-B-A) z3_incr =", z3_incr, "(B-A-A) right_z2_incr =", right_z2_incr  
        print "(B-A-B) z4_incr =", z4_incr, "(B-B-A) left_z5_incr =", left_z5_incr 
        print "(B-B-B) z6_incr =", z6_incr, "(A-B-B) right_z5_incr =", right_z5_incr  
        print
                                    
                                                
# Update the total z'i values:                                                                                                                                                                 
    z1_total = z1_total + z1_incr
    left_z2_total = left_z2_total + left_z2_incr 
    right_z2_total = right_z2_total + right_z2_incr      
    z3_total = z3_total + z3_incr
    z4_total = z4_total + z4_incr
    left_z5_total = left_z5_total + left_z5_incr
    right_z5_total = right_z5_total + right_z5_incr  
    z6_total = z6_total + z6_incr
    z5_total = left_z5_total + right_z5_total
    z2_total = left_z2_total + right_z2_total
    
# Debug section: Print totals for right-downwards-then-upwards triplets
    if not ZDebugPrintOff:
        print "Totals for all triplets, after adding in the second wrap-around"
        print "(A-A-A) z1_total =", z1_total, "(A-A-B) left_z2_total =", left_z2_total 
        print "(A-B-A) z3_total =", z3_total, "(B-A-A) right_z2_total =", right_z2_total  
        print "(B-A-B) z4_total =", z4_total, "(B-B-A) left_z5_total =", left_z5_total 
        print "(B-B-B) z6_total =", z6_total, "(A-B-B) right_z5_total =", right_z5_total



# This concludes computation of the z'i totals for a complete pass through a zigzag chain


###################################################################################################
#
# Assign the computed configuration variables to elements of the configVarsList, 
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



def computeConfigZOddUpperToLower (arraySizeList, unitArray, top_row):

    arrayLength = arraySizeList [0]
    arrayLayers = arraySizeList [1]
    unit_array = unitArray
    
# Create the array to hold the partial (the increments in the) z'i's, and populate it with zeros
    zPartialArray = np.zeros((arrayLength), dtype=np.int)

  
    z1_partial = left_z2_partial = right_z2_partial = z3_partial = 0
    z4_partial = left_z5_partial = right_z5_partial = z6_partial = 0


    next_row = top_row + 1
    bottomRow = arrayLayers
    if next_row == bottomRow: next_row = 0


# Start counting through the array elements, L->R.
    for j in range(0, arrayLength-1):
        U = unit_array[top_row,j]
        NN = unit_array[next_row,j+1]
        NNN = unit_array[top_row, j+1]

        TripletValueList = computeSpecificTripletZVariable (U, NN, NNN)

    
    # Debug print statements
        # Debug print statements
        if not detailedDebugPrintOff:
            print ' '
            print 'Debug printing: computeConfigZEvenUpperToLower'  #debugPrintOff false        
            print "Returning from computeSpecificTripletZVariable" 
            print "Unpacking the specific triplet value found"
    
        z1_incr = TripletValueList[0]
        left_z2_incr = TripletValueList[1]
        right_z2_incr = TripletValueList[2]
        z3_incr = TripletValueList[3]
        z4_incr = TripletValueList[4] 
        left_z5_incr = TripletValueList[5]
        right_z5_incr = TripletValueList[6] 
        z6_incr = TripletValueList[7]
            
        # Debug print statements
        if not detailedDebugPrintOff:
            if j == 0:
                print ' '
                print 'Debug printing: computeConfigZOddUpperToLower, incrementing the z(i) in for loop:'  #debugPrintOff false
                print ' j = ', j, '   U = ', U, '   NN = ', NN, '   NNN = ', NNN
                debugPrintZ (z1_incr, left_z2_incr, right_z2_incr, z3_incr, z4_incr, left_z5_incr, right_z5_incr, z6_incr)
    
                                    
        z1_partial = z1_partial + z1_incr 
        left_z2_partial = left_z2_partial + left_z2_incr
        right_z2_partial = right_z2_partial + right_z2_incr 
        z3_partial = z3_partial + z3_incr 
        z4_partial = z4_partial + z4_incr 
        left_z5_partial = left_z5_partial + left_z5_incr 
        right_z5_partial = right_z5_partial + right_z5_incr 
        z6_partial = z6_partial + z6_incr 
      
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



def computeConfigZOddLowerToUpper (arraySizeList, unitArray, top_row):

    arrayLength = arraySizeList [0]
    arrayLayers = arraySizeList [1]
    unit_array = unitArray
    
# Create the array to hold the partial (the increments in the) z'i's, and populate it with zeros
    zPartialArray = np.zeros((arrayLength), dtype=np.int)

  
    z1_partial = left_z2_partial = right_z2_partial = z3_partial = 0
    z4_partial = left_z5_partial = right_z5_partial = z6_partial = 0


    next_row = top_row + 1
    bottomRow = arrayLayers
    if next_row == bottomRow: next_row = 0

    if not ZDebugPrintOff:
        if top_row == 0: 
            print ' ' 
            print ' In OddToEven'
            print ' top_row = ', top_row
            print ' bottomRow = ', bottomRow
            print ' next_row = ', next_row
            print ' '   
        
# NOTE: We are computing the SECOND row of triplets in a zigzag chain,
#   going from the lower row to the top
# Start counting through the array elements, L->R.
    for j in range(0, arrayLength-1):
        U = unit_array[next_row,j]
        NN = unit_array[top_row,j]
        NNN = unit_array[next_row, j+1]
        if not detailedDebugPrintOff:
            if top_row == 15:
                print ' For top_row = ', top_row, ' and j = ', j, ', then U = ', U, ' and next_row = ', next_row, ' and j+1 is ', j+1, ' and NN = ', NN 
        TripletValueList = computeSpecificTripletZVariable (U, NN, NNN)

    
  


    
    # Debug print statements
        if not detailedDebugPrintOff:
            print ' '
            print "Computing the SECOND row of a zigzag chain for leading unit ", j
            print "Returning from computeSpecificTripletZVariable" 
            print "Unpacking the specific triplet value found"
    
        z1_incr = TripletValueList[0]
        left_z2_incr = TripletValueList[1]
        right_z2_incr = TripletValueList[2]
        z3_incr = TripletValueList[3]
        z4_incr = TripletValueList[4] 
        left_z5_incr = TripletValueList[5]
        right_z5_incr = TripletValueList[6] 
        z6_incr = TripletValueList[7]

        # Debug print statements
        if not detailedDebugPrintOff:
            print ' '
            print 'Debug printing: computeConfigZEvenLowerToUpper, incrementing the z(i) in for loop:'  #debugPrintOff false
            debugPrintZ (z1_incr, left_z2_incr, right_z2_incr, z3_incr, z4_incr, left_z5_incr, right_z5_incr, z6_incr)
    
                                    
        z1_partial = z1_partial + z1_incr 
        left_z2_partial = left_z2_partial + left_z2_incr
        right_z2_partial = right_z2_partial + right_z2_incr 
        z3_partial = z3_partial + z3_incr 
        z4_partial = z4_partial + z4_incr 
        left_z5_partial = left_z5_partial + left_z5_incr 
        right_z5_partial = right_z5_partial + right_z5_incr 
        z6_partial = z6_partial + z6_incr 

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
    bottom_row = arrayLayers - 1
    if top_row == bottom_row: next_row = 0
    print ' *************************'
    print ' '     
    print 'top_row = ', top_row, ' next_row = ', next_row
    print 'Row', top_row, ':', blnkspc, 
    print (blnkspc),
    for j in range(0,arrayLength):
        print unit_array[top_row,j], blnkspc,
    print 

    print 'Row ', next_row, ':', blnkspc,
    for j in range(0,arrayLength):
        print unit_array[next_row,j], blnkspc,
    print 
    print ' *************************'
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



def computeConfigZVariablesOddToEven (arraySizeList, unitArray, topRow):


####################################################################################################
# This section unpacks the input variable arraySizeList
####################################################################################################

    arrayLength = arraySizeList [0]
    arrayLayers = arraySizeList [1]
    unit_array = unitArray

    top_row = topRow
    next_row = topRow + 1
    if top_row == arrayLayers-1: next_row = 0

# Debug print statements
    if not detailedDebugPrintOff:
        print ' '
        print "Just entered computeConfigZVariables: Odd-to-Even"

    if not ZDebugPrintOff:
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
    if not detailedDebugPrintOff:
        print ' '
        print "Calling computeSpecificTripletZVariable"

    U  = NN = NNN = 0
    
    TripletValueList = list()  
 
# commenting out for debug   
#for i in range(0,arrayLayers-1):
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
    zPartialArray = computeConfigZOddUpperToLower (arraySizeList, unitArray, top_row) 
             
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
    if not debugPrintOff:
        print ' '
        print "Subtotals so far (downward-right-then_upwards-right-pointing triplets)"
        print "Before any wrap-arounds:"
        print "(A-A-A) z1_total =", z1_total, "(A-A-B) left_z2_total =", left_z2_total 
        print "(A-B-A) z3_total =", z3_total, "(B-A-A) right_z2_total =", right_z2_total  
        print "(B-A-B) z4_total =", z4_total, "(B-B-A) left_z5_total =", left_z5_total 
        print "(B-B-B) z6_total =", z6_total, "(A-B-B) right_z5_total =", right_z5_total



###################################################################################################
#
# Compute the triplet values z(i) for the first wrap-around on an ODD-to-EVEN zigzag
# Start at last unit on top row, use downward-pointing-diagonal for last unit on next row
# Then wrap-around to pick up first unit on top row.
#
###################################################################################################

  

# Debug print statements
    if not detailedDebugPrintOff:
        print ' '
        print "Computing first wrap-around triplet"
        print "Calling computeSpecificTripletZVariable"

    U  = NN = NNN = 0

# Note that when we are working an ODD-to-EVEN zigzag chain, the
#   first wrap-around triplet works on different units than when 
#   we are working an EVEN-to-ODD chain. 

    if not detailedDebugPrintOff:
        print '  '
        print '  debug in first wrap-around triplet in an ODD-to-EVEN zigzag computation'



    U = unit_array[top_row,arrayLength-1]
    NN = unit_array[next_row,0]
    NNN = unit_array[top_row, 0]

    TripletValueList = computeSpecificTripletZVariable (U, NN, NNN)

    
    # Debug print statements
    if not detailedDebugPrintOff:
        print ' '
        print "Returning from computeSpecificTripletZVariable" 
        print "Unpacking the specific triplet value found"
    
    z1_incr = TripletValueList[0]
    left_z2_incr = TripletValueList[1]
    right_z2_incr = TripletValueList[2]
    z3_incr = TripletValueList[3]
    z4_incr = TripletValueList[4] 
    left_z5_incr = TripletValueList[5]
    right_z5_incr = TripletValueList[6] 
    z6_incr = TripletValueList[7]

    # Debug print statements
    if not debugPrintOff:
        print ' '         
        print "(A-A-A) z1_incr =", z1_incr, "(A-A-B) left_z2_incr =", left_z2_incr 
        print "(A-B-A) z3_incr =", z3_incr, "(B-A-A) right_z2_incr =", right_z2_incr  
        print "(B-A-B) z4_incr =", z4_incr, "(B-B-A) left_z5_incr =", left_z5_incr 
        print "(B-B-B) z6_incr =", z6_incr, "(A-B-B) right_z5_incr =", right_z5_incr  
        print
                                    
                                                
# Update the total z'i values:                                                                                                                                                                 
    z1_total = z1_total + z1_incr
    left_z2_total = left_z2_total + left_z2_incr 
    right_z2_total = right_z2_total + right_z2_incr      
    z3_total = z3_total + z3_incr
    z4_total = z4_total + z4_incr
    left_z5_total = left_z5_total + left_z5_incr
    right_z5_total = right_z5_total + right_z5_incr  
    z6_total = z6_total + z6_incr
    z5_total = left_z5_total + right_z5_total
    z2_total = left_z2_total + right_z2_total
    
# Debug section: Print totals for right-downwards-then-upwards triplets
    if not debugPrintOff:
        print ' '
        print "Subtotals so far (downward-right-then_upwards-right-pointing triplets)"
        print "After adding in the top-layer wrap-around:"
        print "(A-A-A) z1_total =", z1_total, "(A-A-B) left_z2_total =", left_z2_total 
        print "(A-B-A) z3_total =", z3_total, "(B-A-A) right_z2_total =", right_z2_total  
        print "(B-A-B) z4_total =", z4_total, "(B-B-A) left_z5_total =", left_z5_total 
        print "(B-B-B) z6_total =", z6_total, "(A-B-B) right_z5_total =", right_z5_total



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
    
    zPartialArray = computeConfigZOddLowerToUpper (arraySizeList, unitArray, top_row)
             
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
    if not ZDebugPrintOff:
        print "Subtotals so far (adding in upwards-right-then_downwards-right-pointing triplets)"
        print "Before the last wrap-around:"
        print "(A-A-A) z1_total =", z1_total, "(A-A-B) left_z2_total =", left_z2_total 
        print "(A-B-A) z3_total =", z3_total, "(B-A-A) right_z2_total =", right_z2_total  
        print "(B-A-B) z4_total =", z4_total, "(B-B-A) left_z5_total =", left_z5_total 
        print "(B-B-B) z6_total =", z6_total, "(A-B-B) right_z5_total =", right_z5_total


                 
                                                
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
    if not ZDebugPrintOff:
        print "Computing second wrap-around triplet"
        print "Calling computeSpecificTripletZVariable"
    U  = NN = NNN = 0

# Note that the triplet wraparound for the second row of an ODD-to-EVEN zigzag chain
#   uses different units as compared with this same triplet on an EVEN-to-ODD zigzag
    U = unit_array[next_row,arrayLength-1]
    NN = unit_array[top_row,arrayLength-1]
    NNN = unit_array[next_row, 0]

    TripletValueList = computeSpecificTripletZVariable (U, NN, NNN)

    
    # Debug print statements
    if not detailedDebugPrintOff:
        print "Returning from computeSpecificTripletZVariable" 
        print "Unpacking the specific triplet value found"
    
    z1_incr = TripletValueList[0]
    left_z2_incr = TripletValueList[1]
    right_z2_incr = TripletValueList[2]
    z3_incr = TripletValueList[3]
    z4_incr = TripletValueList[4] 
    left_z5_incr = TripletValueList[5]
    right_z5_incr = TripletValueList[6] 
    z6_incr = TripletValueList[7]

    # Debug print statements
    if not ZDebugPrintOff:
        print          
        print "(A-A-A) z1_incr =", z1_incr, "(A-A-B) left_z2_incr =", left_z2_incr 
        print "(A-B-A) z3_incr =", z3_incr, "(B-A-A) right_z2_incr =", right_z2_incr  
        print "(B-A-B) z4_incr =", z4_incr, "(B-B-A) left_z5_incr =", left_z5_incr 
        print "(B-B-B) z6_incr =", z6_incr, "(A-B-B) right_z5_incr =", right_z5_incr  
        print
                                    
                                                
# Update the total z'i values:                                                                                                                                                                 
    z1_total = z1_total + z1_incr
    left_z2_total = left_z2_total + left_z2_incr 
    right_z2_total = right_z2_total + right_z2_incr      
    z3_total = z3_total + z3_incr
    z4_total = z4_total + z4_incr
    left_z5_total = left_z5_total + left_z5_incr
    right_z5_total = right_z5_total + right_z5_incr  
    z6_total = z6_total + z6_incr
    z5_total = left_z5_total + right_z5_total
    z2_total = left_z2_total + right_z2_total
    
# Debug section: Print totals for right-downwards-then-upwards triplets
    if not ZDebugPrintOff:
        print "Totals for all triplets, after adding in the second wrap-around"
        print "(A-A-A) z1_total =", z1_total, "(A-A-B) left_z2_total =", left_z2_total 
        print "(A-B-A) z3_total =", z3_total, "(B-A-A) right_z2_total =", right_z2_total  
        print "(B-A-B) z4_total =", z4_total, "(B-B-A) left_z5_total =", left_z5_total 
        print "(B-B-B) z6_total =", z6_total, "(A-B-B) right_z5_total =", right_z5_total



# This concludes computation of the z'i totals for a complete pass through a zigzag chain



# Lots of steps need to happen next

####################################################################################################
#
# Assign the computed configuration variables to elements of the configVarsList, 
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

def computeConfigZVariables (arraySizeList, unitArray):

# Initialize the z'i variables

    z1 = z2 = z3 = z4 = z5 = z6 = 0

#    z1_total = left_z2_total = right_z2_total = z3_total = 0   
#    z4_total = left_z5_total = right_z5_total = z6_total = 0   
#    y1_total = left_y2_total = right_y2_total = y3_total = 0  

    if not debugPrintOff:
        print ' '
        print '  Starting to compute Z variables'
        print '  Total number of pairs of zigzag chains is: ', pairs
        print ' ' 
    for i in range (0, pairs):
        topRow = 2*i
        if not ZDebugPrintOff:
            print '  Row: ', topRow
        # Obtain the z(i) values from the first even-to-odd zigzag chain (0 to 1, running top-to-bottom)
        configVarsZListEvenToOdd = computeConfigZVariablesEvenToOdd (arraySizeList, unitArray, topRow)
        # Assign the returned results to the local sum for each of the z(i) triplets
        z1 = z1+configVarsZListEvenToOdd[0]
        z2 = z2+configVarsZListEvenToOdd[1]
        z3 = z3+configVarsZListEvenToOdd[2]
        z4 = z4+configVarsZListEvenToOdd[3]
        z5 = z5+configVarsZListEvenToOdd[4]
        z6 = z6+configVarsZListEvenToOdd[5]


# Debug section: Print totals for right-downwards-then-upwards triplets
        if not ZDebugPrintOff:
            print ' '
            print 'Starting for loop with i = ', i
            print ' -----------'
            print ' '
            print "Totals for all triplets, after completing Row: ", topRow, "downwards-to-upwards"
            print "             (A-A-A) z1_total =", z1
            print "(A-A-B) plus (B-A-A) z2_total =", z2
            print "             (A-B-A) z3_total =", z3   
            print "             (B-A-B) z4_total =", z4
            print "(B-B-A) plus (A-B-B) z5_total =", z5 
            print "             (B-B-B) z6_total =", z6
            print ' '
            print ' -----------'
            print ' '
                
        # Start working on the next zigzag chain        
        topRow = 2*i+1
        # Obtain the z(i) values from the first odd-to-even-to zigzag chain (1 to 2, running top-to-bottom)

        configVarsZListOddToEven = computeConfigZVariablesOddToEven (arraySizeList, unitArray, topRow)                    
                                        
        # Add the returned results to the local sum for each of the z(i) triplets
        z1 = z1+configVarsZListOddToEven[0]
        z2 = z2+configVarsZListOddToEven[1]
        z3 = z3+configVarsZListOddToEven[2]
        z4 = z4+configVarsZListOddToEven[3]
        z5 = z5+configVarsZListOddToEven[4]
        z6 = z6+configVarsZListOddToEven[5]


# Debug section: Print totals for right-upwards-then-downwards triplets
        if not ZDebugPrintOff:        
            print ' -----------'
            print ' '
            print "Totals for all triplets, after completing Row: ", topRow, "upwards-to-downwards"
            print "             (A-A-A) z1_total =", z1
            print "(A-A-B) plus (B-A-A) z2_total =", z2
            print "             (A-B-A) z3_total =", z3   
            print "             (B-A-B) z4_total =", z4
            print "(B-B-A) plus (A-B-B) z5_total =", z5 
            print "             (B-B-B) z6_total =", z6
            print ' '
            print ' -----------'
            print ' '
            print 'Closing a pass through for loop with i = ', i     
            print ' '    
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


def computeConfigVariables (arraySizeList, unitArray, h):
    
# Define all the configuration variables (x, y, w, and z) as elements of their respective lists,
#   and assign them their equilibrium values when all enthalpy parameters are set to 0.
#   Then, obtain the actual values for the configuration variables by calling procedures
#   specific to each configuration variable type (x, y, w, and z). 
#   Assign the returned configuration variables to list elements in the main
#   configVarsList, which is returned to teh calling procedure. 

    if not debugPrintOff:
        print "In computeConfigVariables"

#   Initialize the empty list for the full set of configuration variables
    configVarsList = list() # empty list
    
    configXVarsList = list() # empty list
    configXVarsList = computeConfigXVariables (arraySizeList, unitArray)
    x1 = configXVarsList[0]
    x2 = configXVarsList[1]    
    
    configYVarsList = list() # empty list

    configYVarsList = computeConfigYVariables (arraySizeList, unitArray)
    y1 = configYVarsList[0]
    y2 = configYVarsList[1]
    y3 = configYVarsList[2]     

    configWVarsList = list() # empty list    

    configWVarsList = computeConfigWVariables (arraySizeList, unitArray)
    w1 = configWVarsList[0]
    w2 = configWVarsList[1]
    w3 = configWVarsList[2]     

    configZVarsList = list() # empty list 
    configZVarsList = computeConfigZVariables (arraySizeList, unitArray)    
    z1 = configZVarsList[0]
    z2 = configZVarsList[1]
    z3 = configZVarsList[2]  
    z4 = configZVarsList[3]
    z5 = configZVarsList[4]
    z6 = configZVarsList[5]     

# Create the master Configuration Variables List; configVarsList, assign configuration variable
#   values, and return the list to the calling procedure 
               
    configVarsList = (x1, x2, y1, y2, y3, w1, w2, w3, z1, z2, z3, z4, z5, z6, unitArray)   
         
    if not debugPrintOff:
        print "In computeConfigVariables, about to return to Main"                
    return (configVarsList)    

####################################################################################################
####################################################################################################


def computeConfigVariablesDelta (h):
    
# Define all the configuration variables (x, y, w, and z) as elements of their respective lists,
#   and assign them their equilibrium values when all enthalpy parameters are set to 0.
#   Then, obtain the actual values for the configuration variables by calling procedures
#   specific to each configuration variable type (x, y, w, and z). 
#   Assign the returned configuration variables to list elements in the main
#   configVarsList, which is returned to teh calling procedure. 


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
    
# Create the master Configuration Variables List; configVarsList, assign configuration variable
#   values, and return the list to the calling procedure 
               
    configVarsList = (x1, x2, y1, y2, y3, w1, w2, w3, z1, z2, z3, z4, z5, z6)  
                
    return (configVarsList)    


    


####################################################################################################
####################################################################################################


def computeConfigVariablesAnalytic (h, configVarsList):
    
# Define all the configuration variables (x, y, w, and z) as elements of their respective lists,
#   and assign them their equilibrium values when all enthalpy parameters are set to 0.
#   Then, obtain the actual values for the configuration variables by calling procedures
#   specific to each configuration variable type (x, y, w, and z). 
#   Assign the returned configuration variables to list elements in the main
#   configVarsList, which is returned to teh calling procedure. 

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

# Create the master Configuration Variables List; configVarsList, assign configuration variable
#   values, and return the list to the calling procedure 
               
    configVarsListDelta = (x1Delta, x2Delta, y1Delta, y2Delta, y3Delta, 
    w1Delta, w2Delta, w3Delta, 
    z1Delta, z2Delta, z3Delta, z4Delta, z5Delta, z6Delta)   
             
    return (configVarsListDelta)    


####################################################################################################
####################################################################################################
#
# Procedure to print out the final set of configuration variables x'i, y'i, w'i, and z'i
#
####################################################################################################
####################################################################################################


def printConfigVarsComparison(h, configVarsList, configVarsListDelta):
  
            
    x1 = configVarsList[0]
    x2 = configVarsList[1]
    y1 = configVarsList[2]
    y2 = configVarsList[3]
    y3 = configVarsList[4] 
    w1 = configVarsList[5]
    w2 = configVarsList[6]
    w3 = configVarsList[7]        
    z1 = configVarsList[8]
    z2 = configVarsList[9]
    z3 = configVarsList[10] 
    z4 = configVarsList[11]
    z5 = configVarsList[12]
    z6 = configVarsList[13] 

    x1Delta = configVarsListDelta[0]
    x2Delta = configVarsListDelta[1]
    y1Delta = configVarsListDelta[2]
    y2Delta = configVarsListDelta[3]
    y3Delta = configVarsListDelta[4] 
    w1Delta = configVarsListDelta[5]
    w2Delta = configVarsListDelta[6]
    w3Delta = configVarsListDelta[7]        
    z1Delta = configVarsListDelta[8]
    z2Delta = configVarsListDelta[9]
    z3Delta = configVarsListDelta[10] 
    z4Delta = configVarsListDelta[11]
    z5Delta = configVarsListDelta[12]
    z6Delta = configVarsListDelta[13] 
        
    sumX = x1+x2 
    sumXDelta =  x1Delta+x2Delta 
    print ' '
    print '  For h =', h, 'The configuration variables are: '
    print ' '    
    print '                x1 = %.4f'  % x1,   '   x1Delta  = %.4f' % x1Delta 
    print '                x2 = %.4f'  % x2,   '   x2Delta  = %.4f' % x2Delta
    print '   Sum of the x(i) = %.4f'  % sumX, ' sumXDelta  = %.4f' % sumXDelta         
                        
    sumZ = z1+2.0*z2+z3 +z4+2.0*z5+z6
    sumZDelta = z1Delta+2.0*z2Delta+z3Delta +z4Delta+2.0*z5Delta+z6Delta    
    print ' '
    print "Totals for the z(i) variables:"
    print '        (A-A-A) z1 = %.4f'  % z1,   '   z1Delta  = %.4f' % z1Delta  
    print '(A-A-B & B-A-A) z2 = %.4f'  % z2,   '   z2Delta  = %.4f' % z2Delta   
    print '        (A-B-A) z3 = %.4f'  % z3,   '   z3Delta  = %.4f' % z3Delta  
    print '        (B-A-B) z4 = %.4f'  % z4,   '   z4Delta  = %.4f' % z4Delta 
    print '(A-B-B & B-B-A) z5 = %.4f'  % z5,   '   z5Delta  = %.4f' % z5Delta 
    print '        (B-B-B) z6 = %.4f'  % z6,   '   z6Delta  = %.4f' % z6Delta         
    print '   Sum of the z(i) = %.4f'  % sumZ, ' sumZDelta  = %.4f' % sumZDelta         
                        
             
                                     
    sumY = y1+2.0*y2+y3
    sumYDelta = y1Delta+2.0*y2Delta+y3Delta        
    print ' '
    print 'Totals for the y(i) variables:'
    print '          (A-A) y1 = %.4f'  % y1,   '   y1Delta  = %.4f' % y1Delta   
    print '    (A-B & B-A) y2 = %.4f'  % y2,   '   y2Delta  = %.4f' % y2Delta    
    print '          (B-B) y3 = %.4f'  % y3,   '   y3Delta  = %.4f' % y3Delta  
    print ' Multiplying by degeneracy factors:'
    print '   Sum of the y(i) = %.4f'  % sumY, ' sumYDelta  = %.4f' % sumYDelta 

    
    sumW = w1+2.0*w2+w3
    sumWDelta = w1Delta+2.0*w2Delta+w3Delta      
    print ' '
    print 'Totals for the w(i) variables:'
    print '        (A---A) w1 = %.4f'  % w1,   '   w1Delta  = %.4f' % w1Delta 
    print '(A---B & B---A) w2 = %.4f'  % w2,   '   w2Delta  = %.4f' % w2Delta   
    print '        (B---B) w3 = %.4f'  % w3,   '   w3Delta  = %.4f' % w3Delta 
    print ' Multiplying by degeneracy factors:'
    print '   Sum of the w(i) = %.4f'  % sumW, ' sumWDelta  = %.4f' % sumWDelta
    print         

   
                     


    return

####################################################################################################
#
# Function to increase the x1 value in the array
#
####################################################################################################
def adjustMatrixX1Up (arraySizeList, unitArray, configVarsList, h):
    
    totalUnits = float(arrayLength*arrayLayers)
    totalUnitsTimesTwo = totalUnits*2.0
                       
    x1 = float(configVarsList[0])/totalUnits
    x2 = float(configVarsList[1])/totalUnits 

    print ' '
    print ' In adjustMatrixXUp'    
# Randomly select a unit; if it is 1, change to 0
    unitRow = randrange(0, arrayLength)       
    unitCol = randrange(0, arrayLayers)        

    if unitArray[unitRow, unitCol] == 0: 
        unitArray[unitRow, unitCol] = 1
        print ' successfully changed unit in [', unitRow, ',', unitCol, '] to 1'
    else: print ' selected unit is already 1'  
            
    return unitArray


####################################################################################################
#
# Function to decrease the x1 value in the array
#
####################################################################################################

def adjustMatrixX1Down (arraySizeList, unitArray, configVarsList, h):
    
    totalUnits = float(arrayLength*arrayLayers)
    totalUnitsTimesTwo = totalUnits*2.0
                       
    x1 = float(configVarsList[0])/totalUnits
    x2 = float(configVarsList[1])/totalUnits 
    
    print ' '
    print ' In adjustMatrixXDown'   
        
# Randomly select a unit; if it is 1, change to 0
    unitRow = randrange(0, arrayLength)       
    unitCol = randrange(0, arrayLayers)        

    if unitArray[unitRow, unitCol] == 1: 
        unitArray[unitRow, unitCol] = 0
        print ' successfully changed unit in [', unitRow, ',', unitCol, '] to 0'
    else: print ' selected unit is already 0'        
                        
                                                                        
    return unitArray


####################################################################################################
#
# Function to adjust the array and bring x1 and x2 closer togetehr
#
####################################################################################################
    
def adjustMatrix (arraySizeList, unitArray, h, jrange, maxXDif):

    configVarsList = computeConfigVariables (arraySizeList, unitArray, h)
        
    totalUnits = float(arrayLength*arrayLayers)
    totalUnitsTimesTwo = totalUnits*2.0
                    

    for j in range (0, jrange, 1):
        x1 = float(configVarsList[0])/totalUnits
        x1DeltaPos = x1 - 0.5
        if x1DeltaPos > maxXDif: # x1 is too large
            print '  Improving the unit array for j = ', j, ' x1 is too large'
            unitArray = adjustMatrixX1Down (arraySizeList, unitArray, configVarsList, h)
            configVarsList = computeConfigVariables (arraySizeList, unitArray, h)
        x1DeltaNeg = 0.5 - x1
        if x1DeltaNeg > maxXDif: # x1 is too small    
            print '  Improving the unit array for j = ', j, ' x1 is too small'
            unitArray = adjustMatrixX1Up (arraySizeList, unitArray, configVarsList, h)
            configVarsList = computeConfigVariables (arraySizeList, unitArray, h)

    return unitArray               
                
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


def computeEntropyV2(configVarsList):
    print ' '
    print '  In computeEntropy '
    print ' '

    totalUnits = float(arrayLength*arrayLayers)
#    printConfigVars (configVarsList)      
                        
    x1 = float(configVarsList[0])/totalUnits
    x2 = float(configVarsList[1])/totalUnits

    sumX = x1+x2  
    print ' ' 
    print '                x1 = %.4f'  % x1
    print '                x2 = %.4f'  % x2
    print '   Sum of the x(i) = %.4f'  % sumX
    print ' ' 
    
    y1 = float(configVarsList[2])/(2*totalUnits)
    y2 = float(configVarsList[3])/(2*totalUnits)
    y3 = float(configVarsList[4])/(2*totalUnits) 
    w1 = float(configVarsList[5])/(2*totalUnits)
    w2 = float(configVarsList[6])/(2*totalUnits)
    w3 = float(configVarsList[7])/(2*totalUnits)        
    z1 = float(configVarsList[8])/(2*totalUnits)
    z2 = float(configVarsList[9])/(2*totalUnits)
    z3 = float(configVarsList[10])/(2*totalUnits) 
    z4 = float(configVarsList[11])/(2*totalUnits)
    z5 = float(configVarsList[12])/(2*totalUnits)
    z6 = float(configVarsList[13])/(2*totalUnits) 

    sumY = y1+y2+y3    
    print
    print "Totals for the y(i) variables:"
    print '          (A-A) y1 = %.4f'  % y1 
    print '    (A-B & B-A) y2 = %.4f'  % y2   
    print '          (B-B) y3 = %.4f'  % y3 
    print '   Sum of the y(i) = %.4f'  % sumY
    print 
    
    sumW = w1+w2+w3 
    print
    print "Totals for the w(i) variables:"
    print '        (A---A) w1 = %.4f'  % w1 
    print '(A---B & B---A) w2 = %.4f'  % w2   
    print '        (B---B) w3 = %.4f'  % w3 
    print '   Sum of the w(i) = %.4f'  % sumW
    print         

    sumZ = z1+z2+z3 +z4+z5+z6
    print
    print "Totals for the z(i) variables:"
    print '        (A-A-A) z1 = %.4f'  % z1 
    print '(A-A-B & B-A-A) z2 = %.4f'  % z2   
    print '        (A-B-A) z3 = %.4f'  % z3 
    print '        (B-A-B) z4 = %.4f'  % z4 
    print '(A-B-B & B-B-A) z5 = %.4f'  % z5 
    print '        (B-B-B) z6 = %.4f'  % z6         
    print '   Sum of the z(i) = %.4f'  % sumZ
    print                    
                                                
    Lfx = LfFunc(x1)+LfFunc(x2)
    print ' lnFunc (x) = %.4f'  % (Lfx)

    Lfy = LfFunc(y1) + 2*LfFunc(y2) + LfFunc(y3)
    print ' lnFunc (y) = %.4f'  % (Lfy)  
      
    Lfw = LfFunc(w1) + 2*LfFunc(w2) + LfFunc(w3)
    print ' lnFunc (w) = %.4f'  % (Lfw)

    Lfz = LfFunc(z1) + 2*LfFunc(z2) + LfFunc(z3) + LfFunc(z4) + 2*LfFunc(z5) + LfFunc(z6)
    print ' lnFunc (z) = %.4f'  % (Lfz)
        
    s = -(2*Lfy+Lfw-Lfx-2*Lfz)

    
    print ' '
    print ' The entropy is: %.4f' % (s)  
  
   
     
    return (s)

    
####################################################################################################
####################################################################################################
#
# Function to compute the entropy, enthalpy, and free energy for a smooth calculation through 
#    a stepped range of h values
#
####################################################################################################
####################################################################################################


def computeThermodynamicVars(h, configVarsList):

    totalUnits = float(arrayLength*arrayLayers)
    totalUnitsTimesTwo = totalUnits*2.0
                       
    x1 = float(configVarsList[0])/totalUnits
    x2 = float(configVarsList[1])/totalUnits 
    
    y1 = float(configVarsList[2])/totalUnitsTimesTwo
    y2 = float(configVarsList[3])/totalUnitsTimesTwo
    y3 = float(configVarsList[4])/totalUnitsTimesTwo 
    
    sumY = y1 + y2 + y3

    w1 = float(configVarsList[5])/totalUnitsTimesTwo
    w2 = float(configVarsList[6])/totalUnitsTimesTwo
    w3 = float(configVarsList[7])/totalUnitsTimesTwo        

    sumW = w1 + w2 + w3

    z1 = float(configVarsList[8])/totalUnitsTimesTwo
    z2 = float(configVarsList[9])/totalUnitsTimesTwo
    z3 = float(configVarsList[10])/totalUnitsTimesTwo 
    z4 = float(configVarsList[11])/totalUnitsTimesTwo
    z5 = float(configVarsList[12])/totalUnitsTimesTwo
    z6 = float(configVarsList[13])/totalUnitsTimesTwo 
              

    sumZ = z1 + z2 + z3 + z4 + z5 + z6

    print ' '
    print ' x1 = %.4f'  % (x1), ' x2 = %.4f'  % (x2)
    print ' y1 = %.4f'  % (y1), ' y2 = %.4f'  % (y2), ' y3 = %.4f'  % (y3) #, ' sumY = %.4f'  % (sumY)                                                 
    print ' w1 = %.4f'  % (w1), ' w2 = %.4f'  % (w2), ' w3 = %.4f'  % (w3) #, ' sumW = %.4f'  % (sumW) 

    print ' '
    print ' z1 = %.4f'  % (z1), ' z2 = %.4f'  % (z2), ' z3 = %.4f'  % (z3)
    print ' z6 = %.4f'  % (z6), ' z5 = %.4f'  % (z5), ' z4 = %.4f'  % (z4)  #, 'sumZ = %.4f'  % (sumW)     
        
                                                                                                                                                    
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

    print ' '
    print ' The computed thermodynamic quantities:'            
    print '    h = %.4f' % (h), '    epsilon1 = %.4f' % (epsilon1) 
    print '    negative entropy s = %.4f' % (negS)
    print '    enthalpy = %.4f' % (enthalpy1)      
    print '    free energy = %.4f' % (freeEnergy)  
                                  
    return (sysValsList)   
         
    

        
####################################################################################################
####################################################################################################
#
# Code Documentation - top-down: 
# The MAIN module comprising of calls to:
#  (1) Welcome
#  (2) Obtain array size specifications for an MxN array of 0 or 1 units (currently pre-defined patterns)
#  (3) computeConfigVariables: Compute the configuration variables for the entire grid
#    -- NOTE: This is the major work-horse function in this program; see documentation below
#   
#  Documentation for computeConfigVariables and its supporting functions:  
#  (3) computeConfigVariables: Computes the configuration variables for the grid; x, y, w, and z.
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

    welcome()

    global debugPrintOff
    global detailedDebugPrintOff
    global ZDebugPrintOff
    global blnkspc
    global arrayLength
    global arrayLayers
    global evenLayers
    global pairs
    evenLayers = True
            
    arraySizeList = list() # empty list
    arraySizeList = obtainArraySizeSpecs ()
    arrayLength = arraySizeList[0]
    arrayLayers = arraySizeList [1]
    sysVarsList = list() 
    maxXDif = 0.015    # maximum difference between x1 from randomly-generated array and the
                      # equiprobable distribution value of x1 = 0.5  
    jrange = 3        # maximal number of steps allowed to improve the x1 distribution

    blnkspc=' '

    if arrayLayers % 2 == 0: evenLayers == True #then an even number of layers

    # Determine the total number of PAIRS of zigzag chains
    pairsLayers = arrayLayers/2
    pairs = int(pairsLayers) 

    debugPrintOff = True
    detailedDebugPrintOff = True
    ZDebugPrintOff = True
    
    if not debugPrintOff:
        print ' '
        print 'Debug printing is on'  #debugPrintOff false
    else: 
        print 'Debug printing is off' #debugPrintOff true    
            





####################################################################################################
# Compute configuration variables
####################################################################################################                

#    if not debugPrintOff:
#        print ' '   
#        print "In __main__, about to call computeConfigVariables"
                                   

                    
#    debugPrintOff = False
#    if not debugPrintOff:
#        print ' '
#        print 'Debug printing is on'  #debugPrintOff false
#    else: 
#        print 'Debug printing is off' #debugPrintOff true
      
                        
    eps0 = 0.0  # The enthalpy for single unit activation is set to zero for this code

# Determine if we are probabilistically generating a pattern (patternProb = 0)
#   or if we are selecting a pre-stored pattern (patternProb = 1 ... N)
    patternProb = 0
# Pick a starting value for h; it should be less than 1; it should be in the realm of 0.7 - 0.8    
    h0 = 0.73
# Pick the number of steps (increasing 0.01) for increasing h
    hrange = 3
# Pick the increment for increasing h
    incr = 0.01 
       
             
 
    xArray = np.zeros(hrange, dtype=np.float)
    negSArray = np.zeros(hrange, dtype=np.float)
    fEps0Array = np.zeros(hrange, dtype=np.float)         
    fEps1Array = np.zeros(hrange, dtype=np.float)  
    fEnergyArray = np.zeros(hrange, dtype=np.float) 
         
    print 
    print 'Initializing a matrix of M x L units, where:'
    print '  M (the arraylength is):', arrayLength 
    print '  L (the layers is):', arrayLayers   
    
    h = h0
    for i in range (0, hrange, 1):        
        h = h + incr
        configVarsList = list ()  # empty list
        unitArray    = initializeMatrix (arraySizeList, patternProb, h)        
# debug print: 
#        print unitArray


# adjust the unit array so that x1 approximately = x2
        unitArray = adjustMatrix (arraySizeList, unitArray, h, jrange, maxXDif)

# obtain the configuration variables (this step was also done while adjusting the array)       
        configVarsList = computeConfigVariables (arraySizeList, unitArray, h)
                                                                
# obtain the thermodynamic variables 
        sysVarsList = computeThermodynamicVars (h, configVarsList)   

# store the thermodynamic variables to plot later
        xArray[i]=h
        negSArray[i]    = sysVarsList[0] 
        fEps0Array[i]   = sysVarsList[1]
        fEps1Array[i]   = sysVarsList[2]-0.6  
        fEnergyArray[i] = sysVarsList[3]-0.6                                                
                                                                                                                              

    print ' '
    print ' The thermodynamic quantities plot vs. h' 
    print '   The negEntropy is in blue, '
    print '   The per-unit enthalpy is zero, and is not shown, '
    print '   The interaction enthalpy (eps1*y2) is in maroon, and is shifted by 0.6,' 
    print '     and the free energy is in red, also shifted by 0.6.'   
    pylab.figure(1)
    pylab.plot (xArray,negSArray)    
    pylab.plot (xArray,fEps1Array,'m')
    pylab.plot (xArray,fEnergyArray,'r')
    pylab.show()  
 
    
          
####################################################################################################
# Conclude specification of the MAIN procedure
####################################################################################################                
    
if __name__ == "__main__": main()

####################################################################################################
# End program
####################################################################################################                
                                                                                                                                                                                                                                                                    