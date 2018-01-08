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
# Jan. 7, 2018: Beginning work on perturbations analysis: Goal is to show that in the realm of 
#  h > 1 (like-with-like nearest neighbor pairings encouraged), it will be possible to return to 
#  a prior state even after perturbation, so that there is greater likelihood of pattern 
#  persistence as compared with systems where the interaction energy is 0 (h = 1). 
#
#  Steps are: 
#  1) Create unitArray as with previous versions of code; this will be for a single set of 
#     x1 and h values; the unitArray will adjust to achieve free energy minimum,
#  2) Perturb the unit array by randomly changing a prescribed amount of nodes,
#  3) Allow the unitArray to come to a free energy minimum again; compare actual node (on/off) 
#     patterns vs. original
# Jan. 4, 2018: Found/fixed bug in communicating updated unitArray to the next calling function
#
# Dec. 17, 2017: This program generates a random population of nodes for a 2-D CVM according to a 
#   pre-specified ratio of x1 to x2. 
#   To adjust the ratio, tweak the parameter x-ratio in main. 
#
#
# Dec. 17, 2017: This particular version of the program is simpler than the documentation below 
#   suggests. It DOES NOT use the h-values, except as a means for stepping through a loop. 
#   All that it does is randomly populate a 2-D CVM grid, with a random generation of  
#   "A" and "B" (x1 and x2) nodes, in a ratio that can be specified inside **main**. 
#   Fairly often, this random generation does not yield the precise ratio of x1 & x2, so it does an 
#   "adjustment step." 
#   During this adjustment (currently set to a max jrange, specified in **main**), it will look for 
#   a node to flip to alter the x1 / x2 ratio. If it succeeds in finding a "flippable" node (one that 
#   brings the x1 and x2 proportions closer to desired ratio), then it keeps the flip. 
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
# At equilibrium, when eps1 = 0, z1 = z6 = z3 = z4 = 0.125; z2 = z5 = 0.125, however,
#  they are "degenerate" so their total values each = 0.25. 
#
# I am using a revised enthalpy equation (cf the original 2016 Brain Sciences paper); 
#  the new equation uses: 2y2-y1-y3. 
#
# The following equations are taken from the 2-D CVM analtyic solution. 
#   They are valid ONLY for the equiprobable case; x1 = x2. 
#   These equations are put in here for reference; they are not used in the current version of code. 
#
# The equations are taken from my 2014 Technical Report on the 2-D CVM. 
# THEY NEED REVISION to account for the different enthalpy equation that I'm using now. 
#    hSquared = h*h
#    hFourth  = hSquared*hSquared
#    denom    = 8.0*(hFourth - 6.0*hSquared + 1.0)  
#    z3Analytic = (hSquared - 3.0)*(hSquared + 1.0)/denom      #  App. B, Eqn. 29
#    z1Analytic = (1.0 - 3.0*hSquared)*(hSquared + 1.0)/denom  #  App. B, Eqn. 30
#    y1Analytic = z1Analytic + 0.5*(0.5 - z1Analytic - z3Analytic)
#    y2Analytic = z3Analytic + 0.5*(0.5 - z1Analytic - z3Analytic)
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
    print 'Version 1.2, 01/07/2018, A.J. Maren'
    print 'This version computes the behavior of a perturbed unitArray,'
    print  ' based on minimizing the free energy both before and after perturbation.'
    print ' ' 
    print 'By changing parameters in the main code, the user can select:'
    print '   x1 value - relative proportion of A (on) vs. B (off) nodes'
    print '   h value - controls the interaction energy, thus affects the '
    print '     y, w, and z configuration values for a given x, and'
    print '   perturbPrcnt - the percentage of nodes in the unitArray that'
    print '     will be flipped from one state to another.'
    print ' ' 
    print 'For comments, questions, or bug-fixes, contact: alianna.maren@northwestern.edu'
    print 'Alternate email address: alianna@aliannajmaren.com'
    print ' '
    print '  NOTE: In these calculations, x1 = A (units are at value 1),'
    print '                           and x2 = B (units are at value 0).'
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


####################################################################################################
####################################################################################################
#
# Function to randomly-generate an array, and then permute it to achieve the desired 
#    z1 & z3 values.
#
#    Inputs:    arraySizeList: a list of two integers; arrayLength and layers
#               h: the interaction enthalpy parameter
#    Return: the matrix unit_array, a matrix of 0's and 1's.
#
####################################################################################################
####################################################################################################

def initializeGeneratedMatrix (arraySizeList, h, x1TargetVal):
    
# The current version of the program has the h-value passed into this routine.
# It is not used, except to compute the "analytic" values for the configuration variables - 
#   which are relevant only in the equiprobable case. 
# Later versions of this program will accept h in order to adjust the configuration variables by 
#   changing the grid population.   

    localArrayLength = arraySizeList[0]
    localArrayLayers = arraySizeList[1]

# Create the matrix 'unit_array' so that it has a random population of 0's and 1's.
#    unit_array = np.random.choice([0, 1],size=(localArrayLayers,localArrayLength)) # Create an array filled with random values
# Note: this function can be used to create proportional distributions: np.random.choice([0, 1], size=(10,), p=[1./3, 2./3])

    x2TargetVal = 1.-x1TargetVal
    unitArray = np.random.choice([0, 1],size=(localArrayLayers,localArrayLength), p=[x2TargetVal, x1TargetVal])
    
    return unitArray


####################################################################################################
####################################################################################################
#
# Procedure to initialize the matrix with a randomly-generate an array, and then permute it to 
#   achieve the desired z1 & z3 values. This is NOT part of the "permutation" experiments;
#   it is getting the first unitArray to be one which is at equilibrium for a given set of 
#   x and h-values. 
#
#    Inputs:    arraySizeList: a list of two integers; arrayLength and layers
#               h: the interaction enthalpy parameter
#               x1TargetVal: the desired x1 in the unitArray; the first rough pass
#                 is randomly generated, then refined until the fraction x1 is 
#                 within a tolerance (specified in **main** as maxXDif) of the desired x1.
#               maxXDif: the maximum difference allowed for the x1 from the randomly-generated
#                 array and the desired x1.
#    Return: the matrix unit_array, a matrix of 0's and 1's.
#
####################################################################################################
####################################################################################################

 
def initializeMatrix (arraySizeList, h, x1TargetVal, maxXDif):
         
    localArrayLength = arraySizeList[0]
    localArrayLayers = arraySizeList[1]


# Note: The passed value patternProb is used to determine if we are returning a stored pattern, or
#       are probabilistically-generating our data.
#       If patternProb = 0: probabilistic generation, dependent on h 
#       If patternProb > 1: select one of the N stored patterns (1 ... N)
        

# This is when we will randomly generate a matrix, with specified x-proportions    

    unitArray = initializeGeneratedMatrix (arraySizeList,h, x1TargetVal)


    if not debugPrintOff:
        print ' '    
        print ' The maximum difference in x from the target values is: ', maxXDif
        
    #    if not beforeAndAfterAdjustedMatrixPrintOff:
    #        print
    #        print 'The distribution among states A and B (x1 and x2) units is:'    
    #        print "  ( A ) x1_total =", x1_total, "( B ) x2_total =", x2_total  
    #        print '  For a total of ', x1_total + x2_total, ' units.'
    #        print ' '         
          
    # Print the unit_array so that it shows as a zigzag chain (or in the 2-D case, layers of zigzag chains).
        print 
        print 'The L x M array of units, where M (across) =', localArrayLength, 'and L (layers) =', localArrayLayers
        print 

# Determining "pairs" - the total number of pairs of zigzag chains - is done in __main__; "pairs" is a global variable


    if not debugPrintOff:
        print ' '        
        for i in range (0, pairs):
            actualEvenRowNum = 2*i
            print 'Row', actualEvenRowNum, ':', blnkspc, 
            for j in range(0,localArrayLength):
                print unitArray[actualEvenRowNum,j], blnkspc,
            print 
            actualOddRowNum = 2*i+1
            print 'Row ', actualOddRowNum, ':', blnkspc,
            print (blnkspc),
            for j in range(0,localArrayLength):
                print unitArray[actualOddRowNum,j], blnkspc,
            print 
        print ' '
    

    #    print ' ' 
    #     print ' This is an effort to build a -, X representation of the unit array'
    #     easyGridArray = createEasyToReadGrid (arraySizeList, unitArray)
    #     print easyGridArray 
    #     print ' ' 
    
    return unitArray



####################################################################################################
####################################################################################################
#
# Procedure to compute configuration variables x'i and return as elements of list configXVarsList
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
    if not detailedAdjustMatrixPrintOff:
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
    # END (function to compute the Y configuration variables for an EVEN row; this 
    #  includes the wrap-around)


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
    # END (function to compute the Y configuration variables for an ODD row; this 
    #  includes the wrap-around)


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
    # END (high-level function calling the EVEN and ODD functions to compute the Y
    #  configuration variables)




####################################################################################################
####################################################################################################
#
# Procedure to compute horizontal configuration variables w'i and return as elements of list configWVarsList
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
    # END (high-level function calling the HORIZONTAL and VERTICAL functions to compute the W
    #  configuration variables)




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


def computeConfigVariables (arraySizeList, unitArray):
    
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
    X1 = configXVarsList[0]
    X2 = configXVarsList[1]    
    
    configYVarsList = list() # empty list

    configYVarsList = computeConfigYVariables (arraySizeList, unitArray)
    Y1 = configYVarsList[0]
    Y2 = configYVarsList[1]
    Y3 = configYVarsList[2]     

    configWVarsList = list() # empty list    

    configWVarsList = computeConfigWVariables (arraySizeList, unitArray)
    W1 = configWVarsList[0]
    W2 = configWVarsList[1]
    W3 = configWVarsList[2]     

    configZVarsList = list() # empty list 
    configZVarsList = computeConfigZVariables (arraySizeList, unitArray)    
    Z1 = configZVarsList[0]
    Z2 = configZVarsList[1]
    Z3 = configZVarsList[2]  
    Z4 = configZVarsList[3]
    Z5 = configZVarsList[4]
    Z6 = configZVarsList[5]     

# Create the master Configuration Variables List; configVarsList, assign configuration variable
#   values, and return the list to the calling procedure 
               
    configVarsList = (X1, X2, Y1, Y2, Y3, W1, W2, W3, Z1, Z2, Z3, Z4, Z5, Z6, unitArray)   
         
    if not debugPrintOff:
        print "In computeConfigVariables, about to return to **main**"                
    return (configVarsList)    



####################################################################################################
####################################################################################################


def computeConfigVariablesDelta (h):
    
# This is an analytic function; not used for the random generation and adjustment steps in the current
#   code; it is used in other code instances

# Define all the configuration variables (x, y, w, and z) as elements of their respective lists,
#   and assign them their equilibrium values when all enthalpy parameters are set to 0.
#   Then, obtain the actual values for the configuration variables by calling procedures
#   specific to each configuration variable type (x, y, w, and z). 
#   Assign the returned configuration variables to list elements in the main
#   configVarsList, which is returned to teh calling procedure. 

# NOTE: These equations are based on an enthalpy of just 2y2;
#   they need to be replaced for with newly-derived equations where
#   the enthalpy is 2y2 - y1 - y3. 

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
    
# This is exploratory code - written prior to generating actual unitArrays and assessing their
#  configuration variables, etc;
# Its purpose was to give a mathematically-plausible set of configuration variables 
#  for which the thermodynamic variables could then be computed;
# There is no knowledge at all if these values are even near equilibrium ... 
# Actually modifying the unitArray to achieve an equilibrium configuration came much later ... 


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
#
# Function to increase the x1 value in the array
#
####################################################################################################

def adjustMatrixX1Up (arraySizeList, unitArray, configVarsList, h):
    
    totalUnits = float(arrayLength*arrayLayers)
    totalUnitsTimesTwo = totalUnits*2.0
                       
    x1 = float(configVarsList[0])/totalUnits
    x2 = float(configVarsList[1])/totalUnits 

    if not detailedAdjustMatrixPrintOff:
        print ' '
        print ' In adjustMatrixXUp'    
# Randomly select a unit; if it is 1, change to 0
    unitRow = randrange(0, arrayLength)       
    unitCol = randrange(0, arrayLayers)        

    if unitArray[unitRow, unitCol] == 0: 
        unitArray[unitRow, unitCol] = 1
        if not detailedAdjustMatrixPrintOff:
            print ' successfully changed unit in [', unitRow, ',', unitCol, '] to 1'
    else: 
        if not detailedAdjustMatrixPrintOff:
            print ' selected unit is already 1'  
            
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

    if not detailedAdjustMatrixPrintOff:    
        print ' '
        print ' In adjustMatrixXDown'   
        
# Randomly select a unit; if it is 1, change to 0
    unitRow = randrange(0, arrayLength)       
    unitCol = randrange(0, arrayLayers)        

    if unitArray[unitRow, unitCol] == 1: 
        unitArray[unitRow, unitCol] = 0
        if not detailedAdjustMatrixPrintOff:
            print ' successfully changed unit in [', unitRow, ',', unitCol, '] to 0'
    else: 
        if not detailedAdjustMatrixPrintOff:
            print ' selected unit is already 0'        
                        
                                                                        
    return unitArray


####################################################################################################
#
# Function to adjust the array and bring randomly-generated x1 and the target x1 closer together
#
####################################################################################################
    
def adjustMatrix (arraySizeList, unitArray, h, jrange, maxXDif, x1TargetVal, beforeAndAfterAdjustedMatrixPrintOff):
     
    configVarsList = computeConfigVariables (arraySizeList, unitArray)
        
    totalUnits = float(arrayLength*arrayLayers)
    totalUnitsTimesTwo = totalUnits*2.0 

    # Test to see if matrix adjustment is needed
    x1 = float(configVarsList[0])/totalUnits
    x1DeltaPos = x1 - x1TargetVal
    x1DeltaNeg = x1TargetVal - x1 
    

    if not debugPrintOff:
        print ' '
        print '    The initial value for x1 is %.4f' % (x1)
                                              
    if x1DeltaPos < maxXDif: # x1 is not too large
        x1DeltaPosAccept = True
    else:
        x1DeltaPosAccept = False        
        
    if x1DeltaNeg < maxXDif: # x1 is not too small
        x1DeltaNegAccept = True
    else:
        x1DeltaNegAccept = False                 

    x1Accept = False 
    if x1DeltaPosAccept:
        if x1DeltaNegAccept:
            x1Accept = True                
       

    if not debugPrintOff:
        print ' '                        
        if x1Accept:
                print '  '
                print '  Matrix adjustment is not needed'
                print '    The initial value for x1 is %.4f' % (x1) 
                print '    This is less than the maximum difference of ', maxXDif 
                print ' '
        if x1DeltaPosAccept:
                print '    The difference between x1 and the x1 target value is %.4f'% (x1DeltaPos) 
        else:
            if x1DeltaNegAccept:
                print '    The difference between x1 and the x1 target value is %.4f'% (x1DeltaNeg)              
                                                            

    if not x1DeltaPosAccept: # x1 IS too large; begin adjusting matrix                                                                                                                      
        for j in range (0, jrange, 1):
            x1 = float(configVarsList[0])/totalUnits
            x1DeltaPos = x1 - x1TargetVal            
            if x1DeltaPos < maxXDif:  # x1 is not too large; exit FOR loop
                if not detailedAdjustMatrixPrintOff: 
                    print '  '
                    print '  Matrix adjustment completed in ', j+1, 'steps'
                    print ' '
                break              
            if not detailedAdjustMatrixPrintOff:
                print '  Improving the unit array for j = ', j, ' x1 is too large'
            unitArray = adjustMatrixX1Down (arraySizeList, unitArray, configVarsList, h)
            configVarsList = computeConfigVariables (arraySizeList, unitArray)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
    if not x1DeltaNegAccept: # x1 IS too small; begin adjusting matrix                                                                                                                      
        for j in range (0, jrange, 1):
            x1 = float(configVarsList[0])/totalUnits
            x1DeltaNeg = x1TargetVal - x1 
            if x1DeltaNeg < maxXDif:  # x1 is not too large; exit FOR loop
                if not detailedAdjustMatrixPrintOff: 
                    print '  '
                    print '  Matrix adjustment completed in ', j+1, 'steps'
                    print ' '
                break              
            if not detailedAdjustMatrixPrintOff:
                print '  Improving the unit array for j = ', j, ' x1 is too large'
            unitArray = adjustMatrixX1Up (arraySizeList, unitArray, configVarsList, h)
            configVarsList = computeConfigVariables (arraySizeList, unitArray)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
                                
    return unitArray
    # END adjustMatrix  
                 

####################################################################################################
#
# Function to find an x1 node, which will be a candidate for flipping to lower the total FE value
#
####################################################################################################

def findCandidateX1node (arraySizeList, unitArray, maxRange, findFEMinimumValsDetailsBoolOff):
    
    arrayLength = arraySizeList[0]
    arrayLayers = arraySizeList[1]
    
    candidateX1RowColList = (0,0,0)    
                       
    successBool = 0
    for k in range (0, maxRange, 1):
        unitRow = randrange(0, arrayLength)       
        unitCol = randrange(0, arrayLayers)

        if unitArray[unitRow,unitCol] == 1:
            candidateX1RowColList = [unitRow, unitCol,k]
            numCandidates = k
            successBool = 1
            if not findFEMinimumValsDetailsBoolOff:
                print ' '
                print '  For k = ', k, 'Candidate x1 found at Row: ', unitRow, 'Col: ', unitCol
                numCandidates = k
                successBool = 1
            break
        # if not findFEMinimumValsDetailsBoolOff:
    
    if successBool == 0:           
        if not findFEMinimumValsDetailsBoolOff:
            print ' ' 
            print ' In findCandidateX2node' 
            print ' Did not find a suitable x2 candidate after ', maxRange, ' attempts'
            print ' '  
                                    
    return candidateX1RowColList
    # END findCandidateX1 node


####################################################################################################
#
# Function to find an x2 node, which will be a candidate for flipping to lower the total FE value
#
####################################################################################################

def findCandidateX2node (arraySizeList, unitArray, maxRange, findFEMinimumValsDetailsBoolOff):
    
    arrayLength = arraySizeList[0]
    arrayLayers = arraySizeList[1]
    
    candidateX2RowColList = (0,0,0)    

    successBool = 0                                              
    for k in range (0, maxRange, 1):
        unitRow = randrange(0, arrayLength)       
        unitCol = randrange(0, arrayLayers)

        if unitArray[unitRow,unitCol] == 0:
            candidateX2RowColList = [unitRow, unitCol,k]
            successBool = 1
            numCandidates = k            
            if not findFEMinimumValsDetailsBoolOff:
                print ' '
                print '  For k = ', k, 'Candidate x2 found at Row: ', unitRow, 'Col: ', unitCol
            break           
        # if not findFEMinimumValsDetailsBoolOff:            
        
    if successBool == 0: 
        if not findFEMinimumValsDetailsBoolOff:
            print ' ' 
            print ' In findCandidateX2node' 
            print ' Did not find a suitable x2 candidate after ', maxRange, ' attempts'
            print ' '  
                                                                                        
    return candidateX2RowColList
    # END findCandidateX2 node



####################################################################################################
#
# Function to fill a new unit matrix so it is identical to the starting matrix
#
####################################################################################################
    
def createIdenticalUnitArray (arraySizeList, unitArray):
    
    localArrayLength = arraySizeList[0]
    localArrayLayers = arraySizeList[1]
    newUnitArray = np.random.choice([0, 1],size=(localArrayLayers,localArrayLength))

#  How to copy a list: b = [x for x in a] 
#  Another way: b = list(a)   
         
# We have an  L x M array of units, where M (across) =', localArrayLength, 'and L (layers) =', localArrayLayers
             
    for l in range (0, localArrayLayers, 1):
        for m in range (0, localArrayLayers, 1): 
            newUnitArray[m,l] =  unitArray[m,l]              
                        
    return newUnitArray





####################################################################################################
#
# Function to adjust the array (while keeping x1, x2 const) to change the other config variables
#   and bring the FE value to a minimum
#
####################################################################################################
    
def adjustMatrixFEMinimum (arraySizeList, unitArray, h, maxRange):

# Note:     sysValsList = (negS, enthalpy0, enthalpy1, freeEnergy): returned from computeThermValues          
# Note:     configVarsList = (x1, x2, y1, y2, y3, w1, w2, w3, z1, z2, z3, z4, z5, z6, unitArray); returned from computeConfigVars     

    maxRange = 30
    
    findFEMinimumValsBoolOff = True
    findFEMinimumValsDetailsBoolOff = True 
    
    totalTrials = 200
    step = 1 
    
    x1ValsArray   = np.zeros(totalTrials, dtype=np.float)
    y1ValsArray   = np.zeros(totalTrials, dtype=np.float)
    y2ValsArray   = np.zeros(totalTrials, dtype=np.float)
    y3ValsArray   = np.zeros(totalTrials, dtype=np.float)
    w1ValsArray   = np.zeros(totalTrials, dtype=np.float)
    w2ValsArray   = np.zeros(totalTrials, dtype=np.float)
    w3ValsArray   = np.zeros(totalTrials, dtype=np.float)
    z1ValsArray   = np.zeros(totalTrials, dtype=np.float)
    z3ValsArray   = np.zeros(totalTrials, dtype=np.float)       
    negSValsArray = np.zeros(totalTrials, dtype=np.float)
    enthalpy0Array= np.zeros(totalTrials, dtype=np.float)
    enthalpy1Array= np.zeros(totalTrials, dtype=np.float)
    freeEnergyArray=np.zeros(totalTrials, dtype=np.float)  
    
    newUnitArray = createIdenticalUnitArray (arraySizeList, unitArray)            

    # Obtain fractional values for the configuration variables (previously retrieved total counts)        
    totalUnits = float(arrayLength*arrayLayers)
    totalUnitsTimesTwo = totalUnits*2.0   

    # Obtain the starting values for the unitArray configuration variables and associated thermodynamic variables        
    configVarsListOrig = computeConfigVariables (arraySizeList, unitArray)
    x1Orig = float(configVarsListOrig[0])/totalUnits
    y1Orig = float(configVarsListOrig[2])/totalUnitsTimesTwo 
    y2Orig = float(configVarsListOrig[3])/(totalUnitsTimesTwo*2.0) 
    y3Orig = float(configVarsListOrig[4])/totalUnitsTimesTwo        
    z1Orig = float(configVarsListOrig[8])/totalUnitsTimesTwo      
    z3Orig = float(configVarsListOrig[10])/totalUnitsTimesTwo
    sysValsListOrig = computeThermodynamicVars(h, configVarsListOrig)
    negEntropyOrig  = sysValsListOrig[0]
    enthalpy1Orig   = sysValsListOrig[2]  
    FEValueOrig     = sysValsListOrig[3]                  
                                
    successfulFlips = 0 # Count the total number of successful flips during the totalTrials

    if not findFEMinimumValsBoolOff:
        print ' ' 
        print ' In adjustMatrixFEMinimum with h = ', h
        print ' Starting to adjust config vars to minimize free energy'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
        print '   Trial parameters: x1(starting) =  %.3f' % (x1Orig), ' h =  %.2f' % (h)
        print ' '
    
            
    for i in range (0, totalTrials, step): 
        x1CandidateRowColList = findCandidateX1node (arraySizeList, unitArray, maxRange, findFEMinimumValsDetailsBoolOff)
        x1Row   = x1CandidateRowColList[0]
        x1Col   = x1CandidateRowColList[1]
        x1Flips = x1CandidateRowColList[2]

        x2CandidateRowColList = findCandidateX2node (arraySizeList, unitArray, maxRange, findFEMinimumValsDetailsBoolOff)
        x2Row   = x2CandidateRowColList[0]
        x2Col   = x2CandidateRowColList[1]
        x2Flips = x2CandidateRowColList[2]  
             
        if not findFEMinimumValsDetailsBoolOff:
            print ' ' 
            print '  For i = ', i, ' the candidate x1 value is at Row: ', x1Row, ' Column: ', x1Col, ' after ', x1Flips, 'flips.'
            print '  For i = ', i, ' the candidate x2 value is at Row: ', x2Row, ' Column: ', x2Col, ' after ', x2Flips, 'flips.'

        # Swap the unit values in newUnitArray: keeping x1 the same, but changing the other config variable values
        newUnitArray[x1Row,x1Col] = 0
        newUnitArray[x2Row,x2Col] = 1
                                                                                                      
        # Obtain the configuration variable values from EACH of the old and new unitArrays
        configVarsListOld = computeConfigVariables (arraySizeList, unitArray)
        configVarsListNew = computeConfigVariables (arraySizeList, newUnitArray)
    
        # Obtain the thermodynamic variables corresponding to each of the old and new unitArrrays
        sysValsListOld = computeThermodynamicVars(h, configVarsListOld)
        negEntropyOld  = sysValsListOld[0]
        enthalpy1Old   = sysValsListOld[2]  
        FEValueOld     = sysValsListOld[3]    
   
        sysValsListNew= computeThermodynamicVars(h, configVarsListNew) 
        negEntropyNew  = sysValsListNew[0]
        enthalpy1New   = sysValsListNew[2] 
        FEValueNew     = sysValsListNew[3]      
                
        x1Old = float(configVarsListOld[0])/totalUnits
        x1New = float(configVarsListNew[0])/totalUnits
        y2Old = float(configVarsListOld[3])/(totalUnitsTimesTwo*2.)
        y2New = float(configVarsListNew[3])/(totalUnitsTimesTwo*2.)
        z1Old = float(configVarsListOld[8])/totalUnitsTimesTwo
        z1New = float(configVarsListNew[8])/totalUnitsTimesTwo         
        z3Old = float(configVarsListOld[10])/totalUnitsTimesTwo
        z3New = float(configVarsListNew[10])/totalUnitsTimesTwo

        x1ValsArray[i]   = x1New
        y2ValsArray[i]   = y2New 
        z1ValsArray[i]   = z1New
        z3ValsArray[i]   = z3New      
        negSValsArray[i]   = negEntropyNew
        enthalpy1Array[i]   = enthalpy1New
        freeEnergyArray[i]   = FEValueNew 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
        # Print the results of the single flip: Note thermodynamic variables
        if not findFEMinimumValsDetailsBoolOff:                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
            print ' '
            print ' ' 
            print ' The results of a single flip in two unitArray positions, at trial number', i, ' of ', totalTrials, ' total trials' 
            print ' '        
            print '          x1      y2       z1       z3       negS    enthalpy1   FE '
            print ' Old:   %.4f' % (x1Old), '  %.4f' % (y2Old), '  %.4f' % (z1Old) , '  %.4f' % (z3Old) , '  %.4f' % (negEntropyOld),   '  %.4f' % (enthalpy1Old),  '  %.4f' % (FEValueOld)       
            print ' New:   %.4f' % (x1New), '  %.4f' % (y2New), '  %.4f' % (z1New) , '  %.4f' % (z3New), '  %.4f' % (negEntropyNew),   '  %.4f' % (enthalpy1New),  '  %.4f' % (FEValueNew)    
            print ' ' 

        # If the flip was successful, keep the values (assign to unitArray)
        #  Otherwise, revert newUnitArray back to its earlier state
        successBool = 0
        if FEValueNew < FEValueOld: successBool = 1
        if successBool:
#            print ' ' 
#             print ' Original unitArray at Row = ', x1Row, ' Col = ', x1Col, ' is ', unitArray[x1Row,x1Col]
#             print ' Original unitArray at Row = ', x2Row, ' Col = ', x2Col, ' is ', unitArray[x2Row,x2Col]
            unitArray[x1Row,x1Col] = 0
            unitArray[x2Row,x2Col] = 1
#             print ' After flip:'
#             print '      New unitArray at Row = ', x1Row, ' Col = ', x1Col, ' is ', unitArray[x1Row,x1Col]
#             print '      New unitArray at Row = ', x2Row, ' Col = ', x2Col, ' is ', unitArray[x2Row,x2Col]
            successfulFlips = successfulFlips + 1
            x1End = x1New
            y2End = y2New
            z1End = z1New
            z3End = z3New
            negEntropyEnd = negEntropyNew
            enthalpy1End = enthalpy1New
            FEValueEnd = FEValueNew                           
            if not findFEMinimumValsDetailsBoolOff:                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
                print ' '
                print ' Successful flip: free energy reduced, keeping the change'         
        else: 
            newUnitArray[x1Row,x1Col] = 1
            newUnitArray[x2Row,x2Col] = 0 
            x1End = x1Old
            y2End = y2Old
            z1End = z1Old
            z3End = z3Old
            negEntropyEnd = negEntropyOld
            enthalpy1End = enthalpy1Old
            FEValueEnd = FEValueOld 
            if not findFEMinimumValsDetailsBoolOff:                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
                print ' '
                print ' Unsuccessful flip: free energy increased, NOT keeping the change'        

    print ' '
    print ' In adjustMatrixFEMinimum, completed all computations, for h = ', h
    print ' '
    print ' The total number of successful flips out of ', totalTrials, 'is: ', successfulFlips 
    print ' ' 
    print ' The total change in the configuration and thermodynamic variables:'
    print '          x1      y2       z1       z3       negS    enthalpy1   FE '
    print ' Orig:   %.4f' % (x1Orig), '  %.4f' % (y2Orig), '  %.4f' % (z1Orig) , '  %.4f' % (z3Orig) , '  %.4f' % (negEntropyOrig),   '  %.4f' % (enthalpy1Orig),  '  %.4f' % (FEValueOrig)       
    print ' End:    %.4f' % (x1End), '  %.4f' % (y2End), '  %.4f' % (z1End) , '  %.4f' % (z3End), '  %.4f' % (negEntropyEnd),   '  %.4f' % (enthalpy1End),  '  %.4f' % (FEValueEnd)    
    print ' '
    print ' '                                                                                                          

    if not findFEMinimumValsBoolOff:
        printUnitArrayModificationResults (x1ValsArray, y2ValsArray, 
            z1ValsArray, z3ValsArray, negSValsArray, enthalpy1Array, freeEnergyArray, h, totalTrials, findFEMinimumValsBoolOff)             
                                     
       
    return unitArray
    # END adjustMatrixFEMinimum      
                
                                                
####################################################################################################
#
# Function to compute a term (Lf(v)) needed to compute the entropy
#
####################################################################################################

def LfFunc(x):
    Lfval = x*log(x)-x
    return (Lfval)


####################################################################################################
####################################################################################################
#
# Function to compute the entropy
# This was used exclusively for computing the predicted (analytic) values for the entropy,
#  in the case where the configuration variables were analytically defined in terms of probability
#  distribution for h = 1. 
#
####################################################################################################
####################################################################################################


def computeAnalyticConfigAndThermVars (x1TargetVal):

#    print ' '
#    print '  In AnalyticConfigAndThermVars '
#    print ' '

    totalUnits = float(arrayLength*arrayLayers)
#    printConfigVars (configVarsList)      
                        
    x1 = x1TargetVal
    x2 = 1.0 - x1TargetVal
    sumX = x1+x2  

#    print ' ' 
#    print '                x1 = %.4f'  % x1
#    print '                x2 = %.4f'  % x2
#    print '   Sum of the x(i) = %.4f'  % sumX
#    print ' ' 
    
    y1 = x1*x1
    y2 = x1*x2
    y3 = x2*x2
    w1 = x1*x1
    w2 = x1*x2
    w3 = x2*x2        
    z1 = x1*x1*x1
    z2 = x1*x1*x2
    z3 = x1*x2*x1 
    z4 = x2*x1*x2
    z5 = x2*x2*x1
    z6 = x2*x2*x2

    sumY = y1+y2+y3    
#    print
#    print "Totals for the y(i) variables:"
#    print '          (A-A) y1 = %.4f'  % y1 
#    print '    (A-B & B-A) y2 = %.4f'  % y2   
#    print '          (B-B) y3 = %.4f'  % y3 
#    print '   Sum of the y(i) = %.4f'  % sumY
#    print 
    
    sumW = w1+w2+w3 
#    print
#    print "Totals for the w(i) variables:"
#    print '        (A---A) w1 = %.4f'  % w1 
#    print '(A---B & B---A) w2 = %.4f'  % w2   
#    print '        (B---B) w3 = %.4f'  % w3 
#    print '   Sum of the w(i) = %.4f'  % sumW
#    print         

    sumZ = z1+z2+z3 +z4+z5+z6

#    print
#    print "Totals for the z(i) variables:"
#    print '        (A-A-A) z1 = %.4f'  % z1 
#    print '(A-A-B & B-A-A) z2 = %.4f'  % z2   
#    print '        (A-B-A) z3 = %.4f'  % z3 
#    print '        (B-A-B) z4 = %.4f'  % z4 
#    print '(A-B-B & B-B-A) z5 = %.4f'  % z5 
#    print '        (B-B-B) z6 = %.4f'  % z6         
#    print '   Sum of the z(i) = %.4f'  % sumZ
#    print                    
                                                
    Lfx = LfFunc(x1)+LfFunc(x2)
#    print ' lnFunc (x) = %.4f'  % (Lfx)

    Lfy = LfFunc(y1) + 2*LfFunc(y2) + LfFunc(y3)
#    print ' lnFunc (y) = %.4f'  % (Lfy)  
      
    Lfw = LfFunc(w1) + 2*LfFunc(w2) + LfFunc(w3)
#    print ' lnFunc (w) = %.4f'  % (Lfw)

    Lfz = LfFunc(z1) + 2*LfFunc(z2) + LfFunc(z3) + LfFunc(z4) + 2*LfFunc(z5) + LfFunc(z6)
#    print ' lnFunc (z) = %.4f'  % (Lfz)
        
    s = -(2*Lfy+Lfw-Lfx-2*Lfz)
  
      
    return (y2, z1, z3, s)

    
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
    y2 = float(configVarsList[3])/(totalUnitsTimesTwo*2.0)
    y3 = float(configVarsList[4])/totalUnitsTimesTwo 
    
    sumY = y1 + 2.0*y2 + y3 # y1, y3 at equiprobable = 0.25
                        # y2 at equiprobable = 0.5
                        # That y2 = A-B + B-A.
                        # sumY should always = 1.0 

    w1 = float(configVarsList[5])/totalUnitsTimesTwo
    w2 = float(configVarsList[6])/(totalUnitsTimesTwo*2.0)
    w3 = float(configVarsList[7])/totalUnitsTimesTwo        

    sumW = w1 + 2.0*w2 + w3


    z1 = float(configVarsList[8])/totalUnitsTimesTwo
    z2 = float(configVarsList[9])/(totalUnitsTimesTwo*2.0)
    z3 = float(configVarsList[10])/totalUnitsTimesTwo 
    z4 = float(configVarsList[11])/totalUnitsTimesTwo
    z5 = float(configVarsList[12])/(totalUnitsTimesTwo*2.0)
    z6 = float(configVarsList[13])/totalUnitsTimesTwo 
              

    sumZ = z1 + 2.0*z2 + z3 + z4 + 2.0*z5 + z6

    # This is V&V code, to confirm that the sums of the different configuration variables is correct
    if not debugPrintOff:
        print ' '         
        print ' '
        print ' x1 = %.4f'  % (x1), '      x2 = %.4f'  % (x2)
        print ' y1 = %.4f'  % (y1), ' y2Total = %.4f'  % (2.*y2), ' y3 = %.4f'  % (y3) #, ' sumY = %.4f'  % (sumY)                                                 
        print ' w1 = %.4f'  % (w1), ' w2Total = %.4f'  % (2.0*w2), ' w3 = %.4f'  % (w3) #, ' sumW = %.4f'  % (sumW) 

        print ' '
        print ' z1 = %.4f'  % (z1), ' z2Total = %.4f'  % (2.0*z2), ' z3 = %.4f'  % (z3)
        print ' z6 = %.4f'  % (z6), ' z5Total = %.4f'  % (2.0*z5), ' z4 = %.4f'  % (z4)  #, 'sumZ = %.4f'  % (sumW)     
        
                                                                                                                                                    
    Lfx = LfFunc(x1)+LfFunc(x2)
#    print ' lnFunc (x) = %.4f'  % (Lfx)

    Lfy = LfFunc(y1) + 2.0*LfFunc(y2) + LfFunc(y3)
#    print ' lnFunc (y) = %.4f'  % (Lfy)  
      
    Lfw = LfFunc(w1) + 2.0*LfFunc(w2) + LfFunc(w3)
#    print ' lnFunc (w) = %.4f'  % (Lfw)

    Lfz = LfFunc(z1) + 2.0*LfFunc(z2) + LfFunc(z3) + LfFunc(z4) + 2.0*LfFunc(z5) + LfFunc(z6)
#    print ' lnFunc (z) = %.4f'  % (Lfz)
        
    negS = -(2.*Lfy+Lfw-Lfx-2.*Lfz)

   
#  If we were to attempt to derive h from the (equilibrium) configuration variables: 
#    h = ((z1*y2Actual)/(z3*y1))**(0.5)
#  However, we are conducting these experiments given a starting h-value

#  Notice, in the above equations, that z2 and z5 are both doubles; this allows sumZ = 1.0
#    However, we need to be clear that (e.g.) z2 = A-A-B + B-A-A as counted above, etc. 
#  For free energy  (entropy) calculations, we need to use the y2Actual, z2Actual, and z5Actual;
#    these have the degeneracy removed, and it shows up as the multiplying factor of 2 in front 
#   of the LfFunc. 

#  The extra multiplying factor of 2 in the entropy equation above (in front of the y and z terms)
#    comes from how the grid is assembled. See the 2014 Technical Report on the 2-D CVM for details.  
    
 
    
    epsilon1 = 4*log(h)         
    epsilon0 = 0.0    

    enthalpy0 = epsilon0*x1
    enthalpy1 = epsilon1*(2.0*y2-y1-y3)
    # Note that the actual equation has: Enthalpy = 2N*epsilon1*(-y1+2y2-y3)
    # Equivalently: Enthalpy = 2N*epsilon1*(z2+z4+z3+z5)  
    # I've simplified the Enthalpy to reflect only the y2 interactions, which in this 
    #   case would be the SINGLE y2 (A-B or B-A), so that
    #   2y2 = the y2 that is approximately 0.5, where y1 + y2 + y3 = 1.0    
    
    freeEnergy = enthalpy0 + enthalpy1 + negS          

        
    sysValsList = (negS, enthalpy0, enthalpy1, freeEnergy)

    if not debugPrintOff:
        print ' ' 
        print ' '
        print ' The computed thermodynamic quantities:'            
        print '    h = %.4f' % (h), '    epsilon1 = %.4f' % (epsilon1) 
        print '    negative entropy s = %.4f' % (negS)
        print '    enthalpy = %.4f' % (enthalpy1)      
        print '    free energy = %.4f' % (freeEnergy)
        print ' ' 
        print ' '   
                                  
    return (sysValsList)   
         
    
####################################################################################################
####################################################################################################
#
# Function to compute a set of configuration and thermodynamic variables given specific target
#  values for x1, along with parameter h. 
#
# For prior experiments investigating how the thermodynamic variables behaved in correlation with
#  the configuration variables, I use a loop to repeat the same calculations for a "numTrials" 
#  total set of runs; the averages were then reported for a given set of x, h values. 
#
# For the perturbation experiments, since the goal is to work with a unitArray that is perturbed,
#  and then allowed to cascade again into a free energy minimum, the original FOR loops are being 
#  kept, but numTrials is assigned (in **main**) to have a value of 1. 
#
####################################################################################################
####################################################################################################


def computeConfigAndThermVars(arraySizeList, h, x1TargetVal, maxXDif, numTrials, jrange, maxRange):

    totalUnits = float(arrayLength*arrayLayers)
    totalUnitsTimesTwo = totalUnits*2.0
    
 
    xArray = np.zeros(numTrials, dtype=np.float)
    x1Array = np.zeros(numTrials, dtype=np.float)  
    y1Array = np.zeros(numTrials, dtype=np.float)
    y2Array = np.zeros(numTrials, dtype=np.float)        
    y3Array = np.zeros(numTrials, dtype=np.float)
    w1Array = np.zeros(numTrials, dtype=np.float)
    w2Array = np.zeros(numTrials, dtype=np.float)        
    w3Array = np.zeros(numTrials, dtype=np.float)        
    z1Array = np.zeros(numTrials, dtype=np.float)
    z2Array = np.zeros(numTrials, dtype=np.float) 
    z3Array = np.zeros(numTrials, dtype=np.float)
    z4Array = np.zeros(numTrials, dtype=np.float)     
    z5Array = np.zeros(numTrials, dtype=np.float)
    z6Array = np.zeros(numTrials, dtype=np.float)         
    negSArray = np.zeros(numTrials, dtype=np.float)
    EnthEps0Array = np.zeros(numTrials, dtype=np.float)         
    EnthEps1Array = np.zeros(numTrials, dtype=np.float)  
    freeEnergyArray = np.zeros(numTrials, dtype=np.float) 

    if not debugPrintOff:
        print ' '     
        print ' The thermodynamic quantities plot vs. x1' 
        print '   The negEntropy is in blue, '
        print '   The y2 value is in maroon.'
        print '   The per-unit enthalpy is non-zero, but cannot be determined at this time, '
        print '   The interaction enthalpy (eps1*y2) is set to zero (eps1 = 0, h = 1)'
   

    sumx1 = 0.0
    sumy1 = 0.0
    sumy2 = 0.0
    sumy3 = 0.0
    sumw1 = 0.0
    sumw2 = 0.0
    sumw3 = 0.0    
    sumz1 = 0.0
    sumz2 = 0.0 
    sumz3 = 0.0
    sumz4 = 0.0     
    sumz5 = 0.0
    sumz6 = 0.0           
    sumNegS = 0.0      
    sumEnthEps0   = 0.0
    sumEnthEps1   = 0.0
    sumFreeEnergy = 0.0

    # This loops to obtain average values for x1, y2, and the negative entropy for 
    #  a specific target x1.    
    for i in range (0, numTrials, 1):        

        configVarsList = list ()  # empty list
        unitArray    = initializeMatrix (arraySizeList, h, x1TargetVal, maxXDif)        
        # debug print: 
        #    print unitArray  

        # adjust the unit array so that the actual x1 approximately = target x1
        unitArray = adjustMatrix (arraySizeList, unitArray, h, jrange, maxXDif, x1TargetVal, 
        beforeAndAfterAdjustedMatrixPrintOff)

        # adjust the unit array so that the configuration variables take values yielding a
        #  free energy minimum for the given h-value        
        unitArrayForFEMinimum = adjustMatrixFEMinimum (arraySizeList, unitArray, h, maxRange)
        

# Obtain the configuration variables once again (this step was also done while adjusting the array) 
# The configuration variables list returned from computeConfigVariables is the list of 
#   RAW COUNTS, meaning that the:
#   -  x vars need to be divided by the total number of units
#   -  y, w, and z vars need to be divided by the 2*total number of units, and
#   -  y2, w2, z2, and z5 vars need to be divided again by 2 to reach their normative values,
#      so that (e.g.), y1 + 2y2 + y3 = 0, etc. 
      
        configVarsList = computeConfigVariables (arraySizeList, unitArrayForFEMinimum)
                                                                
# obtain the thermodynamic variables 
        sysVarsList = computeThermodynamicVars (h, configVarsList)   

# store the thermodynamic variables to plot later
        x1     = float(configVarsList[0])/totalUnits
        y1     = float(configVarsList[2])/totalUnitsTimesTwo          
        y2Real = float(configVarsList[3])/(totalUnitsTimesTwo*2.0)        
        y3     = float(configVarsList[4])/totalUnitsTimesTwo
        w1     = float(configVarsList[5])/totalUnitsTimesTwo          
        w2Real = float(configVarsList[6])/(totalUnitsTimesTwo*2.0)         
        w3     = float(configVarsList[7])/totalUnitsTimesTwo          
        z1     = float(configVarsList[8])/totalUnitsTimesTwo  
        z2Real = float(configVarsList[9])/(totalUnitsTimesTwo*2.0) 
        z3     = float(configVarsList[10])/totalUnitsTimesTwo  
        z4     = float(configVarsList[11])/totalUnitsTimesTwo 
        z5Real = float(configVarsList[12])/(totalUnitsTimesTwo*2.0)   
        z6     = float(configVarsList[13])/totalUnitsTimesTwo 
                                        
        xArray[i]       = i  # x1
        x1Array[i]      = x1  # x1
        y1Array[i]      = y1
        y2Array[i]      = y2Real  # y2
        y3Array[i]      = y3
        w1Array[i]      = w1
        w2Array[i]      = w2Real  # w2
        w3Array[i]      = w3        
        z1Array[i]      = z1  # z1
        z2Array[i]      = z2Real  # z2        
        z3Array[i]      = z3  # z3
        z4Array[i]      = z4  # z4        
        z5Array[i]      = z5Real  # z5
        z6Array[i]      = z6  # z6
        # sysValsList = (negS, enthalpy0, enthalpy1, freeEnergy)
        negSArray[i]    = sysVarsList[0] 
        EnthEps0Array[i]   = sysVarsList[1]         
        EnthEps1Array[i]   = sysVarsList[2] 
        freeEnergyArray[i] = sysVarsList[3]         
#        fEps0Array[i]   = sysVarsList[1]
#        fEps1Array[i]   = sysVarsList[2]-0.6  
#        fEnergyArray[i] = sysVarsList[3]-0.6                                                
        sumx1 = sumx1 + x1
        sumy1 = sumy1 + y1
        sumy2 = sumy2 + y2Real
        sumy3 = sumy3 + y3
        sumw1 = sumw1 + w1
        sumw2 = sumw2 + w2Real
        sumw3 = sumw3 + w3
        sumz1 = sumz1 + z1
        sumz2 = sumz2 + z2Real  
        sumz3 = sumz3 + z3
        sumz4 = sumz4 + z4 
        sumz5 = sumz5 + z5Real
        sumz6 = sumz6 + z6                       
        sumNegS = sumNegS + sysVarsList[0] 
        sumEnthEps0   = sumEnthEps0 + sysVarsList[1]
        sumEnthEps1   = sumEnthEps1 + sysVarsList[2]
        sumFreeEnergy = sumFreeEnergy + sysVarsList[3]                                                                                                                     
# END: FOR loop (to compute x1, y2, and negS for a given target x1)

# Compute the average values for various values
    denom = float(numTrials)
    avgx1 = sumx1/denom
    avgy1 = sumy1/denom
    avgy2 = sumy2/denom
    avgy3 = sumy3/denom    
    avgw1 = sumw1/denom
    avgw2 = sumw2/denom
    avgw3 = sumw3/denom  
    avgz1 = sumz1/denom
    avgz2 = sumz2/denom    
    avgz3 = sumz3/denom
    avgz4 = sumz4/denom 
    avgz5 = sumz5/denom
    avgz6 = sumz6/denom     
    avgNegS = sumNegS/denom
    avgEnthEps0 = sumEnthEps0/denom
    avgEnthEps1 = sumEnthEps1/denom    
    avgFreeEnergy = sumFreeEnergy/denom
    
    if not debugPrintOff:
        print ' '   

    # Print the average results
        print ' '
        print ' For a target x1 of: %.4f' % (x1TargetVal) 
        print '   over a set of ', numTrials, ' trials, the average values are:'
        print '          x1: %.4f' % (avgx1)     
        print '          y1: %.4f' % (avgy1) 
        print '          y2: %.4f' % (avgy2) 
        print '          z1: %.4f' % (avgz1)
        print '          z2: %.4f' % (avgz2) 
        print '          z3: %.4f' % (avgz3)
        print '          z4: %.4f' % (avgz4)
        print '          z5: %.4f' % (avgz5)  
        print '          z6: %.4f' % (avgz6)                                  
        print '         negS: %.4f' % (avgNegS)   
        print '    Enthalpy0: %.4f' % (avgEnthEps0) 
        print '    Enthalpy1: %.4f' % (avgEnthEps1) 
        print '    Free Engy: %.4f' % (avgFreeEnergy)     
    # Plot the results from that FOR loop (multiple tests of a random grid for a given x1 value)          
    #    print '     and the free energy is in red, also shifted by 0.6.'   
        pylab.figure(1)
        pylab.plot (xArray,negSArray)    
        pylab.plot (xArray, y2Array, 'm')
    # #   pylab.plot (xArray,fEps1Array,'m')
    # #   pylab.plot (xArray,fEnergyArray,'r')
        pylab.show()  
    # END: FOR loop (to test for a range of target x1 values) 
 
    
            
    newList = (avgx1, avgy1, avgy2, avgy3, avgw1, avgw2, avgw3, avgz1, avgz2, avgz3, avgz4, avgz5, avgz6, 
               avgNegS, avgEnthEps0, avgEnthEps1, avgFreeEnergy, unitArrayForFEMinimum)
    

                                                                    
    return (newList) 
    # END - return a set of values for a single trial with x1TargetVal   



####################################################################################################
####################################################################################################
#
#
# Function to create a perturbed unitArray, where the fraction perturbFrctn units are flipped
#
#
####################################################################################################
####################################################################################################

def perturb (arraySizeList, unitArray, perturbFrctn):
    
    localArrayLength = arraySizeList[0]
    localArrayLayers = arraySizeList[1]
    
    flipYesProb = perturbFrctn
    flipNoProb  = 1.0 - flipYesProb
    
    perturbedUnitArray = createIdenticalUnitArray (arraySizeList, unitArray)

    expectedFlips = int(perturbFrctn*localArrayLength*localArrayLayers)    
    totalFlips = 0
                                                                                                                                                                                                                                                                                                                                                                                                  
    for i in range (0, localArrayLayers):
        for j in range(0,localArrayLength):
            flipYes = np.random.choice([0, 1], p=[flipNoProb, flipYesProb])   
            if flipYes ==1:
                totalFlips = totalFlips + 1
                if perturbedUnitArray[i,j] ==1: 
                   perturbedUnitArray[i,j] = 0
                else: perturbedUnitArray[i,j] = 1   
    print ' '
    print ' The total number of flips in creating the perturbed unitArray is: ', totalFlips
    print ' With a perturbation fraction of print %.2f' % (perturbFrctn), 'the expected total number of flips was: ', expectedFlips
    print ' '   
    
    return (perturbedUnitArray)



####################################################################################################
####################################################################################################
#
# Procedure to create a visually more intuitive depiction of the "1" and "0" population of a grid.
#   A "1" is shown as an "X" in the print, and a "0" as a "-". 
#
#    Inputs:    arraySizeList, unitArray
#
#
####################################################################################################
####################################################################################################

def prettyPrintArray (arraySizeList, unitArray):

    localArrayLength = arraySizeList[0]
    localArrayLayers = arraySizeList[1]

    # NOTE: pairs is the number of pairs of layers, defined in **main**
    #       as pairs = int(arrayLayers/2.0)
    
    print ' '        
    for i in range (0, pairs):
        actualEvenRowNum = 2*i
        print 'Row', actualEvenRowNum, ':', blnkspc, 
        for j in range(0,localArrayLength):
            if unitArray[actualEvenRowNum,j] ==1:
                print 'X', blnkspc,
            else: print '-', blnkspc,
        print 
        actualOddRowNum = 2*i+1
        print 'Row ', actualOddRowNum, ':', blnkspc,
        print (blnkspc),
        for j in range(0,localArrayLength):
            if unitArray[actualOddRowNum,j] ==1:
                print 'X', blnkspc,
            else: print '-', blnkspc,            
        print ' '
                
                

####################################################################################################
####################################################################################################
#
# Procedure to create a visually more intuitive depiction of the changes in the population of a 
#   perturbed unitArray.
#   A change from 0 to 1 is shown as an "W", and a change from a 1 to a 0 is shown as an "M"  
#
#    Inputs:    arraySizeList, unitArray, perturbedUnitArray
#
#
####################################################################################################
####################################################################################################

def prettyPrintChangesInUnitArray (arraySizeList, unitArray1, unitArray2):

    localArrayLength = arraySizeList[0]
    localArrayLayers = arraySizeList[1]

    # NOTE: pairs is the number of pairs of layers, defined in **main**
    #       as pairs = int(arrayLayers/2.0)
    
    totalChanges = 0
    print ' '        
    for i in range (0, pairs):
        actualEvenRowNum = 2*i
        print 'Row', actualEvenRowNum, ':', blnkspc, 
        for j in range(0,localArrayLength):
            if unitArray1[actualEvenRowNum,j] == 1:
                if unitArray2[actualEvenRowNum,j] == 1:
                    print '-', blnkspc,
                else: 
                    print 'M', blnkspc,
                    totalChanges = totalChanges + 1
            else: # unitArray element == 0
                if unitArray2[actualEvenRowNum,j] == 1:
                    print 'W', blnkspc,
                    totalChanges = totalChanges + 1  
                else: print '-', blnkspc,                 
        print 
        actualOddRowNum = 2*i+1
        print 'Row ', actualOddRowNum, ':', blnkspc,
        print (blnkspc),
        for j in range(0,localArrayLength):
            if unitArray1[actualOddRowNum,j] ==1:
                if unitArray2[actualOddRowNum,j] == 1:
                    print '-', blnkspc,
                else: 
                    print 'M', blnkspc,
                    totalChanges = totalChanges + 1                    
            else: # unitArray element == 0
                if unitArray2[actualOddRowNum,j] == 1:
                    print 'W', blnkspc,
                    totalChanges = totalChanges + 1                      
                else: 
                    print '-', blnkspc,                                               
        print ' '
    print ' Total Changes = ', totalChanges                



####################################################################################################
####################################################################################################
#
# Procedure to print out the final set of configuration variables x'i, y'i, w'i, and z'i
#   This procedure was used for early V&V, to ensure that the configuration variables were
#     as expected; 
#   It is not called in the current configuration of this program.
#
#
####################################################################################################
####################################################################################################


def printConfigVarsComparison(h, configVarsList, configVarsListDelta):

# This print function is not used in this current version of the code;
#   it was handy during earlier debug stages.    
            
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
####################################################################################################
#
# Procedure to print the results of adapting a unitArray so that its configuration (for given 
#   values of x1 & h) achieve free energy minimum. 
#
####################################################################################################
####################################################################################################


def printUnitArrayModificationResults (x1ValsArray, y2ValsArray, 
        z1ValsArray, z3ValsArray, negSValsArray, enthalpy1Array, freeEnergyArray, h, totalTrials, 
        findFEMinimumValsBoolOff):  


    trialNumArray = np.zeros(totalTrials, dtype=np.float)
    
    print ' ' 
    print ' in printUnitArrayModificationResults'
    print '  h = %.3f '  % (h)
        
    print ' Final results, for the case where the starting x1 value is  %.4f' % x1ValsArray[0] 

    print ' ' 
    print '     x1     y2       z1      z3       negS    enth1  freeEnergy'


    for j in range (0, totalTrials, 1):
        trialNumArray[j] = j
        if not findFEMinimumValsBoolOff:               
            print '   %.4f' % (x1ValsArray[j]), '  %.4f' % (y2ValsArray[j]), '  %.4f' % (z1ValsArray[j]), '  %.4f' % (z3ValsArray[j]), '  %.4f' % (negSValsArray[j]), '  %.4f' % (enthalpy1Array[j]), '  %.4f' % (freeEnergyArray[j])          

    # Plot the results from that FOR loop (adjusting the unitArray to achieve free energy minimization)                                             
                                                                                                
    pylab.figure(1)
    pylab.plot (trialNumArray, negSValsArray)          
    pylab.plot (trialNumArray, y2ValsArray-0.8, 'g')     
    pylab.plot (trialNumArray, enthalpy1Array-0.5, 'r')
    pylab.plot (trialNumArray, freeEnergyArray, 'k')
           
    # #   pylab.plot (xArray,fEps1Array,'m')
    # #   pylab.plot (xArray,fEnergyArray,'r')
    print ' ' 
    print ' This plot shows results of adapting a unitArray to minimize free energy.'
    print ' '
    print '  The probabilistic negative entropy is in blue,'
    print '  The probabilistic y2 values (minus 0.8) are in green,' 
    print '  The probabilistic enthalpy(1) (minus 0.5) values are in red. ' 
    print '  The probabilistic free energy values are in black. '        
    pylab.show()                          
                                                
#   END printUnitArrayModificationResults 
                                                                                                                   

####################################################################################################
####################################################################################################
#
# Procedure to print the details of the entropy results; this is useful for understanding the
#   various contributions to the entropy term, especially for cases where h is not equal to 0. 
#
####################################################################################################
####################################################################################################


def printEntropyDetailedResults (x1ValsArray, y1ValsArray, y2ValsArray, y3ValsArray, w1ValsArray, w2ValsArray, w3ValsArray, 
    z1ValsArray, z2ValsArray, z3ValsArray, z4ValsArray, z5ValsArray, z6ValsArray, negSValsArray, enthalpy0Array, 
    enthalpy1Array, freeEnergyArray, x1TotalSteps, step, x1StartingVal, x1TargetIncrement, hArray, hTotalSteps, hStep):  

    x1TargetVal = x1StartingVal
    x1TargetValsArray = np.zeros(x1TotalSteps, dtype=np.float)
    hStart = hArray[0]
            
    print ' ' 
    print ' Entropy results, for the case where:'
    print '  - the starting x1 value is    %.4f' % (x1StartingVal - x1TargetIncrement)
    print '  - and actual x1 value is      %.4f' % (x1ValsArray[0])
    print '  - and the starting h value is %.2f' % (hStart) 
    print ' ' 
    print ' Entropy data for values of h:'
    print ' ' 
    print '     h      xNegS    yNegS     wNegS    zNegS     negS    '               
      

    # Plot the results from that FOR loop (multiple tests of a random grid for a given x1 value)          
    #    print '     and the free energy is in red, also shifted by 0.6.'       

    x2ValsArray     = np.zeros(hTotalSteps, dtype=np.float)
    xNegSArray      = np.zeros(hTotalSteps, dtype=np.float)
    yNegSArray      = np.zeros(hTotalSteps, dtype=np.float)
    wNegSArray      = np.zeros(hTotalSteps, dtype=np.float)
    zNegSArray      = np.zeros(hTotalSteps, dtype=np.float)
        
    for k in range (0, hTotalSteps, hStep):
        x2ValsArray[k]   = 1.0 - x1ValsArray[k]
        xNegSArray[k] = x1ValsArray[k]*log(x1ValsArray[k]) + x2ValsArray[k]*log(x2ValsArray[k])       
        yNegSArray[k]  = y1ValsArray[k]*log(y1ValsArray[k]) + 2.*y2ValsArray[k]*log(y2ValsArray[k]) + y3ValsArray[k]*log(y3ValsArray[k])
        wNegSArray[k]  = w1ValsArray[k]*log(w1ValsArray[k]) + 2.*w2ValsArray[k]*log(w2ValsArray[k]) + w3ValsArray[k]*log(w3ValsArray[k])
        zNegSArray[k]  = z1ValsArray[k]*log(z1ValsArray[k]) + 2.*z2ValsArray[k]*log(z2ValsArray[k]) + z3ValsArray[k]*log(z3ValsArray[k]) \
                       + z4ValsArray[k]*log(z4ValsArray[k]) + 2.*z5ValsArray[k]*log(z5ValsArray[k]) + z6ValsArray[k]*log(z6ValsArray[k]) 
        print '   %.2f' % (hArray[k]), '  %.4f' % (xNegSArray[k]),  '  %.4f' % (yNegSArray[k]), ' %.4f' % (wNegSArray[k]) , ' %.4f' % (zNegSArray[k]),' %.4f' % (negSValsArray[k])  
        

    # Plots have been turned off for this version of code, as it is specialized to deal with perturbation analysis for 
    #  a single h-value                                                                          
#    pylab.figure(1)
#    pylab.plot (hArray, negSValsArray)          
#    pylab.plot (hArray, xNegSArray, 'g')
#    pylab.plot (hArray, yNegSArray, 'c')
#    pylab.plot (hArray, wNegSArray, 'm')  
#    pylab.plot (hArray, zNegSArray, 'r')             

#    print ' ' 
#    print ' ' 
#    print ' Plot of entropy details and the total negative entropy vs. h'
#    print ' '
#    print ' This plot contains the same information as shown in the previous table.'
#    print ' '
#    print '    (When h = 0, we have the only value for which probabilistic results will correspond '
#    print '      with the randomly-generated (and modified) array.)'
#    print ' '
#    print '  The probabilistic negative entropy is in blue,'
#    print '  The probabilistic negative x entropy term  is in green,' 
#    print '  The probabilistic negative y entropy term  is in cyan, '
#    print '  The probabilistic negative w entropy term  is in maroon,' 
#    print '  The probabilistic negative z entropy term  is in red.'  
#    print ' ' 
#    print '  The total negEntropy is given as: -(2*negYS + negWS - negXS - 2*negZS).'
#   from computeThermodynamicVars:  negS = -(2.*Lfy+Lfw-Lfx-2.*Lfz)  
#    print ' '            
#    pylab.show()                          
                                                                    

####################################################################################################
####################################################################################################
#
# Procedure to print final results, for cases where h is not equal to 0. 
#
####################################################################################################
####################################################################################################


def printResults (x1ValsArray, y1ValsArray, y2ValsArray, y3ValsArray, z1ValsArray, z3ValsArray, negSValsArray, enthalpy0Array, 
    enthalpy1Array, freeEnergyArray, x1TotalSteps, step, x1StartingVal, x1TargetIncrement, hArray, hTotalSteps, hStep):  

    x1TargetVal = x1StartingVal
    x1TargetValsArray = np.zeros(x1TotalSteps, dtype=np.float)
    hStart = hArray[0]
    
    print ' '         
    print ' Configuration variable and thermodynamic results for the case where:'
    print '  the starting x1 value is    %.4f' % (x1StartingVal - x1TargetIncrement)
    print '  and actual x1 value is      %.4f' % (x1ValsArray[0])
    print '  and the starting h value is %.2f' % (hStart) 
    print ' ' 
    print ' Data for varying values of h:'
    print ' ' 
    print '     h      y2       z1      z3       negS    delta   epsilon   enth1  freeEnergy'


    deltaArray      = np.zeros(hTotalSteps, dtype=np.float)
    epsilon1Array   = np.zeros(hTotalSteps, dtype=np.float)

    for k in range (0, hTotalSteps, hStep):
        deltaArray[k] = 2.0*y2ValsArray[k] - y1ValsArray[k] - y3ValsArray[k] 
        epsilon1Array[k] = 4.* log(hArray[k])

    for j in range (0, hTotalSteps, hStep):                
        print '   %.2f' % (hArray[j]), '  %.4f' % (y2ValsArray[j]), '  %.4f' % (z1ValsArray[j]), '  %.4f' % (z3ValsArray[j]), '  %.4f' % (negSValsArray[j]), '  %.4f' % (deltaArray[j]), '  %.4f' % (epsilon1Array[j]), '  %.4f' % (enthalpy1Array[j]), '  %.4f' % (freeEnergyArray[j])          

    # Plot the results from that FOR loop (multiple tests of a random grid for a given x1 value)          
    #    print '     and the free energy is in red, also shifted by 0.6.'       

    # Plots have been turned off for this version of code, as it is specialized to deal with perturbation analysis for 
    #  a single h-value    

                        
#    pylab.figure(1)
#    pylab.plot (hArray, negSValsArray)          
#    pylab.plot (hArray, y2ValsArray, 'g')
#    pylab.plot (hArray, deltaArray, 'c')      
#    pylab.plot (hArray, z3ValsArray, 'm')
#    pylab.plot (hArray, z1ValsArray, 'y')
                                
#    pylab.figure(2)
#    pylab.plot (hArray, negSValsArray)          
#    pylab.plot (hArray, y2ValsArray-0.8, 'g')
#    pylab.plot (hArray, deltaArray-0.5, 'c')      
#    pylab.plot (hArray, enthalpy1Array-0.5, 'r')
#    pylab.plot (hArray, freeEnergyArray, 'k')

#    print ' ' 
#    print ' '
#    print ' Figure 1:' 
#    print ' Plot of config. variables and the negative entropy vs. h'
#    print ' '
#    print ' This plot contains the same information as shown in the previous table.'
#    print ' '
#    print '  The probabilistic negative entropy is in blue,'
#    print '  The probabilistic y2 values are in green,' 
#    print '  The probabilistic delta 2.0*y2 - y1 - y3 values are in cyan. '
#    print '  The probabilistic z3 values are in maroon. ' 
#    print '  The probabilistic z1 values are in yellow. ' 
#    print ' '            
                      
                                 
                                                       
#    print ' ' 
#    print ' '
#    print ' Figure 2:' 
#    print ' Plot of deltas in y config. variables, entropy, enthalpy1, and free energy vs. h'
#    print ' '
#    print ' This plot contains the same information as shown in the previous table.'
#    print ' '
#    print '    (This is only value for which probabilistic results will correspond '
#    print '    with the randomly-generated (and modified) array.)'
#    print ' '
#    print '  The probabilistic negative entropy is in blue,'
#    print '  The probabilistic y2 values (minus 0.8) are in green,' 
#    print '  The probabilistic delta 2.0*y2 - y1 - y3 values (minus 0.5) are in cyan. '
#    print '  The probabilistic enthalpy(1) (minus 0.5) values are in red. ' 
#    print '  The probabilistic free energy values are in black. ' 
#    print ' ' 
           
#    pylab.show()                          
                                                
                                                                                                
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
    global detailedAdjustMatrixPrintOff
    global beforeAndAfterAdjustedMatrixPrintOff
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


# Determine if we are probabilistically generating a pattern (patternProb = 0)
#   or if we are selecting a pre-stored pattern (patternProb = 1 ... N)   
    x1TargetVal = 0.5 
    maxXDif = 0.01    # maximum difference between x1 from randomly-generated array and the
                       # desired value of x1 = x1TargetVal  
    jrange = 200        # maximal number of steps allowed to improve the x1 distribution

    blnkspc=' '

    if arrayLayers % 2 == 0: evenLayers == True #then an even number of layers

    # Determine the total number of PAIRS of zigzag chains
    pairsLayers = arrayLayers/2
    pairs = int(pairsLayers) 

    debugPrintOff = True
    detailedDebugPrintOff = True
    ZDebugPrintOff = True
    detailedAdjustMatrixPrintOff = True

# This is a local variable; it will be passed to computeConfigVariables
#  It will determine whether we print the contents of the x-array at the
#  beginning and end of the adjust-matrix step. 
    beforeAndAfterAdjustedMatrixPrintOff = True
        
    if not debugPrintOff:
        print ' '
        print 'Debug printing is on'  #debugPrintOff false
    else: 
        print 'Debug printing is off' #debugPrintOff true    
            

####################################################################################################
# Compute configuration variables
####################################################################################################                
      
                        
    eps0 = 0.0  # The enthalpy for single unit activation is set to zero for this code


# Pick the increment and starting value for decreasing x1 (from x1 = 0.5 ... x1Final)
# Note that the smallest workable value for x1 is around x1 = 0.25;
#   typical test ranges are from x1 = 0.5 ... 0.25, stepping 0.5 (incr = step/20.0); step = 1
    
    x1TotalSteps = 1
    step = 1     
    x1StartingVal = 0.40
    x1TargetVal = x1StartingVal
    x1TargetIncrement = float(step)/20.0     
    numTrials = 1
    
    maxRange =  200 # This controls the number of iterations used to reduce free energy,
                   #  once an adjusted unit matrix has been formed with the desired
                   #  x1, x2 values

    x1ValsArray   = np.zeros(x1TotalSteps, dtype=np.float)
    y2ValsArray   = np.zeros(x1TotalSteps, dtype=np.float)
    z1ValsArray   = np.zeros(x1TotalSteps, dtype=np.float)        
    negSValsArray = np.zeros(x1TotalSteps, dtype=np.float) 
       
       

# This entire loop and the following print procedure is a V&V test to confirm that the 
#   the observed results for a probabilistically-generated unitArray (subject to 
#   selection of x1 = 0.5 .. 0.25, where 0.25 is the smallest possible value for 
#   x1 before causing program crash) conforms with the probabilistic expectations for
#   the configuration values given a starting value for x1. 


# This particular sequence cannot be used for different valeus of h, as the functional loop
#   as well as the print procedure called at the end each compare the analytic probabilistic
#   results with the ones obtained from the unitArray. 

#    h = 1.0    
    # This loops finds the configuration variables and the entropy for various values of x1    
#    for j in range (0, x1TotalSteps, step):     
#        x1TargetVal = x1TargetVal - x1TargetIncrement
        
#        if not debugPrintOff:
#            print ' ' 
#            print '  Entering the FOR loop to test different x1 target values'
#            print '    The current x1TargetVal is: %.4f' % (x1TargetVal) 

#        newArrayList = computeConfigAndThermVars(arraySizeList, h, x1TargetVal, maxXDif, numTrials, jrange, maxRange)
#        x1ValsArray[j]   = newArrayList[0]
#        y2ValsArray[j]   = newArrayList[1]
#        z1ValsArray[j]   = newArrayList[2]
#        negSValsArray[j] = newArrayList[3]                                                                                                                

# For the perturbation experiments, this code was modified (01/05/2018) 
#   so that it ONLY worked with a single choice of x and h, because the 
#   perturbation was expected to take some time ... and because I wanted 
#   full prints of the initial unitArray, the perturbArray, and the 
#   finalArray (after the perturbArray had been allowed to come to 
#   free energy minimum once again. 

# I've removed the h=1 probability and analytic results procedures, calls, and
#   test booleans. 
             
                                       

# This sequence expands on the previous sequence, by allowing for multiple values of h to be 
#   tested for a given x1. This introduces an h-loop into the previous x1-loop. 
#   When h is not equal to 1, we cannot compare our results to a probabilistic generation of 
#   configuration variables. Instead, we are now relying completely on the free energy minimization. 

    # This loops selects a (iteratively-varying) value for x1    
    
    # Redefine certain parameters for running unit tests
    x1TotalSteps = 1
   
    # Define parameters for going through the h-increment FOR loop
    hTotalSteps = 1
    hStep = 1
    hIncrement = 0.1

    x1ValsArray   = np.zeros(hTotalSteps, dtype=np.float)
    y1ValsArray   = np.zeros(hTotalSteps, dtype=np.float)
    y2ValsArray   = np.zeros(hTotalSteps, dtype=np.float)
    y3ValsArray   = np.zeros(hTotalSteps, dtype=np.float)
    w1ValsArray   = np.zeros(hTotalSteps, dtype=np.float)
    w2ValsArray   = np.zeros(hTotalSteps, dtype=np.float)
    w3ValsArray   = np.zeros(hTotalSteps, dtype=np.float)   
    z1ValsArray   = np.zeros(hTotalSteps, dtype=np.float) 
    z2ValsArray   = np.zeros(hTotalSteps, dtype=np.float)       
    z3ValsArray   = np.zeros(hTotalSteps, dtype=np.float) 
    z4ValsArray   = np.zeros(hTotalSteps, dtype=np.float)  
    z5ValsArray   = np.zeros(hTotalSteps, dtype=np.float) 
    z6ValsArray   = np.zeros(hTotalSteps, dtype=np.float)      
    negSValsArray = np.zeros(hTotalSteps, dtype=np.float)
    enthalpy0Array= np.zeros(hTotalSteps, dtype=np.float)
    enthalpy1Array= np.zeros(hTotalSteps, dtype=np.float)
    freeEnergyArray=np.zeros(hTotalSteps, dtype=np.float)
    hArray         =np.zeros(hTotalSteps, dtype=np.float)
    unitArray        = np.zeros(shape=(arrayLayers,arrayLength))
    perturbUnitArray = np.zeros(shape=(arrayLayers,arrayLength))   
    
    print ' ' 
    print ' In **main**, about to start the loop to test for both'
    print '   different x1 values and a range of h values.'
    print ' ' 
    
    for j in range (0, x1TotalSteps, step):     
        x1TargetVal = x1TargetVal - x1TargetIncrement
        
#        if not debugPrintOff:
        print ' ' 
        print '  In **main**, for step ', j, 'out of ', x1TotalSteps , ' total steps'
        print '  In the FOR loop to test different x1 target values'
        print '    The current x1TargetVal is: %.4f' % (x1TargetVal) 
        print '  About to start the loop through h values ' 

        hInitial = 1.4
        h = hInitial - hIncrement
        for hVal in range (0, hTotalSteps, hStep):
            h = h + hIncrement
            print ' '
            print '    In the h loop with h = ', h
            # newList = (avgx1, avgy2, avgz1, avgz3, avgNegS, avgEnthEps0, avgEnthEps1, avgFreeEnergy)
            newArrayList = computeConfigAndThermVars(arraySizeList, h, x1TargetVal, maxXDif, numTrials, jrange, maxRange)
            x1ValsArray[hVal]    = newArrayList[0]
            y1ValsArray[hVal]    = newArrayList[1]
            y2ValsArray[hVal]    = newArrayList[2]
            y3ValsArray[hVal]    = newArrayList[3]            
            w1ValsArray[hVal]    = newArrayList[4]
            w2ValsArray[hVal]    = newArrayList[5]
            w3ValsArray[hVal]    = newArrayList[6]
            z1ValsArray[hVal]    = newArrayList[7]
            z2ValsArray[hVal]    = newArrayList[8]    
            z3ValsArray[hVal]    = newArrayList[9]
            z4ValsArray[hVal]    = newArrayList[10] 
            z5ValsArray[hVal]    = newArrayList[11]
            z6ValsArray[hVal]    = newArrayList[12]             
            negSValsArray[hVal]  = newArrayList[13] 
            enthalpy0Array[hVal] = newArrayList[14] 
            enthalpy1Array[hVal] = newArrayList[15] 
            freeEnergyArray[hVal] = newArrayList[16]
            unitArray     = newArrayList[17]  
            hArray[hVal]          = h 
            # debug prints turned off
            # print ' ' 
            # print ' h = ', h, ' hVal (iteration number) = ', hVal, ' hArray[hVal] = ', hArray[hVal]                          

    
# This procedure prints out detailed components of the entropy term; due to x, y, w, and z configuration variables
    printEntropyDetailedResults (x1ValsArray, y1ValsArray, y2ValsArray, y3ValsArray, w1ValsArray, w2ValsArray, w3ValsArray,
            z1ValsArray, z2ValsArray, z3ValsArray, z4ValsArray, z5ValsArray, z6ValsArray, 
            negSValsArray, enthalpy0Array, enthalpy1Array, freeEnergyArray,           
            x1TotalSteps, step, x1StartingVal, x1TargetIncrement, hArray, hTotalSteps, hStep)

# This procedure prints out the results of both configuration and thermodynamic variables. 
#  (I plan to split it into two procedure calls later.) 
    printResults (x1ValsArray, y1ValsArray, y2ValsArray, y3ValsArray, z1ValsArray, z3ValsArray,
            negSValsArray, enthalpy0Array, enthalpy1Array, freeEnergyArray,           
            x1TotalSteps, step, x1StartingVal, x1TargetIncrement, hArray, hTotalSteps, hStep)
            
            
# This procedure works with the LAST VALUE of h used (which should be the same as the initial, 
#   given that the current parameter setting is for only one step in each of the x and h value
#   iterations).

                        
    print ' ' 
    print '  The original (at-equilibrium) unitArray is:'
    print ' '                                                 
    prettyPrintArray (arraySizeList, unitArray)                                                                                                
    print ' ' 
    
    perturbFrctn = 0.1  # amt of perturbation given as a fraction here; 
    #  if the program were revised to accept the perturbPrcnt from user, it should be taken in 
    #    as an integer and then converted to a fraction
    perturbedUnitArray = perturb (arraySizeList, unitArray, perturbFrctn)

    print ' ' 
    print '  With perturbation of ', perturbFrctn, 'The perturbedUnitArray is:'
    print ' '
    prettyPrintArray (arraySizeList, perturbedUnitArray)
    print ' '  
        
    print ' ' 
    print '  The total changes in the unitArray are (with W for 0=>1, and M for 1=>0):'
    print ' '
    prettyPrintChangesInUnitArray (arraySizeList, unitArray, perturbedUnitArray)               
    print ' '  
    
    # Store a copy of the starting perturbed array
    startingPerturbedArray = createIdenticalUnitArray (arraySizeList, perturbedUnitArray)    
    
    # Bring the perturbed array to equilibrium
    # Note that this will make changes to the starting perturbed array
    #  as well as returning an array called perturbedEqlbrmUnitArray;
    #  these will now be identical - 
    #  this is why we made a copy of the starting perturbed array. 
    perturbedEqlbrmUnitArray = adjustMatrixFEMinimum (arraySizeList, perturbedUnitArray, h, maxRange)

    
    print ' ' 
    print '  The previous perturbedUnitArray is:'
    print ' '
    prettyPrintArray (arraySizeList, startingPerturbedArray)
    print ' '      
    
    print ' ' 
    print '  After achieving equilibrium, the new perturbedUnitArray is:'
    print ' '
    prettyPrintArray (arraySizeList, perturbedEqlbrmUnitArray)
    print ' '  
    
    print ' ' 
    print '  The total changes in the new equilibrium perturbedUnitArray and the previous perturbedUnitArray are:'
    print ' '
    prettyPrintChangesInUnitArray (arraySizeList, startingPerturbedArray, perturbedEqlbrmUnitArray)               
    print ' ' 
    
                                                                                                        
    print ' ' 
    print '  The total changes in the new equilibrium perturbedUnitArray and the original unitArray are:'
    print ' '
    prettyPrintChangesInUnitArray (arraySizeList, unitArray, perturbedEqlbrmUnitArray)               
    print ' '                                                 
                                                                                                
####################################################################################################
# Conclude specification of the MAIN procedure
####################################################################################################                
    
if __name__ == "__main__": main()

####################################################################################################
# End program
####################################################################################################                
                                                                                                                                                                                                                                                                    