# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 09:02:03 2018

@author: APO
"""

# Import libraries and whatnot
from scipy import pi
from scipy import signal
import matplotlib.pyplot as plt
import sympy as sy
from sympy import init_printing
from IPython.display import display

init_printing()
# ----------------------------------------------------------------------


def main():
    # CALL FUNCTIONS HERE
    LTI_1 = signal.lti([67471590909091],
                       [1, 7954545.4545455, 39062500000000])
    LTI_2 = signal.lti([269230769230770],
                       [1, 3076923.0769231, 100000000000000])
    LTI_3 = signal.lti([1.9375, 0, 0],
                       [1, 107323.23232323, 10203040506.071])
    # LTI_4 = signal.lti([1.725], [2.56e-14, 2.0368e-7, 1])

    pig = bodeFromTF(LTI_1, LTI_2, LTI_3)
    evalPoly(pig, printList=True)
    # getButter(2,994718.39)
# ----------------------------------------------------------------------


def evalPoly(polynomial, variable='s', printList=None):
    """
    Gets the coefficients of a polynom and returns them as a list. Uses
    one other function LTI_to_sympy provided in this source code.

    Needed imports:
        none
    args:
        polynomial ... input polynomial fraction, e.g. a transfer function
        variable ... variable in the polynomial, default 's'
        printList ... Wheter to print a the coeffcient list or not
    returns:
        2-d-list of numerator and denominator coefficients
    """

    # Convert our LTI polynom into a sympy representation
    polynomial = LTI_to_sympy(polynomial)
    # Split into numerator and denominator
    num = str(polynomial).split("/")[0]
    den = str(polynomial).split("/")[1]

    # Get rid of the parentheses at either end
    num = num.replace("(", "")
    num = num.replace(")", "")
    den = den.replace("(", "")
    den = den.replace(")", "")

    # Split lists into separate terms
    num = num.split(" + ")
    den = den.split(" + ")

    # Remove whitespaces at either end of the coefficients
    for item in num:
        item = item.strip()
    for item in den:
        item = item.strip()

    # Generate empty lists
    coeffNum = []
    coeffDen = []

    # If there is term with a variable ...
    if not num[0].find("*") == -1:
        # ... split the term to get the power, remove ** and that's our order
        lim = int(num[0].split("*" + variable)[1].replace("**", ""))
    # Else set the iteration limit to 0
    else:
        lim = 0

    # Iterate once more than order (2nd order polynom means 3 iterations)
    for iterator in range(lim+1):
        # Split item into coefficient and power
        try:
            item = num[iterator].split("*" + variable)

            # If there is a power term ...
            if len(item) > 1:
                # ... remove ** from current item power
                curItemPower = item[1].replace("**", "")
                # ... get current power that we should be at:
                # (limit - already performed iterations)
                curIterPower = str(lim-len(coeffNum))

                # If entries match or the coefficient is not empty ...
                if curItemPower == curIterPower or not item[0] == '':
                    # ... add coefficient to list
                    coeffNum.append(item[0])
                # else add 0 to the list
                else:
                    coeffNum.append("0")
            # else add the coefficient right away (last coefficient)
            else:
                if not item[0] == '':
                    coeffNum.append(item[0])
        except:
            coeffNum.append("0")

    # If there is term with a variable ...
    if not den[0].find("*") == -1:
        # Else set the iteration limit to 0
        lim = int(den[0].split("*" + variable)[1].replace("**", ""))
    else:
        lim = 0

    # Iterate once more than order (2nd order polynom means 3 iterations)
    for iterator in range(lim+1):
        try:
            # Split item into coefficient and power
            item = den[iterator].split("*" + variable)

            # If there is a power term ...
            if len(item) > 1:
                # ... remove ** from current item power
                curItemPower = item[1].replace("**", "")
                # ... get current power that we should be at:
                # (limit - already performed iterations)
                curIterPower = str(lim-len(coeffDen))

                # If entries match or the coefficient is not empty ...
                if curItemPower == curIterPower or not item[0] == '':
                    # ... add coefficient to list
                    coeffDen.append(item[0])
                # else add 0 to the list
                else:
                    coeffDen.append("0")
            # else add the coefficient right away (last coefficient)
            else:
                if not item[0] == '':
                    coeffDen.append(item[0])
        except:
            coeffNum.append("0")

    # Create empty list to return
    coefficients = []
    # Add numerator and denominator to list (as lists)
    coefficients.append(coeffNum)
    coefficients.append(coeffDen)

    # If user wants to have list printed
    if printList:
        for list in coefficients:
            print(list)

    # return list
    return coefficients
# ----------------------------------------------------------------------


def LTI_to_sympy(LTI, symplify=True):
    """
    Convert LTI instance to sympy expression

    Needed imports:
        import sympy as sy
    args:
        LTI ... transfer function in linear time invariant form (scipy)
        symplify ... wheter to return a simplified equation or not
    returns:
        a sympy representation of the provided LTI transfer function
    """

    # Set symbol to use as variable
    s = sy.Symbol('s')
    # Convert transfer function
    G = sy.Poly(LTI.num, s) / sy.Poly(LTI.den, s)
    # Either return a simplified version or the original
    return sy.simplify(G) if symplify else G
# ----------------------------------------------------------------------


def sympy_to_LTI(xpr, s=sy.Symbol('s')):
    """
    Convert sympy transfer function polynomial to LTI

    Needed imports:
        import sympy as sy
        from scipy import signal
    args:
        xpr ... the expression to convert from sympy to LTI
        s ... symbol to use as variable after conversion (default s)
    returns:
        a LTI representation of the provided sympy transfer function
    """

    # Get numerator and denominator for LTI tf
    num, den = sy.simplify(xpr).as_numer_denom()
    # Generate polynomials from num and den
    pNumDen = sy.poly(num, s), sy.poly(den, s)
    # Expand polynomials to get nice sum representation
    cNumDen = [sy.expand(p).all_coeffs() for p in pNumDen]
    # Convert coefficients to floats
    lNum, lDen = [sy.lambdify((), c)() for c in cNumDen]
    # Return LTI representation
    return signal.lti(lNum, lDen)
# ----------------------------------------------------------------------


def bodeFromTF(*args):
    """
    Combines all provided transfer functions (in LTI form), plots them in
    a bode diagram, gives the cutoff frequency and the transfer function.
    Uses three other functions LTI_to_sympy, sympy_to_LTI and findCutoff
    provided in this source file.

    Needed imports:
        from scipy import signal
        import matplotlib.pyplot as plt
        from IPython.display import display
    args:
        user can supply as many LTI transfer functions as he likes
    returns:
        returns a the combined transfer function in LTI form
    """

    # Create combined variable and initialize it with the first tf
    combinedTF = LTI_to_sympy(args[0])

    # Loop through all TFs and multiply with the combined variable
    for LTI_TF in args[1:]:
        combinedTF = combinedTF * LTI_to_sympy(LTI_TF)

    # Expand to full sum representation (normalize so max power has coeff 1)
    combinedTF = sy.simplify(combinedTF).expand()

    # Set up symbols and show transfer function as nice fraction
    pH = sy.symbols("H")
    display(sy.Eq(pH, combinedTF))

    # Convert transfer function back to LTI representation
    combinedTF = sympy_to_LTI(combinedTF)
    # Get data for bode plot
    w, mag, phase = signal.bode(combinedTF)
    # Change to Hz instead of rad/s
    w = w/(2*pi)

    # Find cutoff freqeuncy and print it
    cutHigh, cutLow = findCutoff(list(w), list(mag))
    print("Lower cutoff frequency: ", cutLow)
    print("Upper cutoff frequency: ", cutHigh)

    # Set up plot
    plt.figure(figsize=(10, 8))  # Create new figure
    plt.subplot(211)  # Create new subplot
    plt.semilogx(w, mag)  # Bode magnitude plot
    plt.axvline(cutLow, color='green')  # cutoff frequency plot
    plt.axvline(cutHigh, color='green')  # cutoff frequency plot
    plt.margins(0, 0.1)  # Set margins
    plt.title('Filter frequency response')  # Set figure title
    plt.ylabel('Amplitude [dB]')  # Set label for y axis
    plt.grid(which='both', axis='both')  # Activate grid for both axis
    plt.subplot(212)  # Create new subplot
    plt.semilogx(w, phase)  # Bode phase plot
    plt.margins(0, 0.1)  # Set margins
    plt.xlabel('Frequency [degrees]')  # Set label for x axis
    plt.ylabel('Phase [°]')  # Set lavel for y axis
    plt.grid(which='both', axis='both')  # Activate grid for both axis
    plt.show()  # Show plot

    # Return combined transfer function in LTI form
    return combinedTF
# ----------------------------------------------------------------------


def getButter(order, cutoff, filterType=None):
    """
    Generates an analog butterworth filter of provided order with provided
    cutoff frequency. Does not return anything but calls another function
    bodeFramTF I have here that processes the filter further.

    Needed imports:
        from scipy import signal
    args:
        order ... order of the filter
        cutoff ... cutoff frequency in Hz
        filterType ... 'low' or 'high'
    returns:
        none
    """

    # Default filter type is low
    if filterType is None:
        filterType = 'low'

    # Get filter numerator and denominator
    num, den = signal.butter(order, cutoff, filterType, analog=True)

    # Convert to linear time invariant
    LTI = signal.lti(num, den)

    # Show plots 'n shit
    bodeFromTF(True, LTI)
# ----------------------------------------------------------------------


def findCutoff(w, mag):
    """
    Finds the cutoff frequency by looping until the difference betweeen a
    moving average and the current value is greater 3 dB, then calculating
    slope between found point and the point before and finally
    interpolating to find the exact frequency for 3 dB

    Needed imports:
        none
    args:
        NOTE: w and mag must have compatible indices and sort order
        w ... collection of frequencies for the transfer function
        mag ... collection of magnitudes for the transfer function (in dB)
    returns:
        cutoff frequency in Hz
    """

    # Find index of max value
    maxIndex = mag.index(max(mag))

    # Create dummy variables
    sum = 0
    movingAvg = 0
    items = 0

    # Loop through collection of magnitudes
    for i in range(maxIndex, len(mag)):
        value = mag[i]

        # If difference of average and value greater 3 dB ...
        if (movingAvg - value) > 3:
            # ... keep track of the index
            index = mag.index(value)
            # Calculate slope using current point and the point before
            slope = (mag[index] - mag[index-1]) / (w[index] - w[index-1])
            # Calculate constant term
            c = mag[index] - w[index] * slope
            # Calculate cutoff frequency from linear equation
            cutoffUpper = (mag[index-1] - 3 - c) / slope
            # Exit loop
            break

        # Calculate sum of already passed magnitudes
        sum = sum + value
        items += 1
        # Calculate average value of already passed magnitudes
        movingAvg = sum / items

    # Reset dummy variables
    sum = 0
    movingAvg = 0
    items = 0

    # Loop through collection of magnitudes
    for i in range(maxIndex, 0, -1):
        value = mag[i]

        # If difference of average and value greater 3 dB ...
        if (movingAvg - value) > 3:
            # ... keep track of the index
            index = mag.index(value)
            # Calculate slope using current point and the point before
            slope = (mag[index] - mag[index+1]) / (w[index] - w[index+1])
            # Calculate constant term
            c = mag[index] - w[index] * slope
            # Calculate cutoff frequency from linear equation
            cutoffLower = (mag[index+1] - 3 - c) / slope
            # Exit loop
            break

        # Calculate sum of already passed magnitudes
        sum = sum + value
        items += 1
        # Calculate average value of already passed magnitudes
        movingAvg = sum / items

    # Return cutoff frequency
    return cutoffUpper, cutoffLower
# ----------------------------------------------------------------------


main()
