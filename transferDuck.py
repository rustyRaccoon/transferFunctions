# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 09:02:03 2018

@author: APO
"""

from scipy import signal
import matplotlib.pyplot as plt
import sympy as sy
from sympy import init_printing
from IPython.display import display

init_printing()

#----------------------------------------------------------------------
def main():
    # INSERT YOUR TRANSFER FUNCTIONS HERE
    LTI_1 = signal.lti([67471590909091],[1,7954545.4545455,39062500000000])
    LTI_2 = signal.lti([269230769230770],[1,3076923.0769231,100000000000000])
    LTI_3 = signal.lti([1.9375,0,0],[1,107323.23232323,10203040506.071])

    # CALL FUNCTIONS HERE
    bodeFromTF(True, LTI_1)
    bodeFromTF(True, LTI_2)
    bodeFromTF(True, LTI_3)
    pig = bodeFromTF(True, LTI_1, LTI_2, LTI_3)
    
    # GETS POLYNOMIAL COEFFICIENTS
    evalPoly(pig,printList = True)

#----------------------------------------------------------------------
def evalPoly(polynomial, variable='s', printList = None):
    
    polynomial = LTI_to_sympy(polynomial)
    
    num = str(polynomial).split("/")[0]
    den = str(polynomial).split("/")[1]
    
    num = num.replace("(","")
    num = num.replace(")","")
    den = den.replace("(","")
    den = den.replace(")","")    
    
    num = num.split(" + ")
    den = den.split(" + ")
    
    for item in num:
        item = item.strip()
    for item in den:
        item = item.strip()
    
    coeffNum = []
    coeffDen = []
    
    try:
        limit = int(num[0].split("*" + variable)[1].replace("**",""))
    except:
        limit = 0
        
    for iterator in range(limit+1):
        try:
            item = num[iterator].split("*" + variable)
            try:
                
                if item[1].replace("**","") == str(limit-len(coeffNum)) or not item[0] == '':
                    coeffNum.append(item[0])
                else:
                    coeffNum.append("0")
            except:
                if not item[0] == '':
                    coeffNum.append(item[0])
        except:
            coeffNum.append("0")
    
    try:
        limit = int(den[0].split("*" + variable)[1].replace("**",""))
    except:
        limit = 0
        
    for iterator in range(limit+1):
        try:
            item = den[iterator].split("*" + variable)
            try:
                
                if item[1].replace("**","") == str(limit-len(coeffDen)) or not item[0] == '':
                    coeffDen.append(item[0])
                else:
                    coeffDen.append("0")
            except:
                if not item[0] == '':
                    coeffDen.append(item[0]) 
        except:
            coeffDen.append("0")
            
    coefficients = []
    coefficients.append(coeffNum)
    coefficients.append(coeffDen)
    
    if printList == True:
        for list in coefficients:
            print(list)
        
    return coefficients

#----------------------------------------------------------------------
def LTI_to_sympy(lsys, symplify=True):
    """ Convert Scipy's LTI instance to Sympy expression """
    s = sy.Symbol('s')
    G = sy.Poly(lsys.num, s) / sy.Poly(lsys.den, s)
    return sy.simplify(G) if symplify else G

#----------------------------------------------------------------------
def sympy_to_LTI(xpr, s=sy.Symbol('s')):
    """ Convert Sympy transfer function polynomial to Scipy LTI """
    num, den = sy.simplify(xpr).as_numer_denom()  # expressions
    p_num_den = sy.poly(num, s), sy.poly(den, s)  # polynomials
    c_num_den = [sy.expand(p).all_coeffs() for p in p_num_den]  # coefficients
    l_num, l_den = [sy.lambdify((), c)() for c in c_num_den]  # convert to floats
    return signal.lti(l_num, l_den)

#----------------------------------------------------------------------
def bodeFromTF(showTF, *args):
    combinedTF = LTI_to_sympy(args[0])
    
    for LTI_TF in args[1:]:
        combinedTF = combinedTF * LTI_to_sympy(LTI_TF)
      
    combinedTF = sy.simplify(combinedTF).expand()
    
    if showTF:
        pH = sy.symbols("H")
        display(sy.Eq(pH, combinedTF))
    
    combinedTF = sympy_to_LTI(combinedTF)
    w, mag, phase = signal.bode(combinedTF)
    
    plt.figure(figsize=(10,8))
    plt.subplot(211)
    plt.semilogx(w, mag) # Bode magnitude plot
    plt.subplot(212)
    plt.semilogx(w, phase) # Bode phase plot
    plt.show()
    
    return combinedTF

#----------------------------------------------------------------------
main()