#!/usr/bin/env python

from boutdata import collect
import numpy as np

def getVarDict(path, var=''):
    """
    Returns a dictionary of grouped variables chosen by var
    """

    if var == 'ddt':
        theVars = [
                  'ddt(lnN)',\
                  'ddt(jPar)',\
                  'ddt(momDensPar)',\
                  'ddt(vortD)',\
                  ]
    elif var == 'lnN':
        theVars = [
                  'lnNAdv',\
                  'lnNRes',\
                  'gradUEPar',\
                  'lnNUeAdv',\
                  'srcN',\
                  'lnNParArtVisc',\
                  'lnNPerpArtVisc',\
                  ]
    elif var == 'jPar':
        theVars = [
                  'jParAdv',\
                  'uIParAdvSum',\
                  'uEParDoubleAdv',\
                  'jParRes',\
                  'elField',\
                  'muElPressure',\
                  'neutralERes',\
                  'neutralIRes',\
                  'jParParArtVisc',\
                  'jParPerpArtVis',\
                   ]
    elif var == 'momDensPar':
        theVars = [
                  'momDensAdv',\
                  'elPressure',\
                  'densDiffusion',\
                  'momDensParArtVisc',\
                  'momDensPerpArtVisc',\
                  ]
    elif var == 'vortD':
        theVars = [
                  'vortNeutral',\
                  'potNeutral',\
                  'parDerDivUIParNGradPerpPhi',\
                  'vortDAdv,',\
                  'kinEnAdvN',\
                  'divParCur',\
                  'vortDParArtVisc,',\
                  'vortDPerpArtVisc',\
                  'vortDhyperVisc',\
                  ]
    else:
        raise ValueError("Input var not recognized")

    varDict = {}
    for var in theVars:
        varDict[var] = collect(var, xguards=False, yguards=False)

    return varDict

def printVarDict(varDict):
    """
    Print the MaxAbs, mean and variance for the group of variables given
    as varDict.

    Useful for finding stiff variables when using IMEX
    """

    for variable in varDict.keys():
        print("Var: {:<20} MaxAbs: {:<20.3e} Mean: {:<20.3e} Var: {:<20.3e} ".\
                format(variable,\
                       np.max(np.abs(varDict[variable])),\
                       np.mean(varDict[variable]),\
                       np.var(varDict[variable]),\
                        )
              )
