#!/usr/bin/env python

"""
Contains getCollectFieldsAndPlotOrder
"""

#{{{getCollectFieldsAndPlotOrder
def getCollectFieldsAndPlotOrder(fieldPlotType, hyperIncluded):
    #{{{docstring
    """
    Returns a tuple of the fields to collect and how to organize them

    Parameters
    ----------
    fieldPlotType : str
        What predefined fieldPlotType to plot.
    hyperIncluded : bool
        If hyper viscosities are used.

    Returns
    -------
    collectFields : tuple
        The fields to collect
    plotOrder : tuple
        The plot order
    """
    #}}}

    plotOrder = None
    if fieldPlotType == "mainFields":
        collectFields  = ("lnN"       ,\
                          "jPar"      ,\
                          "phi"       ,\
                          "vort"      ,\
                          "vortD"     ,\
                          "momDensPar",\
                          "S"         ,\
                         )
        plotOrder = ("lnN"  , "phi"       ,\
                     "n"    , "vortD"     ,\
                     "jPar" , "vort"      ,\
                     "uIPar", "momDensPar",\
                     "uEPar", "S"         ,\
                    )
    elif fieldPlotType == "mainFieldsBoussinesq":
        collectFields  = ("lnN"       ,\
                          "jPar"      ,\
                          "phi"       ,\
                          "vort"      ,\
                          "momDensPar",\
                          "S"         ,\
                         )
        plotOrder = ("lnN"  , "phi"       ,\
                     "n"    , "vort"      ,\
                     "jPar" , "momDensPar",\
                     "uIPar", "S"         ,\
                     "uEPar",\
                    )
    elif fieldPlotType == "lnN":
        collectFields  = ("ddt(lnN)"      ,\
                          "lnNAdv"        ,\
                          "gradUEPar"     ,\
                          "lnNUeAdv"      ,\
                          "srcN"          ,\
                          "lnNRes"        ,\
                          "lnNPerpArtVisc",\
                          "lnNParArtVisc" ,\
                         )
    elif fieldPlotType == "jPar":
        collectFields  = ("ddt(jPar)"      ,\
                          "jParAdv"        ,\
                          "uEParAdv"       ,\
                          "uIParAdv"       ,\
                          "jParParAdv"     ,\
                          "gradPhiLnN"     ,\
                          "jParRes"        ,\
                          "neutralERes"    ,\
                          "neutralIRes"    ,\
                          "jParPerpArtVisc",\
                          "jParParArtVisc" ,\
                         )
    elif fieldPlotType == "momDensPar":
        collectFields  = ("ddt(momDensPar)"   ,\
                          "momDensAdv"        ,\
                          "uIFluxAdv"         ,\
                          "uIParAdv"          ,\
                          "elPressure"        ,\
                          "densDiffusion"     ,\
                          "neutralIRes"       ,\
                          "neutralEResMu"     ,\
                          "momDensPerpArtVisc",\
                          "momDensParArtVisc" ,\
                         )
    elif fieldPlotType == "vortD":
        collectFields  = ["ddt(vortD)"                ,\
                          "vortNeutral"               ,\
                          "potNeutral"                ,\
                          "parDerDivUIParNGradPerpPhi",\
                          "vortDAdv"                  ,\
                          "kinEnAdvN"                 ,\
                          "divParCur"                 ,\
                          "vortDParArtVisc"           ,\
                          "vortDPerpArtVisc"          ,\
                         ]
        if hyperIncluded:
            collectFields.append("vortDHyperVisc")
        collectFields = tuple(collectFields)
    elif fieldPlotType == "vort":
        collectFields  = ["ddt(vort)"               ,\
                          "vortNeutral"             ,\
                          "potNeutral"              ,\
                          "DDYGradPerpPhiGradPerpUI",\
                          "vortAdv"                 ,\
                          "vortParAdv"              ,\
                          "divParCur"               ,\
                          "divSourcePhi"            ,\
                          "vortParArtVisc"          ,\
                          "vortPerpArtVisc"         ,\
                        ]
        if hyperIncluded:
            collectFields.append("vortHyperVisc")
        collectFields = tuple(collectFields)
    else:
        message = "{} is not implemented".format(fieldPlotType)
        raise NotImplementedError(message)

    return collectFields, plotOrder
#}}}
