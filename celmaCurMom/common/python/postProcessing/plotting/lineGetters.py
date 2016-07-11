#!/usr/bin/env python

"""
Contains drivers for plotting 1D plots
"""

from .organizer import Organizer
from .line import Line

#{{{Field getters
#{{{getMainFields
def getMainFields(path):
    """
    Prepares the main fields
    """

    # Making the orgObj instance
    mainFields = Organizer("mainFields", path = path)
    # Making lines in the pattern name, lable, plotPos
    # Evolved fields
    mainFields.lines.append(Line('lnN'        , r'\ln(n)'          , plotPos=0))
    mainFields.lines.append(Line('jPar'       , r'j_{\parallel}'   , plotPos=8))
    mainFields.lines.append(Line('momDensPar' , r'nu_{i,\parallel}', plotPos=1))
    mainFields.lines.append(Line('vortD'      , r'\Omega^D'        , plotPos=3))
    # Helping field
    mainFields.lines.append(Line('uIPar', r'u_{i,\parallel}', plotPos=6))
    mainFields.lines.append(Line('uEPar', r'u_{e,\parallel}', plotPos=4))
    mainFields.lines.append(Line('phi'  , r'\phi'           , plotPos=7))
    mainFields.lines.append(Line('S'    , r'S'              , plotPos=9))
    mainFields.lines.append(Line('vort' , r'\Omega'         , plotPos=5))
    # Extra lines
    # FIXME: This API is not so intuitive, consider change
    mainFields.extraLines['n'] = Line('n'   , r'n'          , plotPos=2)

    return mainFields
#}}}

#{{{getLnNFields
def getLnNFields(path):
    """
    Prepares the lnN fields
    """

    # Making the orgObj instance
    lnN = Organizer(r"\ln(n)", useCombinedPlot=True, path = path)
    # Making lines in the pattern name, lable, plotPos
    # Evolved fields
    lnN.lines.append(Line('lnNAdv'                             ,\
        r'-\frac{1}{JB}\{\phi,\ln(n)\}'                        ))
    lnN.lines.append(Line('lnNRes'                             ,\
        r'\frac{0.51\nu_{ei}}{\mu}\left(\nabla^2_\perp\ln(n) +'+\
        r'\left[\nabla_\perp\ln(n)\right]^2\right)'            ))
    lnN.lines.append(Line('gradUEPar'                          ,\
        r'-\partial_{\parallel}u_{e,\parallel}'                ))
    lnN.lines.append(Line('lnNUeAdv'                           ,\
        r'-u_{e,\parallel}\partial_\parallel\ln(n)'            ))
    lnN.lines.append(Line('srcN'                               ,\
        r'\frac{S}{n}'                                         ))
    lnN.lines.append(Line('lnNPerpArtVisc'                     ,\
        r'D_{n,\perp} \partial^2_{\perp}\ln(n)'                ))
    lnN.lines.append(Line('lnNParArtVisc'                      ,\
        r'D_{n,\parallel} \partial^2_{\parallel}\ln(n)'        ))

    return lnN
#}}}

#{{{getJParFields
def getJParFields(path):
    """
    Prepares the jPar fields
    """

    # Making the orgObj instance
    jPar = Organizer(r"j_{\parallel}", useCombinedPlot=True, path = path)
    # Making lines in the pattern name, lable, plotPos
    jPar.lines.append(Line('jParAdv'                                    ,\
    r'-\frac{1}{JB}\{\phi,j_{\parallel}\}'                              ))
    jPar.lines.append(Line('uIParAdvSum'                                ,\
    r'- u_{e,\parallel}\partial_{\parallel}'                            +\
    r'\left(n\left[u_{i,\parallel}+u_{e,\parallel}\right]\right)'       ))
    jPar.lines.append(Line('uEParDoubleAdv'                             ,\
    r'2u_{e,\parallel}\partial_{\parallel}\left(nu_{e,\parallel}\right)'))
    jPar.lines.append(Line('jParRes'                                    ,\
    r'-0.51\nu_{ei}j_\parallel'                                         ))
    jPar.lines.append(Line('elPressure'                                 ,\
    r'\mu T_e\partial_{\parallel}n'                                     ))
    jPar.lines.append(Line('elField'                                    ,\
    r'-\mu n \partial_{\parallel}\phi'                                  ))
    jPar.lines.append(Line('neutralERes'                                ,\
    r'n\nu_{en}u_{e,\parallel}'                                         ))
    jPar.lines.append(Line('neutralIRes'                                ,\
    r'-n\nu_{in}u_{i,\parallel}'                                        ))
    jPar.lines.append(Line('jParPerpArtVisc'                            ,\
    r'D_{j_{\parallel}, \perp}\nabla^2_\perp j_{\parallel}'             ))
    jPar.lines.append(Line('jParParArtVisc'                             ,\
    r'D_{j_{\parallel}, \parallel} \partial^2_{\parallel}j_{\parallel}' ))

    return jPar
#}}}

#{{{momDensParFields
def momDensParFields(path):
    """
    Prepares the n*momDensPar fields
    """

    # Making the orgObj instance
    momDensPar = Organizer(r"(nu_{i,\parallel})", useCombinedPlot=True, path = path)
    # Making lines in the pattern name, lable, plotPos
    momDensPar.lines.append(Line('momDensAdv'                      ,\
    r'-\frac{1}{JB}\{\phi, nu_{i,\parallel}\}'                     ))
    momDensPar.lines.append(Line('uIParAdvSum'                     ,\
    r'- u_{e,\parallel}\partial_{\parallel}'                       +\
    r'\left(n\left[u_{i,\parallel}+u_{e,\parallel}\right]\right)'  ))
    momDensPar.lines.append(Line('elPressure'                      ,\
    r'\mu T_e\partial_{\parallel}n'                                ))
    momDensPar.lines.append(Line('neutralIRes'                     ,\
    r'-n\nu_{in}u_{i,\parallel}'                                   ))
    momDensPar.lines.append(Line('densDiffusion'                   ,\
    r'0.51\frac{\nu_{ei}}{\mu}u_{i,\parallel}\nabla_\perp^2n'      ))
    momDensPar.lines.append(Line('momDensPerpArtVisc'              ,\
             r'D_{nu_i,\perp}'                                     +\
             r'\nabla^2_{\perp}\left(nu_{i,\parallel}\right)'      ))
    momDensPar.lines.append(Line('momDensParArtVisc'               ,\
             r'D_{nu_i,\parallel}'                                 +\
             r'\partial^2_{\parallel}nu_{i,\parallel}'             ))

    return momDensPar
#}}}

#{{{getVortDFields
def getVortDFields(path):
    """
    Prepares the vortD fields
    """

    # Making the orgObj instance
    vortD = Organizer(r"\Omega^D", useCombinedPlot=True, path = path)
    # Making lines in the pattern name, lable, plotPos
    vortD.lines.append(Line('vortNeutral'                          ,\
         r'-\nu_{in}n\Omega'                                       ))
    vortD.lines.append(Line('potNeutral'                           ,\
         r'-\nu_{in}\nabla_\perp \phi \cdot \nabla_\perp n'        ))
    vortD.lines.append(Line('vortDAdv'                             ,\
         r'-\frac{1}{JB}\{\phi, \Omega^D\}'                        ))
    vortD.lines.append(Line('kinEnAdvN'                            ,\
         r'-\frac{1}{2J}\{\mathbf{u}_E \cdot\mathbf{u}_E, n\}'     ))
    vortD.lines.append(Line('parDerDivUIParNGradPerpPhi'           ,\
         r'-\partial_\parallel \nabla\cdot\left('                  +\
         r'u_{i,\parallel} '                                       +\
         r'n \nabla_\perp \phi'                                    +\
         r'\right)'                                                ))
    vortD.lines.append(Line('divParCur'                            ,\
        r'\partial_{\parallel}j_{\parallel}'                       ))
    vortD.lines.append(Line('vortDParArtVisc'                      ,\
         r'D_{\Omega^D} \partial^2_{\parallel}\Omega^D'            ))
    vortD.lines.append(Line('vortDPerpArtVisc'                     ,\
         r'D_{\Omega^D, \perp} \nabla_\perp^2\Omega^D'             ))

    return vortD
#}}}
#}}}
