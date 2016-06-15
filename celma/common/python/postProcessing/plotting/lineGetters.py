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
    mainFields.lines.append(Line('lnN'  , r'\ln(n)'         , plotPos=0))
    mainFields.lines.append(Line('uIPar', r'u_{i,\parallel}', plotPos=6))
    mainFields.lines.append(Line('uEPar', r'u_{e,\parallel}', plotPos=4))
    mainFields.lines.append(Line('vortD', r'\Omega^D'       , plotPos=3))
    # Helping field
    mainFields.lines.append(Line('phi'  , r'\phi'           , plotPos=1))
    mainFields.lines.append(Line('S'    , r'S'              , plotPos=7))
    mainFields.lines.append(Line('vort' , r'\Omega'         , plotPos=5))
    # Extra lines
    # FIXME: This API is not so intuitive, consider change
    mainFields.extraLines['jPar'] = Line('jPar', r'j_\parallel', plotPos=8)
    mainFields.extraLines['n'   ] = Line('n'   , r'n'          , plotPos=2)

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
        r'-\frac{1}{J}\{\phi,\ln(n)\}'                         ))
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

#{{{getUEParFields
def getUEParFields(path):
    """
    Prepares the uEPar fields
    """

    # Making the orgObj instance
    uEPar = Organizer(r"u_{e,\parallel}", useCombinedPlot=True, path = path)
    # Making lines in the pattern name, lable, plotPos
    uEPar.lines.append(Line('uEParAdv'                                     ,\
    r'-\frac{1}{J}\{\phi,u_{e,\parallel}\}'                                ))
    uEPar.lines.append(Line('uEParParAdv'                                  ,\
    r'- u_{e,\parallel}\partial_{\parallel}u_{e,\parallel}'                ))
    uEPar.lines.append(Line('gradPhiLnN'                                   ,\
    r'\mu\partial_\parallel\left( \phi - \ln(n)\right)'                    ))
    uEPar.lines.append(Line('uEParRes'                                     ,\
    r'-0.51\nu_{ei}\left(u_{e,\parallel}-u_{i,\parallel} \right)'          ))
    uEPar.lines.append(Line('ueSrc'                                        ,\
    r'-\frac{S u_{e,\parallel}}{n}'                                        ))
    uEPar.lines.append(Line('ueNeutral'                                    ,\
    r'-\nu_{en}u_{e,\parallel}'                                            ))
    uEPar.lines.append(Line('uEParPerpArtVisc'                             ,\
    r'D_{u_{e,\parallel}, \perp}'                                          +\
    r'\frac{\partial^2_{\parallel}u_{e,\parallel}}{n}'                     ))
    uEPar.lines.append(Line('uEParParArtVisc'                              ,\
    r'D_{u_{e,\parallel}, \parallel}'                                      +\
    r'\frac{\partial^2_{\parallel}u_{e,\parallel}}{n}'                     ))
    uEPar.lines.append(Line('uEParPerpArtViscNoN'                          ,\
    r'D_{u_{e,\parallel}, \perp}\nabla^2_\perp u_{e,\parallel}'            ))
    uEPar.lines.append(Line('uEParParArtViscNoN'                           ,\
    r'D_{u_{e,\parallel}, \parallel} \partial^2_{\parallel}u_{e,\parallel}'))

    return uEPar
#}}}

#{{{getUIParFields
def getUIParFields(path):
    """
    Prepares the uIPar fields
    """

    # Making the orgObj instance
    uIPar = Organizer(r"u_{i,\parallel}", useCombinedPlot=True, path = path)
    # Making lines in the pattern name, lable, plotPos
    uIPar.lines.append(Line('uIParAdv'                                    ,\
             r'-\frac{1}{J}\{\phi,u_{i,\parallel}\}'                      ))
    uIPar.lines.append(Line('uIParParAdv'                                 ,\
             r'-u_{i,\parallel}\partial_{\parallel}u_{i,\parallel}'       ))
    uIPar.lines.append(Line('gradPhi'                                     ,\
             r'-\partial_\parallel\phi'                                   ))
    uIPar.lines.append(Line('uIParRes'                                    ,\
             r'-0.51\nu_{ei}\left(u_{i,\parallel}-u_{e,\parallel} \right)'))
    uIPar.lines.append(Line('uiNeutral'                                   ,\
             r'-\nu_{in}u_{i,\parallel}'                                  ))
    uIPar.lines.append(Line('uiSrc'                                       ,\
             r'-\frac{S u_{i,\parallel}}{n}'                              ))
    uIPar.lines.append(Line('uIParPerpArtVisc'                            ,\
             r'D_{u_i,\perp}'                                             +\
             r'\frac{\partial^2_{\parallel}u_{i,\parallel}}{n}'           ))
    uIPar.lines.append(Line('uIParParArtVisc'                             ,\
             r'D_{u_i,\parallel}'                                         +\
             r'\frac{\partial^2_{\parallel}u_{i,\parallel}}{n}'           ))
    uIPar.lines.append(Line('uIParPerpArtViscNoN'                         ,\
             r'D_{u_{i,\parallel}, \perp}\nabla^2_\perp u_{i,\parallel}'  ))
    uIPar.lines.append(Line('uIParParArtViscNoN'                          ,\
             r'D_{u_i,\parallel} \partial^2_{\parallel}u_{i,\parallel}'   ))

    return uIPar
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
    vortD.lines.append(Line('divExBAdvGradPerpPhiN'                ,\
         r'-\nabla\cdot\left('                                     +\
         r'\mathbf{u}_E \cdot\nabla'                               +\
         r'\left[ n \nabla_\perp \phi \right]'                     +\
         r'\right)'                                                ))
    vortD.lines.append(Line('vortDAdv'                             ,\
         r'-\frac{1}{J}\{\phi, \Omega^D\}'                         ))
    vortD.lines.append(Line('kinEnAdvN'                            ,\
         r'-\frac{1}{2J}\{\mathbf{u}_E \cdot\mathbf{u}_E, n\}'     ))
    vortD.lines.append(Line('parDerDivUIParNGradPerpPhi'           ,\
         r'-\partial_\parallel \nabla\cdot\left('                  +\
         r'u_{i,\parallel} '                                       +\
         r'n \nabla_\perp \phi'                                    +\
         r'\right)'                                                ))
    vortD.lines.append(Line('nGradUiUe'                            ,\
         r'n\partial_\parallel (u_{i,\parallel}-u_{e,\parallel})'  ))
    vortD.lines.append(Line('uiUeGradN'                            ,\
        r'(u_{i,\parallel}-u_{e,\parallel})\partial_{\parallel}n'  ))
    vortD.lines.append(Line('divParCur'                            ,\
        r'\partial_{\parallel}(n[u_{i,\parallel}-u_{e,\parallel}])'))
    vortD.lines.append(Line('vortDParArtVisc'                      ,\
         r'D_{\Omega^D} \partial^2_{\parallel}\Omega^D'            ))
    vortD.lines.append(Line('vortDPerpArtVisc'                     ,\
         r'D_{\Omega^D, \perp} \nabla_\perp^2\Omega^D'             ))

    return vortD
#}}}
#}}}
