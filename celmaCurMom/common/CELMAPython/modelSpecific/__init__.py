#!/usr/bin/env python

""" Init for the modelSpecific package """

from .celmaCurMomModel import getOrgObjFromModel

labelNames = ["mainFields",\
              "lnNFields",\
              "jParFields",\
              "momDensParFields",\
              "vortDFields",\
              "vortFields",\
              ]

varAndPlotNames = [
                   ("lnN", r"\ln(n)"),\
                   ("jPar", r"j_\parallel"),\
                   ("phi", r"\phi"),\
                   ("n", r"n"),\
                   ("vort", r"\Omega"),\
                   ("vortD", r"\Omega^D"),\
                   ("uIPar", r"u_{i,\parallel}"),\
                   ("uEPar", r"u_{e,\parallel}"),\
                  ]
