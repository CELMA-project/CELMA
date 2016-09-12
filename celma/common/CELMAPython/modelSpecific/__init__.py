#!/usr/bin/env python

""" Init for the modelSpecific package """

from .celmaModel import getOrgObjFromModel

labelNames = ["mainFields" ,\
              "lnNFields"  ,\
              "uEParFields",\
              "uIParFields",\
              "vortDFields",\
              ]

varAndPlotNames = [
                   ("lnN", r"\ln(n)"),\
                   ("phi", r"\phi"),\
                   ("n", r"n"),\
                   ("vort", r"\Omega"),\
                   ("vortD", r"\Omega^D"),\
                   ("uIPar", r"u_{i,\parallel}"),\
                   ("uEPar", r"u_{e,\parallel}"),\
                  ]
