#!/usr/bin/env python

""" Init for the modelSpecific package """

from .celmaCurMomModel import getOrgObjFromModel

labelNames = ("mainFields",\
              "lnNFields",\
              "jMParFields",\
              "momDensParFields",\
              "vortDFields",\
              )

varAndPlotNames = (
                   ('lnN', r'\ln(n)'),\
                   ("jMPar", r'j^M_\parallel'),\
                   ('phi', r'\phi'),\
                   ("n", r'n'),\
                   ('vort', r'\Omega'),\
                   ('vortD', r'\Omega^D'),\
                   ('uIPar', r'u_{i,\parallel}'),\
                   ('uEPar', r'u_{e,\parallel}'),\
                  )
