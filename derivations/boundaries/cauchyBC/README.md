# cauchyBC

**Warning**: The following methods has only proven first order, as seen in
[../../MES/boundaries/3-cauchyBC/conclusion.txt](../../MES/boundaries/3-cauchyBC/conclusion.txt)

Maybe thinking in the wrong way here?

[cauchyBCSubs0.ipynb](cauchyBCSubs0.ipynb) - Checkme
[cauchyBCSubs04thOrderDeriv.ipynb](cauchyBCSubs04thOrderDeriv.ipynb) - Checkme
[cauchyBCSubsF0.ipynb](cauchyBCSubsF0.ipynb) - Checkme
[cauchyBCTaylor4thOrder.ipynb](cauchyBCTaylor4thOrder.ipynb) - looks like this fellah is trying to derive Robin BCs
[cauchyHelper.mw](cauchyHelper.mw) -  Delme?
    120 - eq3 := 0=-alpha+(5*f0+15*f1-5*f2+f3)/16
    130 - eq4 := 0=-beta+(1/h)*(-f0+f1)
    138 - solve({eq3,eq4},{f0})
    ...
    146 - eq5 := -alpha+(5*f0+15*f1-5*f2+f3)/16 =-beta+(1/h)*(-f0+f1)
    154 - solve({eq5},{f0})


a = 0
b = 0

? => a = b?

yes...because 0=0

if so

=> a-b=0
but also
=> a+b=0
=> a*b=0

as a=b=0
