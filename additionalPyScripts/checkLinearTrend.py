import numpy as np
import matplotlib.pylab as plt
B=0.1
# Values from param dict in analytic growth rates
sigVal =6051.722409/(0.1**2)
bVal   =1.260301*B**2
omsVal =-105172*B

# Dependencies
sig = lambda B: sigVal*B**2
b = lambda B: bVal/(B**2)
oms = lambda B: omsVal/B

# Real expression
def real(B):
    return\
    -(np.sqrt(sig(B))/2)*\
    (16*oms(B)**2 + sig(B)**2*(b(B)**2+2*b(B)+1)**2)**(1/4)*\
    np.cos(0.5*(np.arctan2(4*oms(B),-sig(B)*(b(B)**2+2*b(B)+1))))

# Plot
x = np.linspace(0.02,0.1,100)
y = real(x)
plt.plot(x,y)
plt.show()
