# this is going to be a testing place from which I call teh various
# functions I've decided to make from the origional main.py.


# The way the instance cat works: I think, it is only run once for each catalogue, you make.
import numpy
from weakLensing import WL

myWL = WL()
print myWL.NbinsX
myWL.initialize()
print myWL.NbinsX

# let's make some numpy arrays!

ra = numpy.array([0.0, 0.1, 10.8, 12.2, 23.94, 23.94, 23.94])
dec = numpy.array([-1.73, 1.5, 0.2, -0.2, -1.2, -1.2, -1.1])
z = numpy.array([1.2, 0.8, 1.1, 2.1, 1.5, 1.7, 1.7])

shear1, shear2, conv = myWL.calc(ra, dec, z)

print shear1
