""" TrimGeneration

    ajc@astro June 30 2010
    Methods to generate trim file parameters - including metadata

"""

import numpy
import warnings

from lsst.sims.catalogs.measures.astrometry.Astrometry import *
from lsst.sims.catalogs.measures.photometry.Bandpass import *
from lsst.sims.catalogs.measures.photometry.Sed import *
#from lsst.sims.catalogs.measures.photometry.Magnitudes import *
from lsst.sims.catalogs.measures.instance.CatalogDescription import *
from lsst.sims.catalogs.measures.instance.InstanceCatalog import *
from lsst.sims.catalogs.measures.instance.Metadata import *


#load Slalib DLLs
slalib = ctypes.CDLL("slalsst.so")

def derivedTrimMetadata(instanceCatalog):
    '''Generate metadata needed for Trim files'''

    # calculate Slalib_date from Opsim_expmjd
    year=ctypes.c_int(0)
    month=ctypes.c_int(0)
    day=ctypes.c_int(0)
    fractionDay=ctypes.c_double(0.0)
    status=ctypes.c_int(0)

    filtMap = {'u':0, 'g':1, 'r':2, 'i':3, 'z':4, 'y':5}
    mjd = float(instanceCatalog.metadata.parameters['Opsim_expmjd'])
    #alt = float(instanceCatalog.metadata.parameters['Opsim_altitude'])
    #az = float(instanceCatalog.metadata.parameters['Opsim_azimuth'])
    filt = filtMap[instanceCatalog.metadata.parameters['Opsim_filter']]
    
    # slalib_date
    slalib.slaDjcl.argtypes = [ctypes.c_double,ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_int)]
    slalib.slaDjcl(mjd, year, month, day, fractionDay, status)

    slalibDate = "%4d/%02d/%02d/%.9g" % (year.value,month.value,day.value, fractionDay.value)
    instanceCatalog.metadata.addMetadata("Slalib_date", slalibDate, "Fractional date from Slalib")
    instanceCatalog.metadata.addMetadata("Opsim_filter", filt, "Remapped filter %s to the integer %i"%(instanceCatalog.metadata.parameters['Opsim_filter'], filt), clobber=True)

#Slalib_date  1994/10/12/0.0945639999991
#Slalib_expMjd 49637.094564

    #add new keywords for imsim configuration (3/1/2011)
    instanceCatalog.metadata.addMetadata("SIM_MINSOURCE", 1, "Minimum sources per chip for full_focalplane to consider simulating that chip", clobber=True)
    instanceCatalog.metadata.addMetadata("SIM_TELCONFIG", 0, "Configuration of the telescope.0 mean normal, 1 shutter is closed (i.e. for darks), 2 dome is illuminated", clobber=True)
    instanceCatalog.metadata.addMetadata("SIM_VISTIME", 32., "Effective dwell time of a visit", clobber=True)


if __name__ == '__main__':
    derivedTrimMetadata()
    
