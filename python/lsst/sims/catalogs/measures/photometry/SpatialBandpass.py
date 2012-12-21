""" 
SpatialBandpass - 

ljones@astro.washington.edu

This class calculates a spatially varying bandpass as needed for
running 'sim1' (aka generating a calibration simulation catalog). As such, there are specific
spatially varying terms we're using:
 variation in the temperature of the detector (together with a temperature-dependent QE curve),
  input data - relative QE(T) as ascii file, plus QE(x,y) over CCD as ascii file (base filename),
   plus Temperature offset values per CCD
  interpolate over a CCD, but no further, smooth wavelength response
 variation in the filter throughput curves
  input data format?? input data - transmission curve (lambda) at various mmx/mmy location in focal plane
  filter = fixed with respect to x/y in focal plane
  interpolate over entire focal plane, relatively smooth wavelength dependence
 variation in the transmission curve for telescope/optics
  input data format?? input data - transmission curve (lambda) at various (mm'x/mm'y) locations in sky plane (need rotation
  angle to translate to focal plane)
  interpolate over entire focal plane, relatively smooth wavelength dependence
 variation in the grey-scale flat field - part 1
  input data format? (probably MEF fits file, one per CCD).
  interpolate over a CCD, but no further. No wavelength response. Fixed to focal plane. 
 variation in the grey-scale flat field - part 2a
  input data format? (probably single fits file for entire fov). part of grey flat field fixed to camera. Relates to camera.
  interpolate over entire fov (mmx/mmy)
 variation in the grey-scale flat field - part 2b
  input data format? (probably single fits file for entire fov). part of grey flat field fixed to telescope, must
   translate to focal plane using rotator angle). Relates to spiders, other telescope pieces.
  interpolate over entire fov (mmx'/mmy')
All of the above are fairly smooth in wavelength response .. so treat these separately, so can do fewer interpolations
(at each wavelength). Then multiply by more sharply varying wavelength features (with no x/y variation) in atmosphere.
Also multiply by 'base' transmission curves (lens/optics/filter/detector) if not included in 'variations' above.
(not taking advantage of this feature yet)

Maybe include clouds here too? 
(not included here : clouds and gain .. gain is per amp & no interpolation, clouds ...)

        # Tangent plane is above the telescope. Rotation angle does not matter.
        # Focal plane is in the telescope - includes rotation angle. (so telescope dependent nonuniformities        
        # like from the spiders have to be corrected to remove rotation angle before applying in Spatial Bandpass). 
        # CCD x/y is position inside each chip. Also generate ccdname and ampname.
"""

import os
import warnings
from copy import deepcopy
import numpy
from scipy import interpolate
from Bandpass import Bandpass
from Sed import Sed
from lsst.sims.atmosphere.transmission.modtranCards import ModtranCards
from lsst.sims.atmosphere.clouds.Clouds import Clouds

# The following *wavelen* parameters are default values for gridding wavelen/sb/flambda.
MINWAVELEN = 300
MAXWAVELEN = 1150
WAVELENSTEP = 0.1

# Some other typical default values for LSST.
EXPTIME = 15                      # Default exposure time. 
NEXP = 2                          # Default number of exposures.
EFFAREA = numpy.pi*(6.5*100/2.0)**2  # Default effective area of primary mirror..
GAIN = 2.3                        # Default gain. 
RDNOISE = 5                       # Default value - readnoise electrons or adu per pixel (per exposure)
DARKCURRENT = 0.2                 # Default value - dark current electrons or adu per pixel per second
OTHERNOISE = 4.69                 # Default value - other noise electrons or adu per pixel per exposure
PLATESCALE = 0.2                  # Default value - "/pixel
SEEING = {'u': 0.77, 'g':0.73, 'r':0.70, 'i':0.67, 'z':0.65, 'y':0.63}  # Default seeing values (in ")

class SpatialBandpass:

    def __init__(self, xdim_ccd=4000, ydim_ccd=4072, rad_fov=1.75, 
                 wavelen_min=MINWAVELEN, wavelen_max=MAXWAVELEN, wavelen_step=WAVELENSTEP):
        """Instantiate object. """
        self._xdim_ccd = xdim_ccd
        self._ydim_ccd = ydim_ccd
        self._rad_fov = rad_fov
        self._wavelen_min = wavelen_min
        self._wavelen_max = wavelen_max
        self._wavelen_step = wavelen_step
        return

    def setDefaults(self, bandpass, modtranParams, chipTempOffsets, vignetteFile=None, cloudSeed=42):
        """Convenience function to set spatial bandpass parameters, with the typical inputs expected for 'sim1'.

        This means: in a single filter, using modtranParameters for the atmosphere, """
        self.setupBaseThroughput(bandpass=bandpass, modtranParams=modtranParams)
        self.setDetectorQEVariation(chipTempOffsets)
        self.setVignette(vignetteFile)
        self.setClouds(cloudSeed)
        return


    def calcADU(self, sed, raft_chip, ccd_xy, expTime=EXPTIME, effarea=EFFAREA, gain=GAIN):
        """For a single object, calculate the adu expected at this position for this SED. """
        # Expect ccd_xy 
        # Generate the QE(T), T(x,y) of the detector.
        sb_qe = self.getQEVariation(raft_chip, ccd_xy)
        # Multiply by the base throughput curves.
        self.sb = self.base.sb * sb_qe
        # Generate the gray-scale variations.
        
        # Calculate the ADU values from this SED.
        if (numpy.not_equal(sed.wavelen, self.wavelen)):
            sed.resampleSED(wavelen_match=self.wavelen)
        # Calculate ADU
        wavelen = self.wavelen
        fnu = self.fnu
        # Calculate the number of photons.
        dlambda = wavelen[1] - wavelen[0]
        # Nphoton in units of 10^-23 ergs/cm^s/nm. 
        nphoton = (fnu / wavelen * bandpass.sb).sum()
        adu = nphoton * (expTime * effarea/gain) * (1/ERGSETC2JANSKY) * (1/PLANCK) * dlambda 
        return adu

    # Base throughput curves, constant across FOV.
    
    def setBaseThroughput(self, bandpass='r', baseDataDir=None, baseComponents=None,
                          modtranParams=None, useAtmosphereFile=None):
        """Set up the basic throughput curve applicable to the entire field of view.
    
        This includes the baseline detector curve, lens and mirror curves, filter curve and
        atmosphere transmission curve, which will be read from baseDataDir (which will be gathered
        from 'LSST_THROUGHPUTS_DEFAULT' if not specified).
        If useAtmosphere is not specified, the specified modtran parameters (a dictionary containing
        atmospheric parameters) will be used to create a new one with MODTRAN."""
        # Generate the basic wavelength-dependent throughput curve for this object.
        self.wavelen = numpy.arange(self._wavelen_min, self._wavelen_max+self._wavelen_step, self._wavelen_step)
        # Check the baseDataDir value. 
        if baseDataDir == None:
            baseDataDir = os.getenv('LSST_THROUGHPUTS_DEFAULT')
            if baseDataDir == None:
                raise Exception('baseDataDir not given and unable to access environment variable LSST_THROUGHPUTS_DEFAULT')
        # Check the component list. 
        if baseComponents == None:
            baseComponents = ['lens1.dat', 'lens2.dat', 'lens3.dat', 'm1.dat', 'm2.dat', 'm3.dat',
                              'detector.dat', 'filter_'+bandpass+'.dat']
        # Combine these base components (coming from a file). 
        self.base = Bandpass()
        self.base.readThroughputList(componentList=baseComponents, rootDir=baseDataDir,
                                     wavelen_min=self._wavelen_min, wavelen_max=(self._wavelen_max+self._wavelen_step),
                                     wavelen_step=self._wavelen_step)
        # Now generate the atmosphere (or possibly read atmosphere from file).
        if (useAtmosphereFile != None) & (modtranParams != None):
            warnings.warn('Only one of modtranParam or useAtmosphereFile should be provided. Will use file %s'
                          %(useAtmosphereFile))
        # (Deal with the atmosphere). If atmospheric transmission curve specified by file: 
        if useAtmosphereFile != None:
            atm = Bandpass()
            atm.readThroughput(useAtmosphereFile, self._wavelen_min,
                               self._wavelen_max+self._wavelen_step, self._wavelen_step)
        # (Deal with the atmosphere). If atmospheric transmission curve should be generated by modtran:
        else:
            if modtranParams == None:
                raise Exception('Please provide either the modtran parameters (modtranParam) or an atmosphere filename.')
            if isinstance(modtranParams, dict) == False:
                raise Exception('Modtran parameters should be in the form of a dictionary.')
            mc = ModtranCards()
            # Read default templates and format files.
            mc.setDefaults()
            # Write card to scratch directory. 
            scratchDir = os.getenv('SCRATCH_DIR')
            if scratchDir == None:
                scratchDir = '.'
            modfile = os.path.join(scratchDir, 'tmp')
            mc.writeModtranCards(paramValues=modtranParams, outfileRoot=modfile)
            mc.runModtran(modfile)
            atm = Bandpass()
            atm.readThroughput(modfile+'.plt', self._wavelen_min, self._wavelen_max+self._wavelen_step, self._wavelen_step)
            mc.cleanModtran(modfile)
            del(mc)
        # Combine the atmosphere with the base throughput.
        wavelen, sb = self.base.multiplyThroughputs(atm.wavelen, atm.sb)
        self.base.setBandpass(wavelen, sb, self._wavelen_min, self._wavelen_max+self._wavelen_step, self._wavelen_step)
        return


    # Wavelength-dependent changes. (generally can only be applied a single-object at a time).

    def setDetectorQEVariation(self, chipTempOffsets, tilt=True, qeTempFile=None, tempXYFiles=None):
        """Set up the temperature variation as a function of wavelength, and the chip temperature
        as a function of x/y in each chip.
        Assumes that the QE(T) is specified in data/QE_T.dat and that data also carries expected files
        for temperature variation within a raft, also that 'chip' names (including raft) in input dictionary can
        be translated to chip names for the T(x/y) files.
        if 'tilt' is True, then will use the 'tilted' chip T(x/y) files. """
        # Read QE(T) data.
        if qeTempFile == None:
            dataDir = os.path.join(os.getenv('CATALOGS_MEASURES_DIR'), 'data')
            qeTempFile = os.path.join(dataDir, 'QE_T.dat')
        file = open(qeTempFile, 'r')
        wavelen = []
        temperature = []
        qe_rel = []
        for line in file:
            if line.startswith("#"):
                continue
            values = line.split()
            wavelen.append(values[0])
            temperature.append(values[1])
            qe_rel.append(values[2])
        file.close()
        # Convert to arrays.
        wavelen = numpy.array(wavelen, dtype='float')
        temperature = numpy.array(temperature, dtype='float')
        qe_rel  = numpy.array(qe_rel, dtype='float')        
        self.qeT_interp = interpolate.LinearNDInterpolator((wavelen, temperature), qe_rel)
        # (you can then call this like follows .. )
        # w2 = numpy.arange(850, 1100, 1)
        # t2 = numpy.arange(-135, -80, 1)
        # (wavelen_grid, temperature_grid) = numpy.meshgrid(w2, t2)
        # qe_grid = qe_interp((wavelen_grid, temperature_grid))        
        # Read Temperature(X/Y) for each chip.
        # First find chips which we will be using (might be whole set or subset). 
        self.chiplist = []
        for c in chipTempOffsets.keys():
            # Assume that chiplist is given as RAFTID_CCDID
            # (in cameraGeometry file it is usually R:4,3 S:2,2 so might be some parsing somewhere)
            chip = c.split('_')[1]
            if chip not in self.chiplist:
                self.chiplist.append(chip)
        # Read the files specifying temperature profile variations across each CCD. 
        self.tempProfile_interp = {}
        if tilt == True:
            ccdstyle = 'Tilt'
        else:
            ccdstyle = 'Flat'
        for ccd in self.chiplist:
            if tempXYFiles == None:
                filename = os.path.join(dataDir, 'Txy_' + ccdstyle + '_' + ccd + '.dat')
            else:
                filename = tempXYFiles[ccd]
            try:
                file = open(filename, 'r')
            except IOError:
                raise Exception("Could not find %s file, for chip number %s" % (filename, ccd))
            x = []
            y = []
            temperature = []
            i = 0
            j = 0
            for line in file:
                values = line.split()
                for i in range(len(values)):
                    x.append(i)
                    y.append(j)
                    temperature.append(values[i])
                    i += 1
                j += 1
            file.close()
            x = numpy.array(x, 'float')
            y = numpy.array(y, 'float')
            temperature = numpy.array(temperature, 'float')
            x = x/x.max() * self._xdim_ccd
            y = y/y.max() * self._ydim_ccd
            #self.tempProfile_interp[ccd] = interpolate.LinearNDInterpolator((x, y), temperature)
            self.tempProfile_interp[ccd] = interpolate.Rbf(x, y, temperature)
        # And store pointer to chip offset temperature values internally (although not keeping a copy). 
        self.chipTempOffsets = chipTempOffsets
        return

    def getQEvariation(self, raft_chip, ccd_xy):
        """Given raft and chip name, and CCD x/y location, generate the throughput transmission curve variation
        due to Temperature variation.
        Note this works for a *single* object at a time, but ccd_xy should be an array with second dimension =2."""
        # Apply a temperature offset specified by chipTempOffsets.
        tempOffset = self.chipTempOffsets[raft_chip]
        ccd = raft_chip.split('_')[1]
        # Determine the temperature at this particular x/y in the chip (x/y = 2-d numpy array).
        temp = self.tempProfile_interp[ccd](ccd_xy[0][0], ccd_xy[0][1])
        temperature = temp + tempOffset
        # Calculate the variation in QE as a function of wavelength at this location.
        sb_qe = self.qeT_interp(self.wavelen, temperature)
        return sb_qe

    def setFilter(self, filter_jitter):
        pass

    def getFilter(self, focalxy):
        pass

    def setColorFlat(self):
        pass

    def getColorFlat(self):
        pass


    # Grey-scale changes. (generally can be applied multiple objects at a time). 

    def setVignette(self, vignetteFile=None):
        """Set up the vignetting function, reading form file 'vignetteFile'."""
        if vignetteFile == None:        
            dataDir = os.path.join(os.getenv('CATALOGS_MEASURES_DIR'), 'data')
            vignetteFile = os.path.join(dataDir, 'Vignette.dat')
        file = open(vignetteFile, 'r')
        radius = []
        vignetteFrac = []
        for line in file:
            if line.startswith("#"):
                continue
            values = line.split()
            radius.append(values[0])
            vignetteFrac.append(values[1])
        file.close()
        # Convert to arrays.
        radius = numpy.array(radius, 'float')
        vignetteFrac = numpy.array(vignetteFrac, 'float')
        # Note that vignette radius is in degrees from the center. (i.e. goes from 0-2). 
        self.vignette_interp = interpolate.interp1d(radius, vignetteFrac, kind='cubic', bounds_error=False, fill_value=0.0)
        return

    def getVignette(self, focalxy):
        """Given focal plane x/y location in ...
        focalxy format - numpy array with second dimension = 2 (i.e. shape= [X,2]).
        Note that this function will work for individual objects or multiple. """        
        rad = numpy.sqrt(focalxy[:,0]**2 + focalxy[:,1]**2)        
        vigFrac = self.vignette_interp(rad)
        return vigFrac
        
    def setClouds(self, total_extinction, cloudSeed=42, cloudSampling=240):
        """Set up the cloud extinction."""
        windowsize = numpy.sqrt(2.*(2*self._rad_fov)**2)
        clouds = Clouds(windowsize, cloudSampling)
        # Set up the power spectrum for the clouds.                
        clouds.setupPowerSpectrum()
        # Generate the clouds in 'real space'. 
        clouds.DirectClouds(randomSeed=cloudSeed)
        # Scale clouds by total extinction.
        cloudimage = clouds.clouds * total_extinction
        # The cloud extinction image is just an 'image' numbered from 0/0 to samplesize/samplesize
        #  So we need to add x/y locations appropriate to stellar x/y locations.
        # PROBABLY NEED TO CHANGE THIS SCALING, APPROPRIATELY FOR TANGENT PLANE
        x, y = numpy.mgrid[0:cloudSampling, 0:cloudSampling]
        x = x / (numpy.max(x)/2.0) - 1.0  # if focal plane goes from -1 to 1
        y = y / (numpy.max(y)/2.0) - 1.0   # probably need to change this.
        cloudxy = numpy.column_stack([numpy.ravel(x), numpy.ravel(y)])
        self.cloud_interp = interpolate.NearestNDInterpolator(cloudxy, cloudimage.ravel())
        return

    def getClouds(self, tangentxy):
        """Generate the cloud extinction for points at tangentxy .. [X, 2] array. """        
        cloudExtinction = self.cloud_interp(tangentxy)
        #cloud_extinction = numpy.power(10.0, -0.4*cloud_extinction)
        return cloudExtinction

    
