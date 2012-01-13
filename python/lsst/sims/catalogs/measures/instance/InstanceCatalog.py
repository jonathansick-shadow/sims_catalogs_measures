""" InstanceCatalog Class

    ajc@astro Feb 10 2010
    Class and methods that operate on an InstanceClass

"""

import numpy
import warnings
import sys
import itertools
import gzip
from copy import deepcopy

from lsst.sims.catalogs.measures.astrometry.Astrometry import *
import lsst.sims.catalogs.measures.photometry.photUtils as phot
from lsst.sims.catalogs.measures.photometry.Bandpass import *
from lsst.sims.catalogs.measures.photometry.Sed import *
import lsst.sims.catalogs.measures.photometry.EBV as ebv
import lsst.sims.catalogs.measures.photometry.Variability as variability 
#from lsst.sims.catalogs.measures.instance.CatalogDescription import *
from lsst.sims.catalogs.measures.instance.SiteDescription import *
from lsst.sims.catalogs.measures.instance.Metadata import *
from CatalogDescription import *

class InstanceCatalog (Astrometry):
    """ Class that describes the instance catalog for the simulations. 

    Instance catalogs include a dictionary of numpy arrays which contains 
    core data. Additional arrays can be appended as ancillary data are 
    derived

    Catalog types and Object types are defined in the CatalogDescription class
    catalogType =  TRIM, SCIENCE, PHOTCAL, DIASOURCE, MISC, INVALID
    objectType = Point, Moving, Sersic, Image, Artefact, MISC
    catalogTable is name of the database table queried
    dataArray dictionary of numpy arrays of data

    """

    def __init__(self, configFile=None, cd=None, md=None):
        """Create an InstanceClass

        Instantiate an InstanceClass with the catalog type set to invalid
        """
        if cd is None:
            self.catalogDescription = CatalogDescription(configFile)
        else:
            self.catalogDescription = cd
        self.site = SiteDescription()
        if md is None:
            self.metadata = Metadata(configFile)
        else:
            self.metadata = md

        self.catalogType = None
        self.neighborhoodType = None
        self.objectType = None
        self.catalogTable = ""
        self.dataArray = {}

    # dataArray operations    
    def addColumn(self, array, name):
        """Add a numpy array to dataArray and warn if it already exists """
        if name in self.dataArray:
            warnings.warn("Entry %s exists in dataArray" % name)
        self.dataArray[name] = array
    def deleteColumn(self, name):
        """Delete a  numpy array from dataArray """
        if self.dataArray.has_key(name):
            del self.dataArray[name] 
        else:
            warnings.warn("Entry %s does not exists in dataArray" % name)
            
    def noneCheck(self, sedPaths, sedFile):
        '''Parse any None in the sedFilenames'''
        try:
            return "%s.gz"%sedPaths[sedFile]
        except KeyError:
            return None
            
    def makeFilePaths(self, fileName):
        '''manipulate sedFilename to add path for file on disk'''
        sedPaths = self.catalogDescription.getPathMap()
        self.addColumn([self.noneCheck(sedPaths, x) for x in  self.dataArray[fileName]], fileName)
        
#        print "FILEPATHS"
#        print self.dataArray[fileName]

#        print [x for x in self.dataArray[fileName] if x != None]
#        self.addColumn(map(lambda x: sedPaths[x] + ".gz", self.dataArray[fileName]),fileName)


    # validate that the catalog contains the correct data
    def validateData(self, catalogType):
        """Validate that the class contains the correct attributes
        
        This does not test validity of the data
            This combines searches through the required and derived attributes"""

        # validate required and derived data 
        requiredAttributeList = self.catalogDescription.getRequiredFields(catalogType,
                                                                          self.neighborhoodType,
                                                                          self.objectType)

        for name in requiredAttributeList:
            if ((self.dataArray.has_key(name) == False)):
                raise ValueError("Entry %s does not exist in required data"%name)

        derivedAttributeList = self.catalogDescription.getDerivedFields(catalogType,
                                                                          self.neighborhoodType,
                                                                          self.objectType)
        for name in derivedAttributeList:
            if ((self.dataArray.has_key(name) == False)):
                raise ValueError("Entry %s does not exist in derived data" % name)

        return True


    def formatData(self, conversion):
        '''Format data based on conversion function

        conversion is a tuple of (attribute, function for conversion, array of indices in dataArray)
        '''

        # extract data value for attribute
        x = self.dataArray[conversion[0]][conversion[2]]
#        print x
#        if (x == None):
#            x=0
#        print conversion[0],conversion[1],eval(conversion[1]),conversion[2]
        return eval(conversion[1])
        
        
    # Output of formatted data catalogs
    def writeCatalogData(self, filename, catalogType, newfile = False, compress = False):
        """Write an instanceCatalog dataArray based on the catalog type given

           If the catalogType is TRIM use the objectType to determine the output format
           Provide the option to clobber the file
        """
        # open file
        if compress:
  	    if (newfile == False):
	        outputFile = gzip.open(filename+".gz","a")
	    else:
	        outputFile = gzip.open(filename+".gz","w")
        else:
  	    if (newfile == False):
	        outputFile = open(filename,"a")
	    else:
	        outputFile = open(filename,"w")

        # Determine the catalogType and objectType for printing
        format, attributeList, conversion = self.catalogDescription.getFormat(catalogType, self.objectType)

        #write trim file based on objectType


        #output header info for reference catalog
        if (catalogType == 'REFERENCECATALOG'):
            outputFile.write("# ")
            outputFile.write(" ".join(attributeList))
            outputFile.write("\n")
        
        # add newline to format string - configobj does not interpret these correctly
        format = format +"\n"
        for i in range(len(self.dataArray["id"])):
            # use list comprehension to output all attributes in the given format string
            # 2.6 outputFile.write(formatString,(map(lambda x: self.dataArray[x][i],attributeList)))
            # outputFile.write(format % tuple(map(lambda x: self.dataArray[x][i],attributeList)))
#            print attributeList
#            print conversion
#            print format
            outputFile.write(format % tuple([self.formatData(x) for x in
                                             zip(attributeList,conversion,numpy.ones(len(attributeList),dtype=int)*i)]))


        outputFile.close()

    # Composite astrometry operations
    def makeHelio(self):
        """ Generate Heliocentric coordinates """

        # apply precession
        raOut, decOut = self.applyPrecession(self.dataArray['raJ2000'], self.dataArray['decJ2000'],
                                                   MJD = self.metadata.parameters['Opsim_expmjd'])

        # apply proper motion (note parallax is converted to arcseconds)
        raOut, decOut = self.applyProperMotion(raOut, decOut, self.dataArray['properMotionRa'],
                                               self.dataArray['properMotionDec'], self.dataArray['parallax']/1000.,
                                               self.dataArray['radialVelocity'], MJD = self.metadata.parameters['Opsim_expmjd'])

        # TODO 3/29/2010 convert FK5 to ICRS?
        self.addColumn(raOut, 'raHelio')
        self.addColumn(decOut, 'decHelio')

    def makeApparent(self):
        """ Generate apparent coordinates
        

        This converts from the J2000 coordinates to the position as
        viewed from the center of the Earth and includes the effects
        of light defection (ignored), annual aberration, precession
        and nutation
        """
        # (note parallax is converted to arcseconds)
        raOut, decOut = self.applyMeanApparentPlace(self.dataArray['raJ2000'], self.dataArray['decJ2000'],
                                                    self.dataArray['properMotionRa'],
                                                    self.dataArray['properMotionDec'], self.dataArray['parallax']/1000.,
                                                    self.dataArray['radialVelocity'],
                                                    MJD=self.metadata.parameters['Opsim_expmjd'])

        self.addColumn(raOut, 'raApp')
        self.addColumn(decOut, 'decApp')


    def makeObserved(self):
        """ Generate Observed coordinates

        From the observed coordinates generate the position of the
        source as observed from the telescope site. This includes the
        hour angle, diurnal aberration, alt-az, and refraction. 
        """
        if ((("raApp" in self.dataArray) and 
             ("decApp" in self.dataArray)) != True):
            self.makeApparent()
        raOut, decOut = self.applyMeanObservedPlace(self.dataArray['raApp'], self.dataArray['decApp'],
                                                    MJD=self.metadata.parameters['Opsim_expmjd'])

        self.addColumn(raOut, 'raObs')
        self.addColumn(decOut, 'decObs')

    def calculateUnrefractedAltAz(self):
        '''Calculate the unrefracted AltAz for the telescope given the ra dec opsim pointing'''
        #Calculate pointing of telescope in observed frame and the rotation matrix to transform to this position
        raCenter, decCenter, altCenter, azCenter = self.transformPointingToObserved(
            self.metadata.parameters['Unrefracted_RA'],
            self.metadata.parameters['Unrefracted_Dec'],
            includeRefraction = False)

        #Update the meta data for ALt-Az
        self.metadata.addMetadata("Unrefracted_Altitude", altCenter,\
                                  "Opsim value of the altitude of the observation")
        self.metadata.addMetadata("Unrefracted_Azimuth", azCenter,\
                                  "Opsim value of the azimuth of the observation")

    def makeEBV(self, Rv=3.1):
        '''Populate E(B-V) values from gLon, gLat'''

        glon, glat = self.equatorialToGalactic(self.dataArray['raJ2000'],self.dataArray['decJ2000'])

        datadir = os.environ.get("CAT_SHARE_DATA")

        ebvMapNorth = ebv.EbvMap()
        ebvMapNorth.readMapFits(os.path.join(datadir, "data/Dust/SFD_dust_4096_ngp.fits"))

        ebvMapSouth = ebv.EbvMap()
        ebvMapSouth.readMapFits(os.path.join(datadir, "data/Dust/SFD_dust_4096_sgp.fits"))
			            
        self.addColumn(ebv.calculateEbv(glon, glat, ebvMapNorth, ebvMapSouth, interp = True)*Rv, 'galacticAv')
        self.addColumn(numpy.ones(len(self.dataArray['galacticAv']))*Rv, 'galacticRv')
        self.addColumn(numpy.array(['CCM' for val in
                                    range(len(self.dataArray['galacticAv']))]),
                       'galacticExtinctionModel')

    def makeTrimCoords(self):
        """ Generate TRIM coordinates

        From the apparent coordinates generate the position of the
        source as observed from the telescope site (required for the
        trim files). This includes the hour angle, diurnal aberration,
        alt-az. This does NOT include refraction.
        """
        
        # calculate E(B-V) parameters if Extragalactic
        if (self.neighborhoodType == 'EXTRAGALACTIC'): 
            self.makeEBV()
                
        #Calculate pointing of telescope in observed frame and the rotation matrix to transform to this position
        raCenter, decCenter, altCenter, azCenter = self.transformPointingToObserved(
            self.metadata.parameters['Unrefracted_RA'],
            self.metadata.parameters['Unrefracted_Dec'],
            includeRefraction = False)

        # update file paths for place on disk
        self.makeFilePaths('sedFilename')

        #Update the meta data for ALt-Az
#        self.metadata.addMetadata("Unrefracted_Altitude", altCenter,\
#                                  "Opsim value of the altitude of the observation")
#        self.metadata.addMetadata("Unrefracted_Azimuth", azCenter,\
#                                  "Opsim value of the azimuth of the observation")
        
        
        xyzJ2000 = self.sphericalToCartesian(self.metadata.parameters['Unrefracted_RA'],
                                           self.metadata.parameters['Unrefracted_Dec'])
        xyzJ2000 /= math.sqrt(numpy.dot(xyzJ2000, xyzJ2000))

        xyzObs = self.sphericalToCartesian(raCenter, decCenter)
        xyzObs /= math.sqrt(numpy.dot(xyzObs, xyzObs))

        rotationMatrix = self.rotationMatrixFromVectors(xyzObs, xyzJ2000)

        #convert positions of sources from reference to observed
        if ((("raApp" in self.dataArray) and 
             ("decApp" in self.dataArray)) != True):
            self.makeApparent()
            raOut, decOut = self.applyMeanObservedPlace(self.dataArray['raApp'],
                                                     self.dataArray['decApp'],
                                                     MJD=self.metadata.parameters['Opsim_expmjd'],
                                                     altAzHr=False, includeRefraction=False)


        # correct for pointing of the telescope (so only have differential offsets)
        # and set output as ICRS
        xyz = self.sphericalToCartesian(raOut, decOut).transpose()
        for i,_xyz in enumerate(xyz):
            xyzNew = numpy.dot(rotationMatrix,_xyz)
            raOut[i], decOut[i] = self.cartesianToSpherical(xyzNew)

        self.addColumn(raOut, 'raTrim')
        self.addColumn(decOut, 'decTrim')

    def makeReferenceCoords(self):
        '''Generate reference file attributes'''

        # generate photometry and objid
        datadir = os.path.join(os.environ.get("CAT_SHARE_DATA"),"data")
        if (self.neighborhoodType == 'GALACTIC'):
            # generate file paths for SEDS
            self.makeFilePaths('sedFilename')
            self.calculateStellarMagnitudes(dataDir=datadir)
            #generate isVar - use varsimobjid (check for None), bit shift id<<10 and subtract to get variability number
            variableArray = [x[0] if x[0] is not None else x[1]<<10 for x in zip(self.dataArray['varsimobjid'],
                                                                self.dataArray['id'])]
#            print self.dataArray['varsimobjid'][0:10]
#            print variableArray
#            print self.dataArray['id'][0:10]
            self.addColumn((variableArray  - (self.dataArray['id']<<10)),'isVar')

            #generate ID and isVar - id is simobjid <<10 + appendint
            self.dataArray['id'] = (self.dataArray['id']  << 10) + self.dataArray['appendint']
            

            
        elif (self.neighborhoodType == 'EXTRAGALACTIC'):
            # generate isVar - use varsimobjid (check for None), bit shift id<<10
            # and subtract to get variability number
            variableArray = [x[0] if x[0] is not None else x[1]<<10 for x in zip(self.dataArray['varsimobjid'],
                                                                                     self.dataArray['id'])]
            self.addColumn((variableArray  - (self.dataArray['id']<<10)),'isVar')

            #generate ID for galaxies - id is simobjid <<10 + appendint
            self.dataArray['galtileid'] = (self.dataArray['galtileid']  << 10) + self.dataArray['appendint']

            #generate EBV,Av, values
            #self.makeEBV()
            #derive filepaths:
            self.makeFilePaths('bulgeSedFilename')
            self.makeFilePaths('diskSedFilename')
            self.makeFilePaths('agnSedFilename')
            #generate magnitudes
            self.calculateGalaxyMagnitudes(dataDir=datadir)


    def transformPointingToObserved(self, ra, dec, includeRefraction = False):
        """Take an LSST central pointing and determine is observed position on the sky
        (excluding refraction)"""
        raOut, decOut = self.applyMeanApparentPlace([ra], [dec], [0.], [0.],[0.], [0.],
                                                    MJD=self.metadata.parameters['Opsim_expmjd'])
        #print ra,dec,raOut, decOut,self.metadata.parameters['Opsim_expmjd']
        raObs, decObs, altObs, azObs = self.applyMeanObservedPlace(raOut, decOut, MJD=self.metadata.parameters['Opsim_expmjd'],
                                                                   altAzHr=True, includeRefraction = includeRefraction)

        return raObs[0], decObs[0], altObs[0], azObs[0]


    # Photometry composite methods

    def calculateStellarMagnitudes(self, filterList=('u', 'g', 'r', 'i', 'z', 'y'),
                                   dataDir = None, filterroot='total_'):
        """For stellar sources and a list of bandpass names generate magnitudes"""
    
        # load bandpasses
        bandpassDict = phot.loadBandpasses(filterlist=filterList, dataDir = None)
        
        #load required SEDs
        sedDict = phot.loadSeds(self.dataArray["sedFilename"], dataDir=dataDir)

        # pick one SED to be reference - randomly choose #1 sed
        refsed = sedDict[self.dataArray["sedFilename"][0]]
        if (refsed.needResample(wavelen_match=bandpassDict[filterList[0]].wavelen)):
            refsed.resampleSED(wavelen_match=bandpassDict[filterList[0]].wavelen)

        # need to put all SEDs on same wavelength grid
        for sed in sedDict.values():
            if (sed.needResample(wavelen_match=refsed.wavelen)):
                sed.resampleSED(wavelen_match=refsed.wavelen)

        # Calculate dust parameters for all stars  (a_x/b_x have implicit wavelength dep)
        a_x, b_x = refsed.setupCCMab()

        # generate magnitude arrays and set to zero
        for name,bandpass in bandpassDict.items():
            self.addColumn(numpy.zeros(len(self.dataArray["id"])), name+"Derived")

        # set up for calculating mags in all bandpasses
        phiarray, wavelenstep = phot.setupPhiArray_dict(bandpassDict, filterList)

        # loop through magnitudes and calculate   
        for i,sedName in enumerate(self.dataArray["sedFilename"]):
            sed = deepcopy(sedDict[sedName])
            sed.multiplyFluxNorm(self.dataArray["fluxNorm"][i])
            sed.addCCMDust(a_x, b_x, A_v=self.dataArray["galacticAv"][i], R_v=float(self.dataArray["galacticRv"][i]))
            sed.flambdaTofnu()

            mags = sed.manyMagCalc(phiarray, wavelenstep)
            for j,name in enumerate(filterList):
                self.dataArray[name+"Derived"][i] = mags[j]
            del sed


    def calculateGalaxyMagnitudes(self, filterList=('u', 'g', 'r', 'i', 'z', 'y'),
                                   dataDir = None, filterroot='total_'):
        """For galaxy sources and a list of bandpass names generate magnitudes"""

        # load bandpasses
        bandpassDict = phot.loadBandpasses(filterlist=filterList, dataDir = None)
        
        #load required SEDs ignore None
        sedDict = phot.loadSeds(self.dataArray["bulgeSedFilename"], dataDir=dataDir)
        sedDict.update(phot.loadSeds(self.dataArray["diskSedFilename"], dataDir=dataDir))
        sedDict.update(phot.loadSeds(self.dataArray["agnSedFilename"], dataDir=dataDir))

        # pick one SED to be reference and grid all to common wavelength
        refsed = sedDict.values()[0]
        for sed in sedDict.values():
            if (sed.needResample(wavelen_match=refsed.wavelen)):
                sed.resampleSED(wavelen_match=refsed.wavelen)

        # generate magnitude arrays and set to zero
        for name,bandpass in bandpassDict.items():
            self.addColumn(numpy.zeros(len(self.dataArray["id"])), name+"Derived")

        # set up for calculating mags in all bandpasses
        phiarray, wavelenstep = phot.setupPhiArray_dict(bandpassDict, filterList)

        # loop through magnitudes and calculate   
        for i,id in enumerate(self.dataArray["id"]):
            try:
                sedBulge = deepcopy(sedDict[self.dataArray["bulgeSedFilename"][i]])
                sedBulge.multiplyFluxNorm(self.dataArray["bulgeFluxNorm"][i])
                a_x, b_x = sedBulge.setupCCMab()
                sedBulge.addCCMDust(a_x, b_x, A_v=self.dataArray["internalAv_b"][i])
                sedBulge.redshiftSED(redshift=self.dataArray["redshift"][i], dimming=True)
#                print sedBulge.wavelen.min(),sedBulge.wavelen.max()
                sedBulge.resampleSED(wavelen_match=bandpassDict['u'].wavelen)
#                print "BULGE EXISTS",self.dataArray["bulgeSedFilename"][i]
#                print "DUST", self.dataArray["internalAv_b"][i]
            except KeyError:
                sedBulge = None

            try:
                sedDisk = deepcopy(sedDict[self.dataArray["diskSedFilename"][i]])
                sedDisk.multiplyFluxNorm(self.dataArray["diskFluxNorm"][i])
                a_x, b_x = sedDisk.setupCCMab()
                sedDisk.addCCMDust(a_x, b_x, A_v=self.dataArray["internalAv_d"][i])
                sedDisk.redshiftSED(redshift=self.dataArray["redshift"][i], dimming=True)
#                print sedDisk.wavelen.min(),sedDisk.wavelen.max()
                sedDisk.resampleSED(wavelen_match=bandpassDict['u'].wavelen)
#                print "DISK EXISTS",self.dataArray["diskSedFilename"][i]
#                print "DUST", self.dataArray["internalAv_d"][i]
            except KeyError:
                sedDisk = None

            try:
                sedAgn = deepcopy(sedDict[self.dataArray["agnSedFilename"][i]])
                #IS THIS CORRECT
                if (self.dataArray["agnFluxNorm"][i] > 90.):
                    self.dataArray["agnFluxNorm"][i] = 0.0
                    sedAgn = None
                sedAgn.multiplyFluxNorm(self.dataArray["agnFluxNorm"][i])
                sedAgn.redshiftSED(redshift=self.dataArray["redshift"][i], dimming=True)
                sedAgn.resampleSED(wavelen_match=bandpassDict['u'].wavelen)
#                print "AGN EXISTS",self.dataArray["agnSedFilename"][i]

            except KeyError:
                sedAgn = None
                

            #print "Galactic"
            #print self.dataArray["galacticAv"][i]
            # add disk and agn to bulge - allow for None in data

            sedList = [sedBulge, sedDisk, sedAgn]
            sedListNoNone = [x for x in sedList if x != None] 
            sed = sedListNoNone[0]
#            print "SED LIST", sedList
#            print "SED NoNone LIST",sedListNoNone

            if (sedDisk != None):
                sedDisk.flambdaTofnu()
#                print "DISK MAG", sedDisk.manyMagCalc(phiarray, wavelenstep)
#                print "DISK g MAG", sedDisk.calcMag(bandpassDict['g'])
                
            if (sedBulge != None):
                sedBulge.flambdaTofnu()
#                print "BULGE MAG", sedBulge.manyMagCalc(phiarray, wavelenstep)
#                print "BULGE g MAG", sedBulge.calcMag(bandpassDict['g'])

            if (sedAgn != None):
                sedAgn.flambdaTofnu()
#                print "AGN MAG", sedAgn.manyMagCalc(phiarray, wavelenstep)
#                print "AGN  g MAG", sedAgn.calcMag(bandpassDict['g'])

            #print sed.flambda.mean()
            #print mags
#            print "SUMMING"
#            print "INITIAL",sed.flambda.mean()
            for nextSed in sedListNoNone[1:]:
#                print nextSed.flambda.mean()
                sed.flambda += nextSed.flambda
                sed.fnu += nextSed.fnu
#            print "FINAL",sed.flambda.mean()

            # apply galactic extinction
            #a_mw, b_mw = sed.setupCCMab()
            #sed.addCCMDust(a_mw, b_mw, A_v=self.dataArray["galacticAv"][i])
            #sed.flambdaTofnu()
            
            mags = sed.manyMagCalc(phiarray, wavelenstep)
#            print "SUMMED MAGS", mags

#            print "DB MAGS", self.dataArray["id"][i],self.dataArray["u"][i],self.dataArray["g"][i],self.dataArray["r"][i],self.dataArray["i"][i],self.dataArray["z"][i],self.dataArray["y"][i]
            for j,name in enumerate(filterList):
                self.dataArray[name+"Derived"][i] = mags[j]
            del sed

    def applyVariability(self):
        """Apply variability models to magnitude normalization constants.
        The name of the method defined in the Variability class along with 
        the parameters for applying the variability are stored in the database
        as serialized dictionary objects
        """
        var = variability.Variability(cache=True)
        filters = ['u', 'g', 'r', 'i', 'z', 'y','U','G','R','I','Z','Y']
        #Map to translate filter integer to filter character.
        filterMap = {0:"u", 1:"g", 2:"r", 3:"i", 4:"z", 5:"y"}
        filt = ''
        expmjd = None
        if self.metadata.parameters.has_key('Opsim_filter'):
            filt = self.metadata.parameters['Opsim_filter']
        elif self.metadata.parameters.has_key('filter'):
            filt = self.metadata.parameters['filter']
        else:
            raise Exception("Need filter information in metadata to procede with variability application")

        if self.metadata.parameters.has_key('Opsim_expmjd'):
            expmjd = self.metadata.parameters['Opsim_expmjd']
        elif self.metadata.parameters.has_key('expmjd'):
            expmjd = self.metadata.parameters['expmjd']
        else:
            raise Exception("Need time information in metadata to procede with variability application")

        if filt in filters:
            filt = filt.lower()
        elif filterMap.has_key(filt):
            filt = filterMap[filt]
        else:
            raise Exception("Filter %s does not match the LSST filter list"\
                    %(str(self.metadata.parameters['Opsim_filter'])))
        
        #Apply variability to an entire array of magnitude normalization
        #constants.
        if not self.dataArray.has_key('variabilityParameters'):
            raise Exception("Cannot apply variability without parameters")
        numVar = len(self.dataArray['variabilityParameters'])
        magOffset = numpy.zeros(numVar)
        for ind,dstr in zip(range(numVar),self.dataArray['variabilityParameters']):
            dstr = str(dstr)
            if dstr == 'None' or dstr == '':
                continue
            d = eval(dstr)
            magOffset[ind] = eval("var.%s(d['pars'], \
                %f)['%s']"%\
                (d['varMethodName'],expmjd,filt))
        self.dataArray["magNorm"] += magOffset

""" TODO (2/18/2010) incorporate the precession routines
    def makeMeasured(self):
        raOut, decOut = self.applyPropermotion(self.dataArray['raJ2000'], 
                                    self.dataArray['decJ2000'])
        raOut, decOut = self.applyParallax(raOut, decOut)
        self.addColumn(raOut, 'raMeasured')
        self.addColumn(decOut, 'decMeasured')
    def makeGeo(self):
        if ((("raHelio" in self.dataArray) and 
             ("decHelio" in self.dataArray)) != True):
            self.makeHelio()
        raOut, decOut = self.applyParallax()
        raOut, decOut = self.applyAberration()
        self.addColumn(raOut, 'raGeo')
        self.addColumn(decOut, 'decGeo')
    def makeTopo(self):
        if ((("raGeo" in self.dataArray) and 
             ("decGeo" in self.dataArray)) != True):
            self.makeGeo()
        raOut, decOut = self.applyAbsoluteRefraction(raPar, decPar)
        self.addColumn(raOut, 'raHTopo')
        self.addColumn(decOut, 'decTopo')
     """



