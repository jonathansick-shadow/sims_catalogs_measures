import os
import numpy
import pylab
from lsst.sims.catalogs.measures.photometry.SpatialBandpass import SpatialBandpass

import time
def dtime(time_prev):
   return (time.time() - time_prev, time.time())

# Set up some random test data
npts= 1000
xrange = 2000.0
x = numpy.random.uniform(low=-xrange/2.0, high=xrange/2.0, size=npts) 
y = numpy.random.uniform(low=-xrange/2.0, high=xrange/2.0, size=npts) 
test_xy = numpy.column_stack([x,y])
raft = [1, 2, 3]
chip = [1, 2, 3]
raft_chip = []
for r1 in raft:
   for r2 in raft:
      for c1 in chip:
         for c2 in chip:
            raft_chip.append('%s%s_%s%s' %(r1,r2,c1,c2))
n = numpy.random.random_integers(size=npts, low=0, high=len(raft_chip)-1)
test_raftchip = []
for i in range(npts):
   test_raftchip.append(raft_chip[n[i]])
   

def setBase():
   # Read modtran atmospheric parameters for testing.
   file = open('test_modtranparams', 'r')
   modtranParamsList = []
   keys = []
   for line in file:
      if line.startswith('#'):
         line = line.lstrip('#')
         keys = line.split()
         continue
      values = line.split()
      modtranDict = {}
      for k, v in zip(keys, values):
         modtranDict[k] = v
      modtranParamsList.append(modtranDict)
   file.close()
   # set up base throughput (including modtran atmospheric generation)
   t=time.time()
   spb = SpatialBandpass()
   dt, t = dtime(t)
   print "Instantiating spatial bandpass: %f s" %(dt)
   spb.setBaseThroughput(bandpass='r', modtranParams=modtranParamsList[0])
   dt, t = dtime(t)
   print "Setting base throughput values (including running MODTRAN): %f s" %(dt)
   return spb

def testQE(spb):
   # Set up some temperature variation offset parameters for testing.
   chipTempOffsets = {}   
   for rc in raft_chip:
      chipTempOffsets[rc] = 13.0*numpy.random.uniform()
   # Set QE variations.
   t = time.time()
   spb.setDetectorQEVariation(chipTempOffsets)
   dt, t = dtime(t)
   print "Setting QE variations across chip: %f s" %(dt)
   # Test the QE variability.
   txy = numpy.vsplit(test_xy, len(test_xy))
   for i in range(len(test_xy)):
      sb_qe = spb.getQEvariation(test_raftchip[i], txy[i])
   dt, t = dtime(t)
   print "Calculating QE variations for %s objects: %f s" %(len(test_xy), dt)
   return spb

def plotQE(spb):
   chiplist = []
   for c1 in chip:
      for c2 in chip:
         chiplist.append('%s%s' %(c1,c2))
   chipTempOffsets = {}
   print "Using %d raft/chip combos, %d chips" %(len(raft_chip), len(chiplist))
   for rc in raft_chip:
      chipTempOffsets[rc] = 0.0
   spb.setDetectorQEVariation(chipTempOffsets, tilt=True)
   # Set up grid to plot change in temperature.
   x = numpy.arange(0, 500)
   y = numpy.arange(0, 500)
   ccd_xy = numpy.meshgrid(x, y)
   pylab.figure()
   pylab.subplots_adjust(left=0.2, right=0.8, wspace=0.1, hspace=0.05)
   for c, i in zip(chiplist, range(len(chiplist))):
      pylab.subplot(3,3,i)
      tempgrid = spb.tempProfile_interp[c](ccd_xy[0], ccd_xy[1])
      pylab.imshow(tempgrid, origin='lower') #, vmax=-99.9, vmin=-100.3)
      pylab.axis('off')
   pylab.suptitle('Temperature Variation Across Chips')
   pylab.show()
   return

def testVignette(spb):
   t = time.time()
   spb.setVignette()
   dt, t = dtime(t)
   print "Set vignette function across fov: %f s" %(dt)
   # test the vignette function
   rad = numpy.sqrt(test_xy[:,0]**2 + test_xy[:,1]**2)
   scale = 2.1 / rad.max()
   test_xy[:,0] = test_xy[:,0] * scale
   test_xy[:,1] = test_xy[:,1] * scale
   txy=numpy.vsplit(test_xy, len(test_xy))
   dt, t = dtime(t)
   for i in range(len(test_xy)):
      vigFrac = spb.getVignette(txy[i])
      #print vigFrac
   dt, t = dtime(t)
   print "Calculated individal vignette fractions: %f s" %(dt)
   vigFrac = spb.getVignette(test_xy)
   #print vigFrac
   dt, t = dtime(t)
   print "Calculated all vignette fractions: %f s" %(dt)
   return spb

def plotVignette(spb):
   spb.setVignette()
   x = numpy.arange(-200, 200, 1, dtype='float')
   y = numpy.arange(-200, 200, 1, dtype='float')
   (ix,iy) = numpy.meshgrid(x, y)
   ccd_xy = numpy.column_stack([ix.ravel(), iy.ravel()])
   rad = numpy.sqrt(ccd_xy[:,0]**2 + ccd_xy[:,1]**2)
   scale =  2.2 / rad.max()
   ccd_xy[:,0] = ccd_xy[:,0] * scale
   ccd_xy[:,1] = ccd_xy[:,1] * scale
   vigGrid = spb.getVignette(ccd_xy)
   rad = numpy.sqrt(ccd_xy[:,0]**2 + ccd_xy[:,1]**2)
   vigIm = vigGrid.reshape(len(x), len(y))
   pylab.figure()
   pylab.imshow(vigIm, origin='lower')
   pylab.colorbar()
   pylab.title("Vignetting")
   pylab.figure()
   pylab.plot(rad, vigGrid, 'b-')
   pylab.xlabel("Radius (deg)")
   pylab.ylabel("Vignetting Fraction")
   pylab.title('Vignetting')
   pylab.show()
   return 

def testClouds(spb):
   t = time.time()
   spb.setClouds(total_extinction=1.0)
   dt, t = dtime(t)
   print "Set cloud power spectrum and cloud image: %f s"%(dt)
   rad = numpy.sqrt(test_xy[:,0]**2 + test_xy[:,1]**2)
   scale = 0.9 / rad.max()
   test_xy[:,0] = test_xy[:,0]*scale
   test_xy[:,1] = test_xy[:,1]*scale
   txy=numpy.vsplit(test_xy, len(test_xy))
   dt, t = dtime(t)
   for i in range(len(test_xy)):
      cloudExt = spb.getClouds(txy[i])
      #print cloudExt
   dt, t = dtime(t)
   print "Calculated individual cloud extinction: %f s" %(dt)
   cloudExt = spb.getClouds(test_xy)
   #print cloudExt
   dt, t = dtime(t)
   print "Calculated all cloud extinction: %f s" %(dt)
   return spb


if __name__ == '__main__':
   
   spb = setBase()
   spb = testQE(spb)
   spb = testVignette(spb)
   spb = testClouds(spb)
   
   """
   pylab.figure()
   pylab.plot(spb.base.wavelen, spb.base.sb)
   adir = os.getenv('LSST_THROUGHPUTS_ATMOS')
   spb.setBaseThroughput(bandpass='r', useAtmosphereFile=os.path.join(adir, 'atmos_12.dat' ))
   pylab.figure()
   pylab.plot(spb.base.wavelen, spb.base.sb)
   pylab.show()
   """
   

