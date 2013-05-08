import os
import glob
import numpy
import pylab
from scipy import interpolate
from lsst.sims.catalogs.measures.photometry.Bandpass import Bandpass

def read_files():
    tfiles = glob.glob('sensorQE*.txt')
    temps = []
    for t in tfiles:
        temps.append(float(t.split('_')[2].rstrip('.txt').lstrip('T')))
    temps = numpy.array(temps, dtype='float')
    detector = {}
    for file, t in zip(tfiles, temps):
        detector[t] = Bandpass()
        detector[t].readThroughput(file)
    base_detector = Bandpass()
    base_detector.readThroughput(os.path.join(
        os.getenv('LSST_THROUGHPUTS_DEFAULT'),
        'detector.dat'))
    qe_var = {}
    t_base = 173.1 # numpy.median(temps)
    eps = 1e-20
    for t in temps:
        qe_var[t] = Bandpass()
        qe_var[t].setBandpass(wavelen = detector[t].wavelen,
                              sb = numpy.where(detector[t_base].sb<eps, 0,
                                               detector[t].sb / detector[t_base].sb))
    return temps, t_base, detector, base_detector, qe_var


if __name__== "__main__":
    
    temps, t_base, detector, base_detector, qe_var = read_files()
    
    # Compare Andy's new curves (read above) with standard detector throughput
    # curve (in LSST_THROUGHPUTS, aka docushare version).
    
    pylab.figure()
    for t in temps:
        pylab.plot(detector[t].wavelen, detector[t].sb, label='%.1f' %(t))
        
        pylab.plot(base_detector.wavelen, base_detector.sb, color='k',
                   linestyle=':',
                   linewidth=2, label = 'Standard')
        
    pylab.xlabel('Wavelength (nm)')
    pylab.ylabel('QE')

    # Well, those are not quite the same. Andy's is smoother in places,
    # especially @ the drop near 380nm and has slightly lower throughput
    #  @ many wavelengths.
    # Will send inquiry.
    

    # Compare temperature dependence as function of wavelength to temperature
    #  dependence from PS curves.  (divide by basic throughput curve).

    pylab.figure()
    stemps = numpy.array([-85, -95, -105], 'float')
    stemps = stemps + 273
    for st in stemps:
        # find the closest actual temperature to this value
        i = numpy.where((st - temps) == abs(st-temps).min())[0][0]
        t = temps[i]
        pylab.plot(qe_var[t].wavelen, qe_var[t].sb, label='%.1f (%.1f C)' %(t, t-273))
    pylab.xlabel('Wavelength (nm)')
    pylab.ylabel('Relative QE')
    pylab.legend(numpoints=1, fancybox=True, loc='lower left')
    pylab.xlim(800, 1100)
    pylab.ylim(0.2,1.6)
    pylab.savefig('qe_delta_t.png')

    # The results from this look generally similar to the PS data, although
    # at the 'hot' end (188K = -87C) the PS data do not show as much of an
    #  increase  in QE at the red end as Andy's data. 

    # Generate interpolation method as previously done ..
#    wavelen, temperature = numpy.meshgrid(qe_var[t_base].wavelen, temps)
    
#    qe_interp = interpolate.LinearNDInterpolator((wavelen, temperature), qe_var[t].sb)
    #pylab.show()
   
