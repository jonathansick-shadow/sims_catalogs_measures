import pyfits
import math
import numpy
from darkenergy import *
from interpolation import *
from comoving_distance import *


class WL:
    """ Class for obtaining the shear parameters for a galaxy. 
    The shear parameters are read in from the WL maps, and 
    interpolated to the galaxy position in ra/dec/redshift.
    The WL maps used in this version consist of one original
    map, mirrored and repeated over the whole sky. This is a
    nonsense, cosmologically, but small areas (pointings)
    of the sky are usable (preferably if not crossing a map
    boundary). 
    """

    # some constants.
    NbinsX, NbinsY = 0, 0
    w0, wa, Omega_m, Omega_Lambda, h = 0, 0, 0, 0, 0
    survey_angle = 0
    number_of_maps = 3  # number of redshift bins.
    source_redshift, last_plane, comoving_distance = [], [], []
    shear1map, shear2map, convmap = [], [], []

    def __init__(self):
        data = None

    def initialize(self):
        """ To start with, we should get some of the info required
        from the WL maps, set up some variables and initialize this and that. 
        """

        hdulist = pyfits.open(
            "/astro/net/lsst1/shared/djbard/mirrored/WL-shear1_m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_4096xy_0001r_0029p_0100z_og.gre_X4_smaller.fit")

        hdr = hdulist[0].header
        self.NbinsX = hdr['MAP']
        self.NbinsY = self.NbinsX
        H_0 = hdr['H_0']
        self.h = H_0/100.0
        self.Omega_m = hdr['OMEGA_M']
        self.Omega_Lambda = hdr['OMEGA_L']
        self.w0 = hdr['W_0']
        self.wa = hdr['W_A']
        self.survey_angle = hdr['ANGLE']
        hdulist.close()

        # there are only three redshift planes - hardcode them.
        # Can also access them from the WL map header files.
        self.source_redshift = [1.0, 1.5, 2.0]
        self.last_plane = [29, 38, 46]
        self.comoving_distance = [2370.3, 3152.481, 3759.214]

        # load in the WL maps! this takes some time.

        # can I avoid initializing?
        self.shear1map = []
        self.shear2map = []
        self.convmap = []

        for i in range(0, self.number_of_maps):
            self.shear1map.append([])
            self.shear2map.append([])
            self.convmap.append([])
            for j in range(0, int(self.NbinsX*self.NbinsY)):
                self.shear1map[i].append(0)
                self.shear2map[i].append(0)
                self.convmap[i].append(0)

        for i in range(0, self.number_of_maps):

            shear1name = self.get_filename("shear1", i)
            shear2name = self.get_filename("shear2", i)
            convname = self.get_filename("conv", i)
            shear1hdu = pyfits.open(shear1name)
            temp_shear1map = shear1hdu[0].data
            shear1hdu.close()
            shear2hdu = pyfits.open(shear2name)
            temp_shear2map = shear2hdu[0].data
            shear2hdu.close()
            convhdu = pyfits.open(convname)
            temp_convmap = convhdu[0].data
            convhdu.close()
    # now, the maps have to be in a 1D array! not 2D matrix. bummer.

            for xx in range(0, int(self.NbinsX)):
                for yy in range(0, int(self.NbinsY)):
                    idx = int(self.NbinsY*xx + yy)

                    self.shear1map[i][idx] = temp_shear1map[xx][yy]
                    self.shear2map[i][idx] = temp_shear2map[xx][yy]
                    self.convmap[i][idx] = temp_convmap[xx][yy]
            del(temp_shear1map)
            del(temp_shear2map)
            del(temp_convmap)

        # initialize dark energy.
        initialize_darkenergy(self.w0, self.wa)

    def calc(self, ra, dec, z):
        """ This is where we do all the calculating. 
        I have the ra/dec/z as input, in the form of numpy arrays. 
        I return arrays for shear1, shear2 and conv. 
        """

        shear1 = numpy.empty(len(ra))
        shear2 = numpy.empty(len(ra))
        conv = numpy.empty(len(ra))

        # I'm gonna have to loop through these arrays, irritatingly.
        for g in range(0, len(ra)):
            gal_x, gal_y = self.get_pixel_coordinates(ra[g], dec[g])

            i = 0
            while self.source_redshift[i+1] < z[g] and i < self.number_of_maps-2:
                i += 1

            shear1map_close = self.shear1map[i]
            shear2map_close = self.shear2map[i]
            convmap_close = self.convmap[i]
            shear1map_far = self.shear1map[i+1]
            shear2map_far = self.shear2map[i+1]
            convmap_far = self.convmap[i+1]
            gal_z_close = self.source_redshift[i]
            gal_z_far = self.source_redshift[i+1]
            comoving_distance_close = self.comoving_distance[i]
            comoving_distance_far = self.comoving_distance[i+1]

    # interpolate shear from close and far shear maps

            shear1_close = get_interpolated_value(shear1map_close, self.NbinsX, self.NbinsY, gal_x, gal_y)
            shear2_close = get_interpolated_value(shear2map_close, self.NbinsX, self.NbinsY, gal_x, gal_y)
            conv_close = get_interpolated_value(convmap_close, self.NbinsX, self.NbinsY, gal_x, gal_y)
            shear1_far = get_interpolated_value(shear1map_far, self.NbinsX, self.NbinsY, gal_x, gal_y)
            shear2_far = get_interpolated_value(shear2map_far, self.NbinsX, self.NbinsY, gal_x, gal_y)
            conv_far = get_interpolated_value(convmap_far, self.NbinsX, self.NbinsY, gal_x, gal_y)

            # this is not used for anything. So leave it out for now.
            ##gal_comoving_distance = calculate_comoving_distance( 1.0/(1.0+z[g]), self.Omega_m, self.Omega_Lambda, self.w0, self.wa)

            shear1[g] = weight_shear_2(z[g], gal_z_close, gal_z_far, shear1_close, shear1_far)
            shear2[g] = weight_shear_2(z[g], gal_z_close, gal_z_far, shear2_close, shear2_far)
            conv[g] = weight_shear_2(z[g], gal_z_close, gal_z_far, conv_close, conv_far)

            del(shear1map_close)
            del(shear2map_close)
            del(convmap_close)
            del(shear1map_far)
            del(shear2map_far)
            del(convmap_far)

        #self.shear1map = []
        #self.shear2map = []
        #self.convmap = []
        return shear1, shear2, conv

    def get_filename(self, basename, z):
        """ This is unlikely to be used as originally intended, 
        as we won't be reading in more than one set of WL maps yet. 
        """

        redshift_tag = int(math.floor(100*self.source_redshift[z]))
        filename = "/astro/net/lsst1/shared/djbard/mirrored/WL-"+basename+"_m-512b240_Om0.260_Ol0.740_w-1.000_ns0.960_si0.798_4096xy_0001r_00" + \
            str(self.last_plane[z])+"p_0"+str(redshift_tag)+"z_og.gre_X4_smaller.fit"

        return filename

    def get_pixel_coordinates(self, ra, dec):
        """ Returns pixel coord for WL maps. Note that we repeat the 
        same WL map periodically across the sky - not ideal! 
        """
        MINSIZE = 1.0e-8
        if ra >= 12:
            ra -= 24.0
        xx = ra * (360.0/24.0) * (self.NbinsX/self.survey_angle) + (self.NbinsX/2.0)
        yy = dec * (self.NbinsY/self.survey_angle) + (self.NbinsX/2.0)

 # make it periodic....
        multiple = math.floor(xx/self.NbinsX)
        x = xx - multiple*self.NbinsX
        multiple = math.floor(yy/self.NbinsX)
        y = yy - multiple*self.NbinsX

        # this can happen if an object is *right* on the edge of a map.
        # i.e. within a pixel or two.
        # it's because the bicubic interpolation breaks at the edge of a map.
        # It's not really out of range of the map, though.
        # if x<(1.0-MINSIZE) or x>=(self.NbinsX - 2.0+MINSIZE):
        #    print "galaxy ra position excedes range of WL maps", ra, xx, x
        # if y<(1.0-MINSIZE) or y>=(self.NbinsY - 2.0+MINSIZE):
        #    print "galaxy dec position excedes range of WL maps", dec, yy, y

        if x < 1.0:
            x = 1.0
        if x >= self.NbinsX-2.0+MINSIZE:
            x = self.NbinsX-2.0+MINSIZE
        if y < 1.0:
            y = 1.0
        if y >= self.NbinsY-2.0+MINSIZE:
            y = self.NbinsY-2.0+MINSIZE

        return x, y



