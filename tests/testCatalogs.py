from __future__ import with_statement
import os
import numpy as np
import unittest
from lsst.utils import getPackageDir
import lsst.utils.tests as utilsTests
from lsst.sims.utils import ObservationMetaData
from lsst.sims.catalogs.generation.db import CatalogDBObject, fileDBObject
from lsst.sims.catalogs.measures.instance import InstanceCatalog, compound
from lsst.sims.utils import haversine, observedFromICRS

# a class of catalog that outputs all the significant figures in
# ra and dec so that it can be read back in to make sure that our
# Haversine-based query actually returns all of the points that
# are inside the circular bound desired


class BoundsCatalog(InstanceCatalog):
    catalog_type = 'bounds_catalog'
    refIdCol = 'id'
    column_outputs = ['id', 'raJ2000', 'decJ2000']

    default_formats = {'f': '%.20f'}


def twice_fn(x):
    return 2.0*x


class TransformationCatalog(InstanceCatalog):
    catalog_type = 'transformation_catalog'
    refIdCol = 'id'
    column_outputs = ['id', 'raJ2000', 'decJ2000']
    default_formats = {'f': '%.12f'}
    transformations = {'raJ2000': twice_fn}


class BasicCatalog(InstanceCatalog):
    catalog_type = 'basic_catalog'
    refIdCol = 'id'
    column_outputs = ['id', 'raJ2000', 'decJ2000', 'umag', 'gmag', 'rmag', 'imag',
                      'zmag', 'ymag']

    default_formats = {'f': '%.12f'}


class TestAstMixin(object):

    @compound('ra_corr', 'dec_corr')
    def get_points_corrected(self):
        # Fake astrometric correction
        ra_corr = self.column_by_name('raJ2000')+0.001
        dec_corr = self.column_by_name('decJ2000')+0.001
        return ra_corr, dec_corr


class CustomCatalog(BasicCatalog, TestAstMixin):
    catalog_type = 'custom_catalog'
    refIdCol = 'id'
    column_outputs = ['id', 'raJ2000', 'decJ2000', 'umag', 'gmag', 'rmag', 'imag',
                      'zmag', 'ymag', 'ra_corr', 'dec_corr']


def compareFiles(file1, file2):
    with open(file1) as fh:
        str1 = "".join(fh.readlines())
    with open(file2) as fh:
        str2 = "".join(fh.readlines())
    return str1 == str2


def write_star_file_db(file_name):

    np.random.seed(88)
    nstars = 10000
    ra = np.random.random_sample(nstars)*360.0
    dec = (np.random.random_sample(nstars)-0.5)*180.0
    umag = np.random.random_sample(nstars)*10.0 + 15.0
    gmag = np.random.random_sample(nstars)*10.0 + 15.0
    rmag = np.random.random_sample(nstars)*10.0 + 15.0
    imag = np.random.random_sample(nstars)*10.0 + 15.0
    zmag = np.random.random_sample(nstars)*10.0 + 15.0
    ymag = np.random.random_sample(nstars)*10.0 + 15.0

    with open(file_name, 'w') as output_file:
        for ix, (rr, dd, um, gm, rm, im, zm, ym) in \
                enumerate(zip(ra, dec, umag, gmag, rmag, imag, zmag, ymag)):

            output_file.write('%d %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f\n' %
                              (ix, rr, dd, um, gm, rm, im, zm, ym))

    starDtype = np.dtype([
        ('id', np.int),
        ('raJ2000', np.float),
        ('decJ2000', np.float),
        ('umag', np.float),
        ('gmag', np.float),
        ('rmag', np.float),
        ('imag', np.float),
        ('zmag', np.float),
        ('ymag', np.float)
    ])

    starDB = fileDBObject(file_name, runtable='stars', dtype=starDtype, idColKey='id')
    starDB.raColName = 'raJ2000'
    starDB.decColName = 'decJ2000'

    controlData = np.genfromtxt(file_name, dtype=starDtype)

    return starDB, controlData


class InstanceCatalogTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        cls.starTextName = os.path.join(getPackageDir('sims_catalogs_measures'), 'tests',
                                        'scratchSpace', 'icStarTestCatalog.txt')

        if os.path.exists(cls.starTextName):
            os.unlink(cls.starTextName)

        cls.starDB, cls.starControlData = write_star_file_db(cls.starTextName)

    @classmethod
    def tearDownClass(cls):

        if os.path.exists(cls.starTextName):
            os.unlink(cls.starTextName)

    def setUp(self):
        self.obsMd = ObservationMetaData(boundType = 'circle', pointingRA = 210.0, pointingDec = -60.0,
                                         boundLength=20.0, mjd=52000., bandpassName='r')

    def testStarLike(self):
        """
        Write a couple of catalogs.  Verify that the objects that end up in the catalog fall within the pointing
        and that the objects that do not end up in the catalog fall outside of the pointing
        """

        catName = os.path.join(getPackageDir('sims_catalogs_measures'),
                               'tests', 'scratchSpace', '_starLikeCat.txt')

        if os.path.exists(catName):
            os.unlink(catName)

        # this dtype corresponds to the outputs of the catalog
        dtype = np.dtype([
            ('id', np.int),
            ('raJ2000', np.float),
            ('decJ2000', np.float),
            ('umag', np.float),
            ('gmag', np.float),
            ('rmag', np.float),
            ('imag', np.float),
            ('zmag', np.float),
            ('ymag', np.float),
            ('ra_corr', np.float),
            ('dec_corr', np.float)
        ])

        t = self.starDB.getCatalog('custom_catalog', obs_metadata=self.obsMd)
        t.write_catalog(catName)

        testData = np.genfromtxt(catName, delimiter = ', ', dtype=dtype)

        # make sure that something ended up in the catalog
        self.assertGreater(len(testData), 0)

        # iterate over the lines in the catalog
        # verify that those line exist in the control data
        # also verify that those lines fall within the requested field of view
        for line in testData:
            ic = np.where(self.starControlData['id'] == line['id'])[0][0]
            self.assertAlmostEqual(line['umag'], self.starControlData['umag'][ic], 6)
            self.assertAlmostEqual(line['gmag'], self.starControlData['gmag'][ic], 6)
            self.assertAlmostEqual(line['rmag'], self.starControlData['rmag'][ic], 6)
            self.assertAlmostEqual(line['imag'], self.starControlData['imag'][ic], 6)
            self.assertAlmostEqual(line['zmag'], self.starControlData['zmag'][ic], 6)
            self.assertAlmostEqual(line['ymag'], self.starControlData['ymag'][ic], 6)
            self.assertAlmostEqual(line['raJ2000'], self.starControlData['raJ2000'][ic], 6)
            self.assertAlmostEqual(line['decJ2000'], self.starControlData['decJ2000'][ic], 6)
            self.assertAlmostEqual(line['ra_corr'], line['raJ2000']+0.001, 6)
            self.assertAlmostEqual(line['dec_corr'], line['decJ2000']+0.001, 6)
            dl = haversine(np.radians(line['raJ2000']), np.radians(line['decJ2000']),
                           self.obsMd._pointingRA, self.obsMd._pointingDec)

            self.assertLess(np.degrees(dl), self.obsMd.boundLength)

        # examine the lines that did not fall in the catalog
        lines_not_in_catalog = np.where(self.starControlData['id'] not in testData['id'])[0]

        self.assertGreater(len(lines_not_in_catalog), 0)

        # make sure that those lines are, indeed, outside of the field of view
        for ic in lines_not_in_catalog:
            dl = haversine(self.obsMd._pointingRA, self.obsMd._pointingDec,
                           np.radians(self.starControlData['raJ2000'][ic]),
                           np.radians(self.starControlData['decJ2000'][ic]))

            self.assertGreater(np.degrees(dl), self.obsMd.boundLength)

        if os.path.exists(catName):
            os.unlink(catName)

        # now do the same thing for the basic catalog class

        dtype = np.dtype([
            ('id', np.int),
            ('raJ2000', np.float),
            ('decJ2000', np.float),
            ('umag', np.float),
            ('gmag', np.float),
            ('rmag', np.float),
            ('imag', np.float),
            ('zmag', np.float),
            ('ymag', np.float)
        ])

        t = self.starDB.getCatalog('basic_catalog', obs_metadata=self.obsMd)
        t.write_catalog(catName)

        testData = np.genfromtxt(catName, delimiter = ', ', dtype=dtype)

        # make sure that something ended up in the catalog
        self.assertGreater(len(testData), 0)

        # iterate over the lines in the catalog
        # verify that those line exist in the control data
        # also verify that those lines fall within the requested field of view
        for line in testData:
            ic = np.where(self.starControlData['id'] == line['id'])[0][0]
            self.assertAlmostEqual(line['umag'], self.starControlData['umag'][ic], 6)
            self.assertAlmostEqual(line['gmag'], self.starControlData['gmag'][ic], 6)
            self.assertAlmostEqual(line['rmag'], self.starControlData['rmag'][ic], 6)
            self.assertAlmostEqual(line['imag'], self.starControlData['imag'][ic], 6)
            self.assertAlmostEqual(line['zmag'], self.starControlData['zmag'][ic], 6)
            self.assertAlmostEqual(line['ymag'], self.starControlData['ymag'][ic], 6)
            self.assertAlmostEqual(line['raJ2000'], self.starControlData['raJ2000'][ic], 6)
            self.assertAlmostEqual(line['decJ2000'], self.starControlData['decJ2000'][ic], 6)
            dl = haversine(np.radians(line['raJ2000']), np.radians(line['decJ2000']),
                           self.obsMd._pointingRA, self.obsMd._pointingDec)

            self.assertLess(np.degrees(dl), self.obsMd.boundLength)

        # examine the lines that did not fall in the catalog
        lines_not_in_catalog = np.where(self.starControlData['id'] not in testData['id'])[0]

        self.assertGreater(len(lines_not_in_catalog), 0)

        # make sure that those lines are, indeed, outside of the field of view
        for ic in lines_not_in_catalog:
            dl = haversine(self.obsMd._pointingRA, self.obsMd._pointingDec,
                           np.radians(self.starControlData['raJ2000'][ic]),
                           np.radians(self.starControlData['decJ2000'][ic]))

            self.assertGreater(np.degrees(dl), self.obsMd.boundLength)

        if os.path.exists(catName):
            os.unlink(catName)

    def test_transformation(self):
        """
        Test that transformations are applied to columns in an InstanceCatalog
        """
        catName = os.path.join(getPackageDir('sims_catalogs_measures'), 'tests', 'scratchSpace',
                               'transformation_catalog.txt')

        if os.path.exists(catName):
            os.unlink(catName)

        t = self.starDB.getCatalog('transformation_catalog', obs_metadata=self.obsMd)
        t.write_catalog(catName)

        dtype = np.dtype([
                         ('id', np.int),
                         ('raJ2000', np.float),
                         ('decJ2000', np.float)
                         ])

        testData = np.genfromtxt(catName, delimiter=', ', dtype=dtype)
        self.assertGreater(len(testData), 0)
        for line in testData:
            ic = np.where(self.starControlData['id'] == line['id'])[0][0]
            self.assertAlmostEqual(line['decJ2000'], self.starControlData['decJ2000'][ic], 5)
            self.assertAlmostEqual(line['raJ2000'], 2.0*self.starControlData['raJ2000'][ic], 5)

        if os.path.exists(catName):
            os.unlink(catName)


class boundingBoxTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.starTextName = os.path.join(getPackageDir('sims_catalogs_measures'), 'tests', 'scratchSpace',
                                        'bbStarTestCatalog.txt')

        cls.starDB, cls.starControlData = write_star_file_db(cls.starTextName)

    @classmethod
    def tearDownClass(cls):
        if os.path.exists(cls.starTextName):
            os.unlink(cls.starTextName)

    def setUp(self):

        self.RAmin = 190.
        self.RAmax = 210.
        self.DECmin = -70.
        self.DECmax = -50.

        self.RAcenter = 200.
        self.DECcenter = -60.
        self.radius = 10.0

        self.obsMdCirc = ObservationMetaData(boundType='circle', pointingRA=self.RAcenter, pointingDec=self.DECcenter,
                                             boundLength=self.radius, mjd=52000., bandpassName='r')

        self.obsMdBox = ObservationMetaData(boundType='box', pointingRA=0.5*(self.RAmax+self.RAmin),
                                            pointingDec=0.5*(self.DECmin+self.DECmax),
                                            boundLength=np.array(
                                                [0.5*(self.RAmax-self.RAmin), 0.5*(self.DECmax-self.DECmin)]),
                                            mjd=52000., bandpassName='r')

    def testBoxBounds(self):
        """
        Make sure that box_bound_constraint in sims.catalogs.generation.db.dbConnection.py
        does not admit any objects outside of the bounding box
        """

        catName = os.path.join(getPackageDir('sims_catalogs_measures'), 'tests',
                               'scratchSpace', 'box_test_catalog.txt')

        myCatalog = self.starDB.getCatalog('bounds_catalog', obs_metadata = self.obsMdBox)

        myIterator = myCatalog.iter_catalog(chunk_size=10)

        for line in myIterator:
            self.assertGreater(line[1], self.RAmin)
            self.assertLess(line[1], self.RAmax)
            self.assertGreater(line[2], self.DECmin)
            self.assertLess(line[2], self.DECmax)

        myCatalog.write_catalog(catName)

        # now we will test for the completeness of the box bounds

        dtype = np.dtype([('id', np.int),
                          ('raJ2000', np.float),
                          ('decJ2000', np.float)])

        testData = np.genfromtxt(catName, dtype=dtype, delimiter=', ')

        for line in testData:
            self.assertGreater(line['raJ2000'], self.RAmin)
            self.assertGreater(line['decJ2000'], self.DECmin)
            self.assertLess(line['raJ2000'], self.RAmax)
            self.assertLess(line['decJ2000'], self.DECmax)

        ct = 0
        for line in self.starControlData:
            if line['id'] not in testData['id']:
                ct += 1
                in_bounds = (line['raJ2000'] < self.RAmax) and (line['raJ2000'] > self.RAmin) \
                    and (line['decJ2000'] < self.DECmax) and (line['decJ2000'] > self.DECmin)

                msg = 'violates bounds\nRA: %e < %e <%e\nDec: %e < %e < %e\n' % \
                    (self.RAmin, line['raJ2000'], self.RAmax,
                     self.DECmin, line['decJ2000'], self.DECmax)

                self.assertFalse(in_bounds, msg=msg)

        self.assertGreater(ct, 0)

        if os.path.exists(catName):
            os.unlink(catName)

    def testCircBounds(self):
        """
        Make sure that circular_bound_constraint in sims.catalogs.generation.db.dbConnection.py
        does not admit any objects outside of the bounding circle
        """

        catName = os.path.join(getPackageDir('sims_catalogs_measures'), 'tests', 'scratchSpace',
                               'circular_test_catalog.txt')

        if os.path.exists(catName):
            os.unlink(catName)

        myCatalog = self.starDB.getCatalog('bounds_catalog', obs_metadata = self.obsMdCirc)
        myIterator = myCatalog.iter_catalog(chunk_size=10)

        for line in myIterator:
            rtest = np.degrees(haversine(np.radians(self.RAcenter), np.radians(self.DECcenter),
                                         np.radians(line[1]), np.radians(line[2])))

            self.assertLess(rtest, self.radius)

        myCatalog.write_catalog(catName)

        # now we will test for the completeness of the circular bounds

        dtype = np.dtype([('id', np.int),
                          ('raJ2000', np.float),
                          ('decJ2000', np.float)])

        testData = np.genfromtxt(catName, dtype=dtype, delimiter=', ')

        self.assertGreater(len(testData), 0)

        for line in testData:
            dl = np.degrees(haversine(np.radians(line['raJ2000']), np.radians(line['decJ2000']),
                                      np.radians(self.RAcenter), np.radians(self.DECcenter)))

            self.assertLess(dl, self.radius)

        ct = 0
        for line in self.starControlData:
            if line['id'] not in testData['id']:
                ct += 1
                dl = np.degrees(haversine(np.radians(line['raJ2000']), np.radians(line['decJ2000']),
                                          np.radians(self.RAcenter), np.radians(self.DECcenter)))

                self.assertGreater(dl, self.radius)

        self.assertGreater(ct, 0)

        if os.path.exists(catName):
            os.unlink(catName)


def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(InstanceCatalogTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    suites += unittest.makeSuite(boundingBoxTest)

    return unittest.TestSuite(suites)


def run(shouldExit=False):
    """Run the tests"""
    utilsTests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
