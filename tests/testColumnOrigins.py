import os
import numpy
import unittest
import sqlite3
import json
import lsst.utils.tests as utilsTests
from lsst.sims.catalogs.measures.instance import InstanceCatalog, cached, compound
from lsst.sims.catalogs.generation.db import CatalogDBObject


def makeTestDB(name, size=10, **kwargs):
    """
    Make a test database

    @param [in] name is a string indicating the name of the database file
    to be created

    @param [in] size is an int indicating the number of objects to include
    in the database (default=10)
    """
    conn = sqlite3.connect(name)
    c = conn.cursor()
#    try:
    c.execute('''CREATE TABLE testTable
                 (id int, aa float, bb float, ra float, decl float)''')
    conn.commit()
#    except Exception:
#        raise RuntimeError("Error creating database.")

    for i in xrange(size):

        ra = numpy.random.sample()*360.0
        dec = (numpy.random.sample()-0.5)*180.0

        # insert the row into the data base
        qstr = '''INSERT INTO testTable VALUES (%i, '%f', '%f', '%f', '%f')''' % (i, 2.0*i, 3.0*i, ra, dec)
        c.execute(qstr)

    conn.commit()
    conn.close()


class testDBObject(CatalogDBObject):
    objid = 'testDBObject'
    tableid = 'testTable'
    idColKey = 'id'
    # Make this implausibly large?
    appendint = 1023
    database = 'colOriginsTestDatabase.db'
    driver = 'sqlite'
    raColName = 'ra'
    decColName = 'decl'
    columns = [('objid', 'id', int),
               ('raJ2000', 'ra*%f'%(numpy.pi/180.)),
               ('decJ2000', 'decl*%f'%(numpy.pi/180.)),
               ('aa', None),
               ('bb', None)]

# Below we define mixins which calculate the variables 'cc' and 'dd in different
# ways.  The idea is to see if InstanceCatalog correctly identifies where
# the columns come from in those cases


class mixin1(object):

    @cached
    def get_cc(self):
        aa = self.column_by_name('aa')
        bb = self.column_by_name('bb')

        return numpy.array(aa-bb)

    @cached
    def get_dd(self):
        aa = self.column_by_name('aa')
        bb = self.column_by_name('bb')

        return numpy.array(aa+bb)


class mixin2(object):

    @compound('cc', 'dd')
    def get_both(self):
        aa = self.column_by_name('aa')
        bb = self.column_by_name('bb')

        return numpy.array([aa-bb, aa+bb])


class mixin3(object):

    @cached
    def get_cc(self):
        aa = self.column_by_name('aa')
        bb = self.column_by_name('bb')

        return numpy.array(aa-bb)

# Below we define catalog classes that use different combinations
# of the mixins above to calculate the columns 'cc' and 'dd'


class testCatalogDefaults(InstanceCatalog):
    column_outputs = ['objid', 'aa', 'bb', 'cc', 'dd', 'raJ2000', 'decJ2000']
    default_columns = [('cc', 0.0, float), ('dd', 1.0, float)]


class testCatalogMixin1(InstanceCatalog, mixin1):
    column_outputs = ['objid', 'aa', 'bb', 'cc', 'dd', 'raJ2000', 'decJ2000']
    default_columns = [('cc', 0.0, float), ('dd', 1.0, float)]


class testCatalogMixin2(InstanceCatalog, mixin2):
    column_outputs = ['objid', 'aa', 'bb', 'cc', 'dd', 'raJ2000', 'decJ2000']
    default_columns = [('cc', 0.0, float), ('dd', 1.0, float)]


class testCatalogMixin3(InstanceCatalog, mixin3):
    column_outputs = ['objid', 'aa', 'bb', 'cc', 'dd', 'raJ2000', 'decJ2000']
    default_columns = [('cc', 0.0, float), ('dd', 1.0, float)]


class testCatalogMixin3Mixin1(InstanceCatalog, mixin3, mixin1):
    column_outputs = ['objid', 'aa', 'bb', 'cc', 'dd', 'raJ2000', 'decJ2000']
    default_columns = [('cc', 0.0, float), ('dd', 1.0, float)]


class testCatalogAunspecified(InstanceCatalog, mixin3, mixin1):
    column_outputs = ['objid', 'bb', 'cc', 'dd', 'raJ2000', 'decJ2000']
    default_columns = [('aa', -1.0, float), ('cc', 0.0, float), ('dd', 1.0, float)]


class testColumnOrigins(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.dbName = 'colOriginsTestDatabase.db'
        if os.path.exists(cls.dbName):
            os.unlink(cls.dbName)

        makeTestDB(cls.dbName)

    @classmethod
    def tearDownClass(cls):
        if os.path.exists(cls.dbName):
            os.unlink(cls.dbName)

    def setUp(self):
        self.myDBobject = testDBObject(database=self.dbName)
        self.mixin1Name = '<class \'__main__.mixin1\'>'
        self.mixin2Name = '<class \'__main__.mixin2\'>'
        self.mixin3Name = '<class \'__main__.mixin3\'>'

    def tearDown(self):
        del self.myDBobject
        del self.mixin1Name
        del self.mixin2Name
        del self.mixin3Name

    def testDefaults(self):
        """
        Test case where the columns cc and dd come from defaults
        """
        myCatalog = testCatalogDefaults(self.myDBobject)

        self.assertEqual(myCatalog._column_origins['objid'], 'the database')
        self.assertEqual(myCatalog._column_origins['raJ2000'], 'the database')
        self.assertEqual(myCatalog._column_origins['decJ2000'], 'the database')
        self.assertEqual(myCatalog._column_origins['aa'], 'the database')
        self.assertEqual(myCatalog._column_origins['bb'], 'the database')
        self.assertEqual(myCatalog._column_origins['cc'], 'default column')
        self.assertEqual(myCatalog._column_origins['dd'], 'default column')

    def testMixin1(self):
        """
        Test case where the columns cc and dd come from non-compound getters
        """
        myCatalog = testCatalogMixin1(self.myDBobject)

        self.assertEqual(myCatalog._column_origins['objid'], 'the database')
        self.assertEqual(myCatalog._column_origins['raJ2000'], 'the database')
        self.assertEqual(myCatalog._column_origins['decJ2000'], 'the database')
        self.assertEqual(myCatalog._column_origins['aa'], 'the database')
        self.assertEqual(myCatalog._column_origins['bb'], 'the database')
        self.assertEqual(str(myCatalog._column_origins['cc']), self.mixin1Name)
        self.assertEqual(str(myCatalog._column_origins['dd']), self.mixin1Name)

    def testMixin2(self):
        """
        Test case where the columns cc and dd come from a compound getter
        """
        myCatalog = testCatalogMixin2(self.myDBobject)

        self.assertEqual(myCatalog._column_origins['objid'], 'the database')
        self.assertEqual(myCatalog._column_origins['raJ2000'], 'the database')
        self.assertEqual(myCatalog._column_origins['decJ2000'], 'the database')
        self.assertEqual(myCatalog._column_origins['aa'], 'the database')
        self.assertEqual(myCatalog._column_origins['bb'], 'the database')
        self.assertEqual(str(myCatalog._column_origins['cc']), self.mixin2Name)
        self.assertEqual(str(myCatalog._column_origins['dd']), self.mixin2Name)

    def testMixin3(self):
        """
        Test case where cc comes from a mixin and dd comes from the default
        """
        myCatalog = testCatalogMixin3(self.myDBobject)

        self.assertEqual(myCatalog._column_origins['objid'], 'the database')
        self.assertEqual(myCatalog._column_origins['raJ2000'], 'the database')
        self.assertEqual(myCatalog._column_origins['decJ2000'], 'the database')
        self.assertEqual(myCatalog._column_origins['aa'], 'the database')
        self.assertEqual(myCatalog._column_origins['bb'], 'the database')
        self.assertEqual(str(myCatalog._column_origins['cc']), self.mixin3Name)
        self.assertEqual(str(myCatalog._column_origins['dd']), 'default column')

    def testMixin3Mixin1(self):
        """
        Test case where one mixin overwrites another for calculating cc
        """
        myCatalog = testCatalogMixin3Mixin1(self.myDBobject)

        self.assertEqual(myCatalog._column_origins['objid'], 'the database')
        self.assertEqual(myCatalog._column_origins['raJ2000'], 'the database')
        self.assertEqual(myCatalog._column_origins['decJ2000'], 'the database')
        self.assertEqual(myCatalog._column_origins['aa'], 'the database')
        self.assertEqual(myCatalog._column_origins['bb'], 'the database')
        self.assertEqual(str(myCatalog._column_origins['cc']), self.mixin3Name)
        self.assertEqual(str(myCatalog._column_origins['dd']), self.mixin1Name)

    def testAunspecified(self):
        """
        Test case where aa is not specified in the catalog (and has a default)
        """
        myCatalog = testCatalogAunspecified(self.myDBobject)

        self.assertEqual(myCatalog._column_origins['objid'], 'the database')
        self.assertEqual(myCatalog._column_origins['raJ2000'], 'the database')
        self.assertEqual(myCatalog._column_origins['decJ2000'], 'the database')
        self.assertEqual(myCatalog._column_origins['aa'], 'the database')
        self.assertEqual(myCatalog._column_origins['bb'], 'the database')
        self.assertEqual(str(myCatalog._column_origins['cc']), self.mixin3Name)
        self.assertEqual(str(myCatalog._column_origins['dd']), self.mixin1Name)


class myDummyCatalogClass(InstanceCatalog):

    default_columns = [('sillyDefault', 2.0, float)]

    def get_cc(self):
        return self.column_by_name('aa')+1.0

    @compound('dd', 'ee', 'ff')
    def get_compound(self):

        return numpy.array([
                           self.column_by_name('aa')+2.0,
                           self.column_by_name('aa')+3.0,
                           self.column_by_name('aa')+4.0
                           ])


class myDependentColumnsClass_shouldPass(InstanceCatalog):

    def get_dd(self):

        if 'ee' in self._all_available_columns:
            delta = self.column_by_name('ee')
        else:
            delta = self.column_by_name('bb')

        return self.column_by_name('aa') + delta


class myDependentColumnsClass_shouldFail(InstanceCatalog):

    def get_cc(self):
        return self.column_by_name('aa')+1.0

    def get_dd(self):

        if 'ee' in self._all_available_columns:
            delta = self.column_by_name('ee')
        else:
            delta = self.column_by_name('bb')

        return self.column_by_name('aa') + delta

    def get_ee(self):
        return self.column_by_name('aa')+self.column_by_name('doesNotExist')


class AllAvailableColumns(unittest.TestCase):
    """
    This will contain a unit test to verify that the InstanceCatalog class
    self._all_available_columns contains all of the information it should
    """

    @classmethod
    def setUpClass(cls):
        cls.dbName = 'allGettersTestDatabase.db'
        if os.path.exists(cls.dbName):
            os.unlink(cls.dbName)

        makeTestDB(cls.dbName)

    @classmethod
    def tearDownClass(cls):
        if os.path.exists(cls.dbName):
            os.unlink(cls.dbName)

    def setUp(self):
        self.db = testDBObject(database=self.dbName)

    def testAllGetters(self):
        """
        test that the self._all_available_columns list contains all of the columns
        definedin an InstanceCatalog and its CatalogDBObject
        """
        cat = myDummyCatalogClass(self.db, column_outputs=['aa'])
        self.assertTrue('cc' in cat._all_available_columns)
        self.assertTrue('dd' in cat._all_available_columns)
        self.assertTrue('ee' in cat._all_available_columns)
        self.assertTrue('ff' in cat._all_available_columns)
        self.assertTrue('compound' in cat._all_available_columns)
        self.assertTrue('id' in cat._all_available_columns)
        self.assertTrue('aa' in cat._all_available_columns)
        self.assertTrue('bb' in cat._all_available_columns)
        self.assertTrue('ra' in cat._all_available_columns)
        self.assertTrue('decl' in cat._all_available_columns)
        self.assertTrue('decJ2000' in cat._all_available_columns)
        self.assertTrue('raJ2000' in cat._all_available_columns)
        self.assertTrue('objid' in cat._all_available_columns)
        self.assertTrue('sillyDefault' in cat._all_available_columns)

    def testDependentColumns(self):
        """
        We want to be able to use self._all_available_columns to change the calculation
        of columns on the fly (i.e. if a column exists, then use it to calculate
        another column; if it does not, ignore it).  This method tests whether
        or not that scheme will work.

        I have written two classes of catalogs.  The getter for the column 'dd'
        depends on the column 'doesNotExist', but only if the column 'ee' is defined.
        The class myDependentColumnsClass_shouldPass does not define a getter for
        'ee', so it does not require 'doesNotExist', so the constructor should pass.
        The class myDependentColumnsClass_shouldFail does have a getter for 'ee',
        so any catalog that requests the column 'dd' should fail to construct.
        """

        cat = myDependentColumnsClass_shouldPass(self.db, column_outputs=['dd'])

        # as long as we do not request the column 'dd', this should work
        cat = myDependentColumnsClass_shouldFail(self.db, column_outputs=['cc'])

        # because we are requesting the column 'dd', which depends on the fictitious column
        #'doesNotExist', this should raise an exception
        self.assertRaises(ValueError, myDependentColumnsClass_shouldFail, self.db, column_outputs=['dd'])


def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(testColumnOrigins)
    suites += unittest.makeSuite(AllAvailableColumns)
    return unittest.TestSuite(suites)


def run(shouldExit = False):
    utilsTests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
