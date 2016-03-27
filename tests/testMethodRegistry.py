from __future__ import with_statement
import unittest
import lsst.utils.tests as utilsTests
from lsst.sims.catalogs.measures.instance import register_class, register_method


@register_class
class ClassA(object):

    def call(self, key):
        return self._methodRegistry[key](self)

    @register_method('a')
    def _a_method(self):
        return 'a'


@register_class
class ClassB(ClassA):

    @register_method('b')
    def _b_method(self):
        return 'b'


@register_class
class ClassC(ClassB):

    @register_method('c')
    def _c_method(self):
        return 'c'


@register_class
class ClassD(ClassA):

    @register_method('d')
    def _d_method(self):
        return 'd'


class MethodRegistryTestCase(unittest.TestCase):

    def testMethodInheritance(self):
        """
        Test that the register_class and register_method decorators
        behave appropriately and preserve inheritance.
        """

        aa = ClassA()
        self.assertTrue(aa.call('a') == 'a')

        # below, we test to make sure that methods which
        # should not be in ClassA's _methodRegistry are not
        # spuriously added to the registry
        self.assertRaises(KeyError, aa.call, 'b')
        self.assertRaises(KeyError, aa.call, 'c')
        self.assertRaises(KeyError, aa.call, 'd')

        bb = ClassB()
        self.assertTrue(bb.call('a') == 'a')
        self.assertTrue(bb.call('b') == 'b')
        self.assertRaises(KeyError, bb.call, 'c')
        self.assertRaises(KeyError, bb.call, 'd')

        cc = ClassC()
        self.assertTrue(cc.call('a') == 'a')
        self.assertTrue(cc.call('b') == 'b')
        self.assertTrue(cc.call('c') == 'c')
        self.assertRaises(KeyError, cc.call, 'd')

        dd = ClassD()
        self.assertTrue(dd.call('a') == 'a')
        self.assertTrue(dd.call('d') == 'd')
        self.assertRaises(KeyError, dd.call, 'b')
        self.assertRaises(KeyError, dd.call, 'c')


def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(MethodRegistryTestCase)

    return unittest.TestSuite(suites)


def run(shouldExit=False):
    """Run the tests"""
    utilsTests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
