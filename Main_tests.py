# MIT License
#
# Copyright (c) 2023 Nicolette Shaw
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import os
import unittest

from Main import *

cwd = os.getcwd()
path_periodicTable = cwd + "/PeriodicInfo.csv"
path_totESCSdata = cwd + "/TotalESCS.csv"


class TestLoadData(unittest.TestCase):
    def test_get_totESCS_data(self):
        # Test with a known element, Carbon.
        element = "Carbon"
        expected_output = {1: 5.28E-01, 2: 2.86E-01, 4: 1.50E-01, 8: 7.75E-02, 16: 4.01E-02, 32: 2.11E-02, 64: 1.15E-02,
                           128: 6.71E-03, 256: 4.35E-03}
        self.assertEqual(loadData.get_totESCS_data(element, path_totESCSdata), expected_output)

        # Test with an element that does not exist in the data.
        element = 'X'
        expected_output = "Element X is not available in TotalESCS data."
        with self.assertRaises(SystemExit) as cm:
            loadData.get_totESCS_data(element, path_totESCSdata)
            self.assertEqual(cm.exception, expected_output)

        # Test with an empty string as the element.
        element = ''
        expected_output = "Element given to get_totESCS_data is empty."
        with self.assertRaises(SystemExit) as cm:
            loadData.get_totESCS_data(element, path_totESCSdata)
            self.assertEqual(cm.exception, expected_output)

    def test_get_element_data(self):
        # Test that the function correctly retrieves element data for a given element.
        elem = 'Carbon'
        expected_output = {'Z': 6, 'MW': 12.011}
        self.assertEqual(loadData.get_element_data(elem, path_periodicTable), expected_output)

        # Test with an element that does not exist in the data.
        element = 'X'
        expected_output = "Element X is not available in PeriodicInfo data."
        with self.assertRaises(SystemExit) as cm:
            loadData.get_totESCS_data(element, path_totESCSdata)
            self.assertEqual(cm.exception, expected_output)

        # Test with an empty string as the element.
        element = ''
        expected_output = "Element given to get_element_data is empty."
        with self.assertRaises(SystemExit) as cm:
            loadData.get_totESCS_data(element, path_totESCSdata)
            self.assertEqual(cm.exception, expected_output)


class TestMathHelper(unittest.TestCase):
    def test_TotESCStoTotSCS(self):
        # Test that the function correctly converts totESCS to totSCS.
        totESCS = 2.0
        Znum = 6
        expected_output = 8.0
        self.assertEqual(mathHelper.TotESCStoTotSCS(totESCS, Znum), expected_output)

        # Test that the function handles totESCS = 0 correctly.
        totESCS = 0
        Znum = 6
        expected_output = 0
        self.assertEqual(mathHelper.TotESCStoTotSCS(totESCS, Znum), expected_output)


class TestCompoundCalculationHelper(unittest.TestCase):
    def assertcompoundDictAlmostEqual(self, d1, d2):

        # First layer.
        # check if both inputs are dicts
        self.assertIsInstance(d1, dict, 'First argument is not a dictionary')
        self.assertIsInstance(d2, dict, 'Second argument is not a dictionary')

        # check if both inputs have the same keys
        self.assertEqual(d1.keys(), d2.keys())

        # Check that internal 'totSCS' and 'atomic' dicts are equal.
        for key, value in d1.items():
            # Check periodic info dict.
            self.assertEqual(d1[key]['atomic'], d2[key]['atomic'])
            # Check keys for d1[key]['totSCS'] (energies).
            self.assertEqual(d1[key]['totSCS'].keys(), d2[key]['totSCS'].keys())
            for i in d1[key]['totSCS'].keys():
                self.assertAlmostEqual(d1[key]['totSCS'][i], d2[key]['totSCS'][i], places=5)

    def test_converttotESCStoSCS(self):
        # Test for single element in compound.
        compoundDict = compoundInfo(name='Carbon', density=2.0, molecularFormulaDict={"Carbon": 1})
        compoundCalcHelper = compoundCalculationHelper(compoundDict)
        expected_output = {'Carbon': {
            'totSCS': {1: 2.112, 2: 1.144, 4: 0.6, 8: 0.31, 16: 0.1604, 32: 0.0844, 64: 0.046, 128: 0.02684,
                       256: 0.0174},
            'atomic': {'Z': 6, 'MW': 12.011}}}
        self.assertcompoundDictAlmostEqual(compoundCalcHelper.converttotESCStoSCS(), expected_output)

        # Test for multiple elements in compound.
        compoundDict = compoundInfo(name='Silicon nitride', density=3.44,
                                    molecularFormulaDict={"Silicon": 3, "Nitrogen": 4})
        compoundCalcHelper = compoundCalculationHelper(compoundDict)
        expected_output = {'Silicon': {
            'totSCS': {1: 3.08571428571429, 2: 1.86742857142857, 4: 1.07885714285714, 8: 0.596571428571429,
                       16: 0.322285714285714, 32: 0.173485714285714, 64: 0.0957714285714286, 128: 0.0562285714285714,
                       256: 0.0365714285714286},
            'atomic': {'Z': 14, 'MW': 28.0855}},
            'Nitrogen': {'totSCS': {1: 1.91071428571429, 2: 1.06071428571429, 4: 0.567857142857143,
                                    8: 0.295357142857143, 16: 0.153571428571429, 32: 0.0810714285714286,
                                    64: 0.0442857142857143, 128: 0.0258214285714286, 256: 0.01675},
                         'atomic': {'Z': 7, 'MW': 14.0067}}
        }
        self.assertcompoundDictAlmostEqual(compoundCalcHelper.converttotESCStoSCS(), expected_output)

    # def test_calcNumDensity(self):
    #     compoundDict = {'density': 1.2}
    #     compoundCalcHelper = compoundCalculationHelper(compoundDict)
    #
    #     # Test that the method correctly calculates number density
    #     elementMW = 12.0107
    #     expected_output = 4.988667836854447e+22
    #     self.assertEqual(compoundCalcHelper.calcNumDensity(elementMW), expected_output)


if __name__ == '__main__':
    unittest.main()
