import sys
import unittest
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
from scipy.interpolate import make_interp_spline, BSpline


# Path to Periodic table information.
path_periodicTable = ".../PeriodicInfo.csv"
# Path to total eleastic scattering cross section information for elements in preCalcElements.
# All cross sections are in units of angstroms squared - conversion to nm is done later.
path_totESCSdata = ".../TotalCS_SelectedElements.csv"

# A list of elements where pre-calculated total elastic scattering cross sections are available at set energies.
preCalcElements = ["Carbon", "Nitrogen", "Oxygen"]


def checkSameLists(lofl: dict, msg: str):
    res = True
    # extracting value to compare
    test_val = lofl[0]
    for l in lofl:
        if l != test_val:
            sys.exit(msg)


class loadData:

    @staticmethod
    def loadGivenData(elem: str) -> tuple[list[float], list[int]]:
        """
        :param path: Path to .csv file with total elastic scattering cross sections (in Angstroms squared) where elements are rows and energies are columns (in degrees).
        :param elem: The full name of the element.
        :return: A list of total elastic scattering cross sections and corresponding list of the energies.
        """
        data = pd.read_csv(path_totESCSdata)
        df = pd.DataFrame(data)
        # Retrieves list of energies (in keV) from data table.
        colnames = list(map(int, df.columns[1:]))
        rownum = df[df["Element"] == elem].index[0]
        # Retrieves list of total elastic scattering cross sections (in degrees) from data table.
        elemList = list(df.iloc[rownum][1:])
        return elemList, colnames

    @staticmethod
    def getElementData(elem: str,
                       path: str = path_periodicTable):
        data = pd.read_csv(path)
        df = pd.DataFrame(data)
        # Retrieves list of energies (in keV) from data table.
        rownum = df[df["Element"] == elem].index[0]
        # Retrieves list of total elastic scattering cross sections (in degrees) from data table.
        data = list(df.iloc[rownum][1:])
        return data


class test_loadData(unittest.TestCase):
    def test_loadGivenData(self):
        get1, get2 = loadData.loadGivenData(path_totESCSdata, "Oxygen")
        want1 = [5.27E-01, 3.01E-01, 1.63E-01, 8.62E-02, 4.50E-02, 2.38E-02, 1.30E-02, 7.62E-03,
                 4.94E-03]  # total elastic scattering cross sections (in degrees) of Oxygen.
        want2 = [1, 2, 4, 8, 16, 32, 64, 128, 256]  # energy list in keV
        self.assertEqual(get1, want1)
        self.assertEqual(get2, want2)


class mathHelper:

    @staticmethod
    def TotESCStoTotSCS(totESCS, Znum):
        return totESCS + (18 * totESCS / Znum)

    @staticmethod
    def TotESCSLISTtoTotSCSLIST(totESCSList, Znum):
        new = []
        for i in totESCSList:
            new.append(mathHelper.TotESCStoTotSCS(i, Znum))
        return new


class plotHelper:

    @staticmethod
    def smoothData(xlist: list[float], ylist: list[float], totalpoints: bool = 300):
        # Smooth data.
        xlist_smooth = np.linspace(min(xlist), max(xlist), totalpoints)
        spl = make_interp_spline(xlist, ylist, k=2)  # type: BSpline
        ylist_smooth = spl(xlist_smooth)
        return xlist_smooth, ylist_smooth


class totESCS:

    def __init__(self, element: str):
        """
        :param element: Full name of element.
        :param importFormat: If data is in ESCS or total ESCS format.
        """
        self.element = element

    def getTotESCS(self):
        """
        Given the element and if it is a fitted (totESCS) or to be calculated (ESCS) values, returns the total ESCS for each energy.
        :return: a tuple of [total ESCS values], [corresponding energies].
        """
        # If the element given is not in preCalcElements, TotESCS must be calculated from ESCS.
        totESCSList, energies = loadData.loadGivenData(self.element)
        return totESCSList, energies


class MFP:

    def __init__(self, compoundList):
        self.densityCompound = compoundList[0]
        self.compoundDict = compoundList[1]
        self.Avogadro = 6.02E+23

    def calcNumDensity(self, elementMW):
        return self.Avogadro * self.densityCompound / elementMW

    def calcMFP(self, totSCSDict_1energy):
        denom = 0
        stoichSum = sum(self.compoundDict.values())
        for elem in self.compoundDict.keys():
            elementMW = loadData.getElementData(elem)[1]
            add = self.calcNumDensity(elementMW) * totSCSDict_1energy[elem] * self.compoundDict[elem] / stoichSum
            denom = denom + add
        return 1 / denom * 10000000

class Transmission:

    def __init__(self, compoundList):
        self.compoundList = compoundList
        self.compoundDict = compoundList[1]
        self.Avogadro = 6.02E+23

        energiesLists = []
        totSCSDict = {}
        for element in self.compoundDict:
            totESCSList, energies = totESCS(element=element).getTotESCS()
            # Translate totESCS in Angstroms^2 to totSCS in cm^2.
            totSCSList = mathHelper.TotESCSLISTtoTotSCSLIST(totESCSList, loadData.getElementData(element)[0])
            energiesLists.append(energies)
            totSCSDict[element] = totSCSList

        # Check that energies across elements are consistent.
        checkSameLists(energiesLists, "Elements in compound do not have matching input energies.")
        self.energies = energiesLists[0]
        self.totSCSDict = totSCSDict

    def calcTransmission(self, thickness, energy):
        # Get the energy in list closest to the energy given.
        if energy not in self.energies:
            energy = min(self.energies, key=lambda x: abs(x - energy))
            print("Used closest energy available for calculation: " + str(energy))

        # Switch to integer for energy in list.
        energyInt = self.energies.index(energy)

        # Parse to one-energy dictionary of totSCS.
        totSCSDict_1energy = {}
        for elem in self.totSCSDict:
            totSCSDict_1energy[elem] = self.totSCSDict[elem][energyInt]
        mfp = MFP(self.compoundList).calcMFP(totSCSDict_1energy) * math.pow(10, 16)  # convert to nm
        return math.e ** (- thickness / mfp)


def plotTransmission(compoundList, thickness, compoundName, linetype, energies=[1, 2, 4, 8, 16, 32, 64, 128, 256]):
    y = []
    for e in energies:
        transmission = Transmission(compoundList).calcTransmission(thickness, e)
        y.append(transmission)
        if e == 128:
            print(compoundName + " " + str(thickness) + " " + str(round(transmission * 100,2)))
    energies, y = plotHelper.smoothData(energies, y)
    plt.plot(energies, y, linetype, label=str(thickness) + " nm " + str(compoundName))

# format: [density of compound, {elem: stoichiometric ratio, elem: stoichiometric ratio}]
SiN = [3.44, {"Silicon": 3, "Nitrogen": 4}]
Carbon = [2.0, {"Carbon": 1}]
Graphene = [2.267, {"Carbon": 1}]

## To recreate Fig 8. from:
# Dwyer, J. R.; Harb, M.
# Through a Window, Brightly: A Review of Selected Nanofabricated Thin-Film Platforms for Spectroscopy, Imaging, and Detection.
# Appl. Spectrosc. 2017, 71 (9), 2051–2075.

for i in [5, 10, 20, 40]:
    plotTransmission(SiN, i, "SiN", "k-")
for i in [5, 40]:
    plotTransmission(Carbon, i, "Carbon", "g-.")
for i in [0.335]:
   plotTransmission(Graphene, i, "Graphene bilayer", "b:")


plt.legend()
plt.ylabel("Transmission")
plt.xlabel("Energy (keV)")
plt.show()

if __name__ == '__main__':
    unittest.main()
