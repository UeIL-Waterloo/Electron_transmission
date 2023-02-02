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

import math
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import make_interp_spline
from scipy.optimize import curve_fit

# Path to Periodic table information.
path_periodicTable = "PeriodicInfo.csv"

# Path to total elastic scattering cross-section (totESCS) information for elements in preCalcElements.
# All cross-sections are in units of angstroms squared - conversion to nm is performed in calcTransmission.
path_totESCSdata = "TotalESCS.csv"


# Create compoundDict elements for compounds of interest.
class compoundInfo:
    def __init__(self, name, density, molecularFormulaDict, MW=None):
        self.name = name
        self.density = density
        self.molecularFormulaDict = molecularFormulaDict
        self.MW = MW


class loadData:

    @staticmethod
    def get_totESCS_data(element: str, path: str = path_totESCSdata) -> dict:
        # Read totESCS data from the CSV file.
        data = pd.read_csv(path)

        # Convert data to a Pandas DataFrame.
        df = pd.DataFrame(data)

        # Check that element is found in available data.
        if not element in list(df['Element']):
            sys.exit("Element " + str(element) + " is not available in TotalESCS data.")

        # Check that element is not empty.
        if element == '':
            sys.exit("Element given to get_totESCS_data is empty.")

        # Retrieve the list of energies (in keV) from the data table.
        colnames = list(map(int, df.columns[1:]))

        # Find the row in the data table for the given element.
        elemData = df[df["Element"] == element]

        # Retrieve the dictionary of total elastic scattering cross sections (in degrees) for each energy provided {energy: totESCS}.
        elemTotESSCSDict = elemData.to_dict('records')[0]

        # Remove the first element of the dict (element name).
        del elemTotESSCSDict["Element"]

        # Convert energies (keys) to floats.
        elemTotESSCSDict = {float(energy): elem for energy, elem in elemTotESSCSDict.items()}

        # Output has format {energy: totESCS, ...}, where totESCS is in Angstroms^2.
        return elemTotESSCSDict

    @staticmethod
    def get_element_data(element: str, path: str = path_periodicTable):
        # Read totESCS data from the CSV file.
        data = pd.read_csv(path)

        # Convert data to a Pandas DataFrame.
        df = pd.DataFrame(data)

        # Check that element is found in available data.
        if not element in list(df["Element"]):
            sys.exit("Element " + str(element) + " is not available in PeriodicInfo data.")

        # Check that element is not empty.
        if element == '':
            sys.exit("Element given to get_element_data is empty.")

        # Get row containing data for element.
        rownum = df[df["Element"] == element].index[0]

        # Get row data for element.
        data = list(df.iloc[rownum][1:])

        return {'Z': data[0], 'MW': data[1]}


class mathHelper:

    @staticmethod
    def TotESCStoTotSCS(totESCS: float, Znum: int):
        return totESCS + (18 * totESCS / Znum)

    @staticmethod
    def expFunc(x, a, b, c):
        return a * (1 - math.e ** (b * x + c))


class plotHelper:

    @staticmethod
    def smoothData(xlist: list[float], ylist: list[float], totalpoints: bool = 300):
        # Smooth data.
        xlist_smooth = np.linspace(min(xlist), max(xlist), totalpoints)
        spl = make_interp_spline(xlist, ylist, k=2)  # type: BSpline
        ylist_smooth = spl(xlist_smooth)
        return xlist_smooth, ylist_smooth


class compoundCalculationHelper:

    def __init__(self, compoundDict):
        self.compoundDict = compoundDict
        self.Avogadro = 6.02E+23

    # Convert totESCS dict of each element in compound to totSCS.
    def converttotESCStoSCS(self):
        compoundTotSCSDict = {}
        for element, ratio in self.compoundDict.molecularFormulaDict.items():
            # Get totESCS data for each element.
            elemTotESCSDict = loadData.get_totESCS_data(element=element)

            # Get Z (atomic number) for element.
            Z_MWDict = loadData.get_element_data(element)

            # Convert totESCS to totSCS.
            elemTotSCSDict = {energy: mathHelper.TotESCStoTotSCS(totESCS, Z_MWDict['Z']) for energy, totESCS in
                              elemTotESCSDict.items()}

            # Add two element dictionaries with totSCS data to a compound dictionary with the format {element: {'totSCS': elemTotSCSDict, 'atomic':Z_MWDict}}.
            compoundTotSCSDict[element] = {'totSCS': elemTotSCSDict, 'atomic': Z_MWDict}

        return compoundTotSCSDict

    def calcCompoundNumDensity(self, compoundMW):
        return (self.Avogadro * self.compoundDict.density) / compoundMW

    def calcMFPdenominator(self, compoundMW, ratioNum, totSCS_1energy):
        return self.calcCompoundNumDensity(compoundMW) * ratioNum * totSCS_1energy


class Transmission:

    def __init__(self, compoundDict):
        self.compoundDict = compoundDict
        self.compoundTotSCSDict = compoundCalculationHelper(compoundDict).converttotESCStoSCS()
        self.Avogadro = 6.02E+23
        if compoundDict.MW != None:
            self.MWtot = compoundDict.MW
        else:
            MWtot = 0
            for elem, elemTotSCSDict in self.compoundTotSCSDict.items():
                MW = elemTotSCSDict['atomic']['MW']
                # Add to MFP denominator value.
                ratioNum = self.compoundDict.molecularFormulaDict[elem]
                MWtot = MWtot + (MW * ratioNum)
            self.MWtot = MWtot

    def calcCompoundMFP(self, energy):
        # Get the energy in list closest to the energy given.
        # if energy not in self.energies:
        #     energy = min(self.energies, key=lambda x: abs(x - energy))
        #     print("Used closest energy available for calculation: " + str(energy))
        # TODO: add interpolation function for any energy calculations.

        denom = 0
        for elem, elemTotSCSDict in self.compoundTotSCSDict.items():
            totSCS = elemTotSCSDict['totSCS'][energy]
            # Add to MFP denominator value.
            ratioNum = self.compoundDict.molecularFormulaDict[elem]
            add = compoundCalculationHelper(self.compoundDict).calcMFPdenominator(compoundMW=self.MWtot,
                                                                                  ratioNum=ratioNum,
                                                                                  totSCS_1energy=totSCS)
            denom = denom + add
        return 1 / denom

    def calcMFP(self, energy):
        # Converts to nm.
        return self.calcCompoundMFP(energy) * math.pow(10, 23)

    def calcTransmission(self, thickness, energy):
        # Get MFP in nm.
        mfp = self.calcMFP(energy)
        # Calculate transmission value (fraction).
        return math.e ** (- (thickness / mfp))

    def getmfpVSEnergy(self, energies):
        transmissionClass = Transmission(self.compoundDict)
        mfp_data = []
        for energy in energies:
            mfp = transmissionClass.calcMFP(energy=energy)
            mfp_data.append(mfp)
        return mfp_data


class analyseTransmission:
    def __init__(self, compoundDict, thickness):
        self.compoundDict = compoundDict
        self.thickness = thickness
        self.energy_data = [1, 2, 4, 8, 16, 32, 64, 128, 256]

    def getTransmissionVSEnergy(self):
        transmissionClass = Transmission(self.compoundDict)
        transmission_data = []
        for energy in self.energy_data:
            transmission = transmissionClass.calcTransmission(self.thickness, energy=energy)
            transmission_data.append(transmission * 100)
        return transmission_data

    def plotTransmissionVSEnergy(self, smooth=True, color='red', linestyle='-'):
        y = analyseTransmission(self.compoundDict, self.thickness).getTransmissionVSEnergy()
        if smooth:
            e, y = plotHelper.smoothData(xlist=self.energy_data, ylist=y)

        # Create label.
        label = str(self.thickness) + 'nm ' + str(self.compoundDict.name)
        plt.plot(e, y, label=label, color=color, linestyle=linestyle)


if __name__ == '__main__':
    energies = [1, 2, 4, 8, 16, 32, 64, 128, 256]
    SiN = compoundInfo(name='Silicon nitride', density=3.2, molecularFormulaDict={"Silicon": 3, "Nitrogen": 4})
    Carbon = compoundInfo(name='Carbon', density=2.0, molecularFormulaDict={"Carbon": 1})
    Graphene = compoundInfo(name='Graphene', density=2.0, molecularFormulaDict={"Carbon": 1})
    Water = compoundInfo(name='Liquid water (oxygen)', density=0.997, molecularFormulaDict={"Oxygen": 1}, MW=18.01528)
    Water2 = compoundInfo(name='Liquid water (He2O)', density=0.997, molecularFormulaDict={"Helium": 2, "Oxygen": 1},
                          MW=18.01528)

    # Plot Liquid Cell transmission data.
    plt.figure()
    analyseTransmission(Graphene, thickness=0.335 * 2).plotTransmissionVSEnergy(color='black')
    analyseTransmission(Carbon, thickness=2.5 * 2).plotTransmissionVSEnergy(color='green')
    analyseTransmission(SiN, thickness=60).plotTransmissionVSEnergy(color='orange')
    analyseTransmission(SiN, thickness=10).plotTransmissionVSEnergy(color='orange', linestyle='--')
    plt.ylabel('Transmission (%)')
    plt.xlabel('Energy (keV)')
    plt.legend()
    plt.show()

    # Plot MFP for different compounds.
    plt.figure()
    for compound in [SiN, Carbon, Graphene]:
        mfpVals = Transmission(compound).getmfpVSEnergy(energies=energies)
        plt.plot(energies, mfpVals, 'o', label=str(compound.name))

    WaterMFP = Transmission(Water).getmfpVSEnergy(energies=energies)
    Water2MFP = Transmission(Water2).getmfpVSEnergy(energies=energies)
    WaterAvg = []
    for i in range(0, len(WaterMFP)):
        WaterAvg.append((WaterMFP[i] + Water2MFP[i]) / 2)

    plt.plot(energies, WaterMFP, 'o', label=str(Water.name))
    plt.plot(energies, Water2MFP, 'o', label=str(Water2.name))
    plt.plot(energies, WaterAvg, 'o', label='Water (average)')

    plt.ylabel('Total MFP (nm)')
    plt.xlabel('Energy (keV)')
    plt.ylim(0, 250)
    plt.legend()
    plt.show()

    # Plot comparison of water thickness to equivalent SiN transmission.
    plt.figure()

    equivalentThicknessWater = []
    thicknessSiN = 10  # nm
    for energy in energies:
        SiNMFP = Transmission(SiN).calcMFP(energy)
        Water2MFP = Transmission(Water2).calcMFP(energy)
        equivalentThicknessWater.append(thicknessSiN / SiNMFP * Water2MFP)

    # Fit to exponential.
    params, _ = curve_fit(mathHelper.expFunc, energies, equivalentThicknessWater, p0=[42.6, -0.133, -2.78],
                          method='dogbox')
    fitx = np.linspace(0, 300, num=300)
    fity = []
    for x in fitx:
        fity.append(mathHelper.expFunc(x, params[0], params[1], params[2]))

    plt.plot(energies, equivalentThicknessWater, 'ko')
    plt.plot(fitx, fity)
    plt.ylabel('Thickness of liquid water (oxygen) equivalent to 10nm Si3N4')
    plt.xlabel('Energy (keV)')
    plt.show()
