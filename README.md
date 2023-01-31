
# Electron_transmission

This program performs estimations of electron transmission through different materials by using total elastic scattering cross section data to calculate mean free paths. More information onn how the calculations are performed can be found in the [wiki](https://github.com/UeIL-Waterloo/Electron_transmission/wiki).

This repository performs the calculations to generate Fig. 8 from:
> Dwyer, J. R.; Harb, M. Through a Window, Brightly: A Review of Selected Nanofabricated Thin-Film Platforms for Spectroscopy, Imaging, and Detection. Appl Spectrosc 2017, 71 (9), 2051–2075. https://doi.org/10.1177/0003702817715496.

The data for these calculations are taken from:
> Riley, M. E.; MacCallum, C. J.; Biggs, F. Theoretical Electron-Atom Elastic Scattering Cross Sections. Atomic Data and Nuclear Data Tables 1975, 15 (5), 443–476. https://doi.org/10.1016/0092-640X(75)90012-1.

Total elastic scattering cross section data for selected elements from the reference above are stored in `TotalCS_SelectedElements.csv`.
Moleuclar weight information for all elements are given in `PeriodicInfo.csv`.

The envrionment used to develop and run the program are as follows:
`Python 3.10.9`
> `matplotlib==3.6.3`
> `numpy==1.24.1`
> `pandas==1.5.3`
> `scipy==1.10.0`
