# Phase1
Energy calculation was performed using "Dynamical_Matrix_AR_Kong.py," and the waves were created using one of the WaveCreator files. The other two filed, "Energy_Check.py" and "K_Chkeck.py," are for the post-processing of the data created by wave-packet simulations.

Dynamical_Matrix_Ar_Kong.py uses the input from Argon_Phonon.txt 

"Dispersion_Frequencies" and "energy_per_freqinterval" in the "WaveCreator" and "Energy_Check.py" are created at the end of Dynamical_Matrix_Ar_Kong.py

"dump.PE.allsections" in the "K_Check.py" is created by a LAMMPS run on the Validation system as mentioned in the paper

The Publication related to this repository is the following,

Phonon wave-packet simulations using the quantized definition of energy and a temperature-dependent phonon dispersion relation and phonon density of states
Amir Behbahanian and Nicholas A. Roberts
Phys. Rev. E 103, 043311 – Published 29 April 2021
