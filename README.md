# Tune spread generator - SPS Pb ions

A module to analytically estimate the space charge tune spread up to 3 beam sigmas using pre-estimated potentials from PySCRDT.

This forked version focuses in particular on estimating the SPS tune spread for Pb ions based on Twiss parameters from Xsuite Twiss. 
- Alternatively, a Twiss file can be included in the Twiss_files folder.

Update input_parameters.py with your desired parameters for estimating and plotting the tune spread.

Run tune_spread.py to get the tune footprint in the resonance plot, such as the one below.

![figure.png](./figure.png)
