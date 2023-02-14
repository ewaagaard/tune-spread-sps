#input parameters for SPS Pb ions 
intensity       =       3.5e8
#bunchLength_ns  =       271 #4 sigma in ns
bunchLength_m	=	0.23 #1 sigma in m
emittance_x     =       1.3e-6
emittance_y     =       0.9e-6 
dpp_rms         =       1.0e-3
bF              =       None #0.174 # None for Gaussian
Qh		=	26.30  # to estimate detuning and plot
Qv		=	26.25    # to estimate detuning and plot
plot_range  =   [[25.9,26.4],[25.9,26.4]]   # range in Qh & Qv for the plot
plot_order  =   4   # order of resonances to plot
periodicity =   16  # periodicity of ring for the colorcode of the plot
figure		=	'figure.png'
twiss_file  =   'twiss_files/twiss_PSB'
