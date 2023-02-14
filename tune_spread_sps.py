import pylab as plt
import numpy as np
from PySCRDT.PySCRDT import PySCRDT
from matplotlib import colors
import resonance_lines
import input_parameters
import matplotlib

import acc_lib

import xobjects as xo
import xtrack as xt
import json

# ---------- Initiate Xtrack tracker to obtain Twiss parameters ----------- 
# Define file path for sequence 
fname_line = 'Sequence_and_Twiss_files/SPS_2021_Pb_ions_for_tracking.json'
fig_str = 'Pb_ions'
    
#%% Set up the X-suite contexts
context = xo.ContextCpu()
buf = context.new_buffer()   # using default initial capacity

#%% Load the sequence 
with open(fname_line, 'r') as fid:
     input_data = json.load(fid)
line = xt.Line.from_dict(input_data)
particle_ref = line.particle_ref

#%% Build tracker and perform Twiss on optics
tracker = xt.Tracker(_context=context,
                    line=line)
twiss_xtrack = tracker.twiss()  # optics adapted for sequence w/o SC

# ----------- Initiate PySCRDT object -----------------
s = PySCRDT()

# For SPS we use meters as input 
bunchLength = input_parameters.bunchLength_m
    
s.setParameters(
    intensity=input_parameters.intensity,
    bunchLength=bunchLength,
    emittance_x=input_parameters.emittance_x,
    emittance_y=input_parameters.emittance_y, 
    dpp_rms=input_parameters.dpp_rms,  
    bF=input_parameters.bF,
    ro=particle_ref.get_classical_particle_radius0()
)

s.loadTwissFromXsuite(twissTableXsuite=twiss_xtrack)


# caclulate detuning coefficients using potentials up to 20th order (needed for up to 3 sigma particles)
detuning=[]
# order in x
for i in range(0,21,2):
    print("--------- PySCRDT X order {} ------------- ".format(i))
    # order in y
    for j in range(0,21,2): 
        print("PySCRDT Y order {}".format(j))
        if (i==0) and (j==0):
            pass
        elif i+j<21:
            s.setOrder([int(i),int(j),'any'])
            s.potential()
            s.detuning()
            detuning.append([i,j,s.getDetuning()])
    
detuning=np.array(detuning)

#  initialize grid for calculation
s_N=6
s_max=3
theta_N=5

def initial_xy_polar(s_max, s_N, theta_N):
    return np.array(
        [
            [(s*np.cos(theta), s*np.sin(theta)) for s in np.linspace(0, s_max, s_N+1)]
            for theta in np.linspace(0, np.pi/2., theta_N)
        ])

S = initial_xy_polar(s_max=s_max, s_N=s_N, theta_N=theta_N)

#  estimate tunes from the detuning coefficients


en_x=s.parameters['emittance_x']
en_y=s.parameters['emittance_y']
beta=s.parameters['b']
gamma=s.parameters['g']
J_x=S[:,:,0]**2*en_x/2./beta/gamma
J_y=S[:,:,1]**2*en_y/2./beta/gamma

Qx,Qy=input_parameters.Qh,input_parameters.Qv 

for x_q,y_q,detuning_coef in detuning:
    if x_q:
        Qx+=x_q/2.*detuning_coef*(J_x**(x_q/2.-1))*(J_y**(y_q/2.))
    if y_q:
        Qy+=y_q/2.*detuning_coef*(J_y**(y_q/2.-1))*(J_x**(x_q/2.))

Q = np.dstack(
    (
        [qx.tolist() + [input_parameters.Qh] for qx in Qx],
        [qy.tolist() + [input_parameters.Qv] for qy in Qy],
    )
)

Q[:,:,0] += 0.00
Q[:,:,1] += 0.00

sx = Q.shape[0]-1
sy = Q.shape[1]-1
p1 = Q[:-1, :-1, :].reshape(sx*sy, 2)[:, :]
p2 = Q[1:, :-1, :].reshape(sx*sy, 2)[:]
p3 = Q[1:, 1:, :].reshape(sx*sy, 2)[:]
p4 = Q[:-1, 1:, :].reshape(sx*sy, 2)[:]


# do the plotting
cmap_base = plt.cm.hot
c_indcs = np.int_(np.linspace(0.1,0.6,s_N+1)*cmap_base.N)
cmap = colors.ListedColormap([cmap_base(c_indx) for c_indx in c_indcs])

# Stack endpoints to form polygons
Polygons = np.transpose(np.stack((p1, p2, p3, p4)), (1, 0, 2))
patches = list(map(matplotlib.patches.Polygon, Polygons))
p_collection = matplotlib.collections.PatchCollection(
#     patches, edgecolor='grey', linewidth=1,
    patches, edgecolor='k', linewidth=0.5,
#     facecolors=[],
    facecolors=cmap.colors,
#     facecolors=['SkyBlue'],
    alpha=0.7
)

detuning_y=[i for i in np.where(detuning[:,1]==2)[0] if i in np.where(detuning[:,0]==0)[0]][0]
detuning_x=[i for i in np.where(detuning[:,0]==2)[0] if i in np.where(detuning[:,1]==0)[0]][0]
print('dQx = ', str(detuning[detuning_x,2]))
print('dQy = ', str(detuning[detuning_y,2]))

fig,ax = plt.subplots(1,figsize=(10, 10))
tune_diagram = resonance_lines.resonance_lines(input_parameters.plot_range[0],
            input_parameters.plot_range[1], np.arange(1,input_parameters.plot_order+1), input_parameters.periodicity)
tune_diagram.plot_resonance(figure_object=fig, interactive=False)
ax.get_xaxis().get_major_formatter().set_useOffset(False)
ax.get_yaxis().get_major_formatter().set_useOffset(False)
ax.add_collection(p_collection)
ax.set_aspect('equal')
plt.tight_layout()
plt.savefig(input_parameters.figure, dpi=250)
plt.show()
