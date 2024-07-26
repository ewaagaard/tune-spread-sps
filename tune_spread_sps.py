"""
Main script to plot Q26 proton tune spread
"""
import matplotlib.pyplot as plt
import numpy as np
from PySCRDT.PySCRDT import PySCRDT
from matplotlib import colors
import resonance_lines
import os
import matplotlib
import xtrack as xt

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.size": 20,
        "axes.titlesize": 20,
        "axes.labelsize": 18,
        "xtick.labelsize": 18,
        "ytick.labelsize": 18,
        "legend.fontsize": 15,
        "figure.titlesize": 20,
    }
)


# Update beam parameters
Nb = 1.18e11 
exn = 0.8e-6
eyn = 0.8e-6
sigma_z = 0.22 #m
delta = 1.8e-3
bF = None #0.174 # None for Gaussian
Qh = 26.13  # to estimate detuning and plot
Qv = 26.25    # to estimate detuning and plot
plot_range = [[25.95, 26.35],[25.95, 26.35]]   # range in Qh & Qv for the plot
plot_order = 4   # order of resonances to plot
periodicity = 16  # periodicity of ring for the colorcode of the plot

# Parameters for tune scna
Qv_range = np.arange(26.15, 26.3, step=0.01)

# Load correct sequence
line = xt.Line.from_json('sequences/sps_q26_protons.json')
twiss = line.twiss()

print('Proton beam:')
print(line.particle_ref.show())

# ----------- Initiate PySCRDT object -----------------
s = PySCRDT()
    
s.setParameters(
    intensity=Nb,
    bunchLength=sigma_z,
    emittance_x=exn,
    emittance_y=eyn, 
    dpp_rms=delta,  
    bF=bF,
    ro=line.particle_ref.get_classical_particle_radius0()
)

s.loadTwissFromXsuite(twissTableXsuite=twiss)


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


for Qv in Qv_range:
    print('Qy = {:.2f}'.format(Qv))
    
    Qx, Qy= Qh, Qv 
    
    for x_q,y_q,detuning_coef in detuning:
        if x_q:
            Qx+=x_q/2.*detuning_coef*(J_x**(x_q/2.-1))*(J_y**(y_q/2.))
        if y_q:
            Qy+=y_q/2.*detuning_coef*(J_y**(y_q/2.-1))*(J_x**(x_q/2.))
    
    Q = np.dstack(
        (
            [qx.tolist() + [Qh] for qx in Qx],
            [qy.tolist() + [Qv] for qy in Qy],
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
    
    #detuning_y=[i for i in np.where(detuning[:,1]==2)[0] if i in np.where(detuning[:,0]==0)[0]][0]
    #detuning_x=[i for i in np.where(detuning[:,0]==2)[0] if i in np.where(detuning[:,1]==0)[0]][0]
    #print('dQx = ', str(detuning[detuning_x,2]))
    #print('dQy = ', str(detuning[detuning_y,2]))


    fig, ax = plt.subplots(1,figsize=(6, 6))
    tune_diagram = resonance_lines.resonance_lines(plot_range[0], plot_range[1], np.arange(1, plot_order+1), periodicity)
    tune_diagram.plot_resonance(figure_object=fig, interactive=False)
    ax.get_xaxis().get_major_formatter().set_useOffset(False)
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    ax.plot(Qh * np.ones(len(Qv_range)), Qv_range, ls='None', marker='o', ms=4.5, c='lime')
    ax.add_collection(p_collection)
    ax.set_aspect('equal')
    #ax.text(0.015, 0.64, '(Qx, Qy) =\n({:.2f}, {:.2f})'.format(Qh, Qv), fontsize=13, transform=ax.transAxes)
    ax.text(0.015, 0.64, 'Qy = {:.2f}'.format(Qv), fontsize=13, transform=ax.transAxes)
    fig.tight_layout()
    fig.savefig('plots/Q26_protons_Qx_dot{}_Qy_dot{}.png'.format(int((Qh % 1) * 100), int((Qv % 1) * 100)), dpi=250)
    del fig, ax
    plt.close()

