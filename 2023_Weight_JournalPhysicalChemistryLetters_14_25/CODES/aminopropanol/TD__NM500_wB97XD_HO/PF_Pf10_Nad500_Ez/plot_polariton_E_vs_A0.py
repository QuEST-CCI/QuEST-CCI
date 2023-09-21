import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import subprocess as sp

#A0_list = np.arange( 0.0, 0.5, 0.02 )
A0_list = np.round( np.arange( 0.0, 0.5 + 0.005, 0.005 ) ,4)
wc = 3.0
NM = 501
NF = 10
e = "z"
plot_states = np.arange(20)
DATA_DIR = "data_polariton"
PLOT_DIR = "data_plot"

#####################################

sp.call(f"mkdir -p {PLOT_DIR}",shell=True)

E_POL = np.zeros(( len(A0_list), NM*NF ))
F_POL = np.zeros(( len(A0_list), NM*NF ))
M_POL = np.zeros(( len(A0_list), NM*NF ))

for count, A0 in enumerate( A0_list ):
    print(f'Reading: A0 = {A0}')
    E_POL[count,:] = np.loadtxt( f"{DATA_DIR}/Epol_E{e}_eta{A0}_wc{wc}_Nad{NM-1}.dat" ) # Already in eV ... 630 kcal/mol/a.u. ; 27.2114 eV/a.u.
    F_POL[count,:] = np.loadtxt( f"{DATA_DIR}/Photon_Number_E{e}_eta{A0}_wc{wc}_Nad{NM-1}.dat" )
    M_POL[count,:] = np.loadtxt( f"{DATA_DIR}/Matter_Number_E{e}_eta{A0}_wc{wc}_Nad{NM-1}.dat" )

E_ZERO = np.min( E_POL[0,0] )

### SAVE DATA FOR FIRST 20 STATES ###
np.savetxt(f"{PLOT_DIR}/E_POL.dat", E_POL[:,:20] - E_ZERO)
np.savetxt(f"{PLOT_DIR}/Photon_Number.dat", F_POL[:,:20])
np.savetxt(f"{PLOT_DIR}/Matter_Number.dat", M_POL[:,:20])


#### PLOT ENERGY AS FUNCTION OF A0 COLORED BLACK ####

fig, ax = plt.subplots()

for count, state in enumerate(plot_states):
    ax.plot( A0_list, E_POL[:,state] - E_ZERO, c='black', linewidth=5, alpha=0.7, label=f"P$_{state}$" )

ax.set_xlim( A0_list[0], A0_list[-1] )
ax.set_ylim( A0_list[0] )
ax.set_xlabel("Coupling Strength, A$_0$",fontsize=15)
ax.set_ylabel("Polariton Energy (eV)",fontsize=15)
plt.tight_layout()
plt.savefig(f"{PLOT_DIR}/EPOL_vs_A0__Default.jpg" ,dpi=500 )
plt.clf()



#### PLOT ENERGY AS FUNCTION OF A0 COLORED BY AVERAGE PHOTON NUMBER ####

fig, ax = plt.subplots()
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)

for count, state in enumerate(plot_states):
    if ( count == 0 ):
        im = ax.scatter( A0_list, E_POL[:,state] - E_ZERO, c=F_POL[:,state], s=70,marker="o",edgecolor='none',cmap=plt.cm.viridis, vmin=0, vmax=4 )
    else:
        #ax.plot( A0_list, E_POL[:,state] - E_ZERO, label=f"P$_{state}$" )
        ax.scatter( A0_list, E_POL[:,state] - E_ZERO, c=F_POL[:,state], s=70,marker="o",edgecolor='none',cmap=plt.cm.viridis, vmin=0, vmax=4 )

cbar = fig.colorbar(im, cax=cax, orientation='vertical')
cbar.ax.set_ylabel('Average Photon Number', fontsize=15)#, rotation=270)
ax.set_xlim( A0_list[0], A0_list[-1] )
ax.set_ylim( A0_list[0] )
ax.set_xlabel("Coupling Strength, A$_0$",fontsize=15)
ax.set_ylabel("Polariton Energy (eV)",fontsize=15)
plt.tight_layout()
plt.savefig(f"{PLOT_DIR}/EPOL_vs_A0__Photon_Number.jpg" ,dpi=500 )
plt.clf()





#### PLOT ENERGY AS FUNCTION OF A0 COLORED BY AVERAGE MATTER EXCITATION ####

fig, ax = plt.subplots()
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)

for count, state in enumerate(plot_states):
    if ( count == 0 ):
        im = ax.scatter( A0_list, E_POL[:,state] - E_ZERO, c=M_POL[:,state], s=70,marker="o",edgecolor='none',cmap=plt.cm.hsv, vmin=0, vmax=30 )
    else:
        #ax.plot( A0_list, E_POL[:,state] - E_ZERO, label=f"P$_{state}$" )
        ax.scatter( A0_list, E_POL[:,state] - E_ZERO, c=M_POL[:,state], s=70,marker="o",edgecolor='none',cmap=plt.cm.hsv, vmin=0, vmax=30 )

cbar = fig.colorbar(im, cax=cax, orientation='vertical')
cbar.ax.set_ylabel('Average Photon Number', fontsize=15)#, rotation=270)
ax.set_xlim( A0_list[0], A0_list[-1] )
ax.set_ylim( A0_list[0] )
ax.set_xlabel("Coupling Strength, A$_0$",fontsize=15)
ax.set_ylabel("Polariton Energy (eV)",fontsize=15)
plt.savefig(f"{PLOT_DIR}/EPOL_vs_A0__Matter_Number.jpg" ,dpi=500 )
plt.clf()






