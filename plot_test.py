import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pdb
import h5py

numBladaff  = 1
numEUSaff   = 1
numPAGaff   = 1
numIND      = 1
numHypo     = 1
numINmplus  = 1
numINmminus = 1
numPGN      = 1
numFB       = 1
numIMG      = 1
numMPG      = 1
numEUSmn    = 1
numBladmn   = 1

config_file = "simulation_config.json"

v_file = "./output/cell_vars.h5"
i = h5py.File(v_file,'r')
mem_pot = i['report']['LUT']['data']


# Plot spike raster ---------------------------------------
Blad_gids = np.arange(0,numBladaff)
EUS_gids = Blad_gids + numBladaff
PAG_gids = EUS_gids + numEUSaff
IND_gids = PAG_gids + numPAGaff
Hypo_gids = IND_gids + numIND
INmplus_gids = Hypo_gids + numHypo
INmminus_gids = INmplus_gids + numINmplus
PGN_gids = INmminus_gids + numINmminus
FB_gids = PGN_gids + numPGN
IMG_gids = FB_gids + numFB
MPG_gids = IMG_gids + numIMG
EUSmn_gids = MPG_gids + numMPG
Bladmn_gids = EUSmn_gids + numEUSmn

df = pd.read_csv("output/spikes.csv",delimiter=' ')
rast = df.values


# Plot Bladder afferent, EUS afferent, PAG afferent, IND, and Hypo on one figure
plt.figure()
Bladspkt = rast[np.in1d(rast[:,-1],Blad_gids),:]
plt.plot(Bladspkt[:,0],Bladspkt[:,-1],'b.',label='Bladder afferent')

EUSspkt = rast[np.in1d(rast[:,-1],EUS_gids),:]
plt.plot(EUSspkt[:,0],EUSspkt[:,-1],'r.',label='EUS afferent')

PAGspkt = rast[np.in1d(rast[:,-1],PAG_gids),:]
plt.plot(PAGspkt[:,0],PAGspkt[:,-1],'m.',label='PAG afferent')

plt.xlabel('Time (t) [ms]')
plt.title('Afferent Firing Times')
plt.legend()

# Plot INM+, INM-, IND, FB --------------------------------
plt.figure()
INMpspkt = rast[np.in1d(rast[:,-1],INmplus_gids),:]
plt.plot(INMpspkt[:,0],INMpspkt[:,-1]-40,'b.',label='INm+')

INMmspkt = rast[np.in1d(rast[:,-1],INmminus_gids),:]
plt.plot(INMmspkt[:,0],INMmspkt[:,-1]-40,'r.',label='INM-')

INDspkt = rast[np.in1d(rast[:,-1],IND_gids),:]
plt.plot(INDspkt[:,0],INDspkt[:,-1]-30,'g.',label='IND')

FBspkt = rast[np.in1d(rast[:,-1],FB_gids),:]
plt.plot(FBspkt[:,0],FBspkt[:,-1]-50,'k.',label='FB')

plt.xlabel('Time (t) [ms]')
plt.title('Interneuron Firing Times')
plt.legend()

# Plot SPN(PGN), Hypo, MPG, IMG ---------------------------
plt.figure()
PGNspkt = rast[np.in1d(rast[:,-1],PGN_gids),:]
plt.plot(PGNspkt[:,0],PGNspkt[:,-1]-50,'b.',label='PGN')

IMGspkt = rast[np.in1d(rast[:,-1],IMG_gids),:]
plt.plot(IMGspkt[:,0],IMGspkt[:,-1]-80,'r.',label='IMG')

MPGspkt = rast[np.in1d(rast[:,-1],MPG_gids),:]
plt.plot(MPGspkt[:,0],MPGspkt[:,-1]-70,'c.',label='MPG')

Hypospkt = rast[np.in1d(rast[:,-1],Hypo_gids),:]
plt.plot(Hypospkt[:,0],Hypospkt[:,-1]-40,'m.',label='Hypo')

plt.xlabel('Time (t) [ms]')
plt.title('Ganglion/Preganglionic Firing Times')
plt.legend()

# Motor neurons -------------------------------------------
plt.figure()
EUSmnspkt = rast[np.in1d(rast[:,-1],EUSmn_gids),:]
plt.plot(EUSmnspkt[:,0],EUSmnspkt[:,-1],'r.',label='EUS MN')

Bladmnspkt = rast[np.in1d(rast[:,-1],Bladmn_gids),:]
plt.plot(Bladmnspkt[:,0],Bladmnspkt[:,-1],'g.',label='Bladder MN')

plt.xlabel('Time (t) [ms]')
plt.title('Motor Neuron Firing Times')
plt.legend()

plt.figure()
[H,B] = np.histogram(INMpspkt[:,0],bins=np.arange(0,10000,500))
plt.plot(B[:-1],H,'r',label='INM+')
[H,B] = np.histogram(INMmspkt[:,0],bins=np.arange(0,10000,500))
plt.plot(B[:-1],H,'b',label='INM-')
plt.legend()

plt.figure()
[H,B] = np.histogram(Bladspkt[:,0],bins=np.arange(0,10000,500))
plt.plot(B[:-1],H,'b',label='Bladder aff.')

[H,B] = np.histogram(PGNspkt[:,0],bins=np.arange(0,10000,500))
plt.plot(B[:-1],H,'r',label='PGN')

[H,B] = np.histogram(EUSspkt[:,0],bins=np.arange(0,10000,500))
plt.plot(B[:-1],H,'g',label='EUS')
plt.legend()
plt.show()
