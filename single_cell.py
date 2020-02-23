from bmtk.builder.networks import NetworkBuilder
from bmtk.utils.sim_setup import build_env_bionet
from bmtk.simulator import bionet
from bmtk.simulator.bionet.default_setters.cell_models import loadHOC

net = NetworkBuilder('LUT')
net.add_nodes(N=1, 
              level='high',
              pop_name='IMG',
              model_type='biophysical',
              model_template='hoc:IMG',
              morphology='blank.swc')
              
net.build()
net.save_nodes(output_dir='network')

build_env_bionet(base_dir='single_cell',      # Where to save the scripts and config files 
                 network_dir='network',    # Location of directory containing network files
                 tstop=2000.0, dt=0.1,     # Run a simulation for 2000 ms at 0.1 ms intervals
                 v_init=-54,
                 report_vars=['v','ina','ik','ninf_k1','minf_na1','hinf_na1'], # Tells simulator we want to record membrane potential and calcium traces
                 current_clamp={           # Creates a step current from 500.ms to 1500.0 ms  
                     'amp': 0.5,
                     'delay': 500.0,
                     'duration': 750.0
                 },
                 include_examples=True,    # Copies components files
                 compile_mechanisms=True   # Will try to compile NEURON mechanisms
                )
                
bionet.pyfunction_cache.add_cell_model(loadHOC, directive='hoc', model_type='biophysical')

conf = bionet.Config.from_json('single_cell/simulation_config.json')
conf.build_env()
net = bionet.BioNetwork.from_config(conf)
sim = bionet.BioSimulator.from_config(conf, network=net)
sim.run()



import h5py
import pdb
import numpy as np
import matplotlib.pyplot as plt

t = np.arange(0,2000,0.1)
v_file = './single_cell/output/v_report.h5'
ina_file = './single_cell/output/ina_report.h5'
ik_file = './single_cell/output/ik_report.h5'
ninf_file = './single_cell/output/ninf_k1_report.h5'
minf_file = './single_cell/output/minf_na1_report.h5'
hinf_file = './single_cell/output/hinf_na1_report.h5'

i = h5py.File(v_file,'r')
mem_pot = i['report']['LUT']['data']
i = h5py.File(ina_file,'r')
ina = i['report']['LUT']['data']
i = h5py.File(ik_file,'r')
ik = i['report']['LUT']['data']
i = h5py.File(ninf_file,'r')
ninf = i['report']['LUT']['data']
i = h5py.File(hinf_file,'r')
hinf = i['report']['LUT']['data']
i = h5py.File(minf_file,'r')
minf = i['report']['LUT']['data']

plt.figure(figsize=(10,10))
plt.subplot(3,1,1)
plt.plot(t,mem_pot)
v_diff = np.abs(mem_pot[4000]-np.min(mem_pot))
tau_v = mem_pot[4000]-0.632*v_diff
idx = (np.abs(mem_pot[:10000] - tau_v)).argmin()
tau = t[idx] - 500
plt.plot(t[idx],mem_pot[idx],'r.')
rin = (v_diff/1000)/0.2e-9
plt.ylim(-70,20)

plt.subplot(3,1,2)
plt.plot(t,ina)

plt.subplot(3,1,3)
plt.plot(t,ik)

# Only valid for current injection of -0.2 starting at 500 ms
plt.title('Rin: {} MOhm , tau: {} ms'.format(rin/1e6,tau))

plt.figure(figsize=(10,10))
plt.plot(mem_pot,ninf,'r.',label='ninf')
plt.plot(mem_pot,hinf,'m.',label='hinf')
plt.plot(mem_pot,minf,'b.',label='minf')
plt.legend()
plt.show()
