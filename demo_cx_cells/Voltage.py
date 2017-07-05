"""
No Errors
SingleCell model, fully working. 
Stores the output data to *.dat
"""


import pyNN.neuron as sim
import numpy.random, os
from pyNN.neuron import *
from numpy import *
import matplotlib.pyplot as plt
from pyNN.utility import Timer, plotting
from quantities import nA

setup()
run_time = 1000 #[ms], 1s

cell_params = {'tau_m'      : 20.0,             # [ms]		Membrane time constant
               'tau_syn_E'  : 5.0,		# [ms]		Rise time of excitatory synaptic conductance (alpha function)
               'tau_syn_I'  : 10.0,		# [ms]		Rise time of the inhibitory synaptic conductance (alpha function)
               'tau_refrac' : 2.5,		# [ms]		Duration of refractory period
               'v_rest'     : -60.0,		# [mV]		Resting membrane potential (Leak reversal potential)
               'v_reset'    : -60.0,#-70.0,	# [mV]		Reset value for membrane potential after a spike
               'v_thresh'   : -50.0,		# [mV]		Spike initiation threshold
               'delta_T'    : 2.5,		# [mV]		Slope factor
               'tau_w'      : 600.0,		# [ms]		Adaptation time constant
	     # 'v_spike'    : -40.0,            # added, [mV]	Spike detection threshold, default is -40
               'cm'         : 0.200,            # [nF]		Capacity of the membrane
               'a'          : 0.03e3,           # [nS], article: 0.001 uS = 0.001e3
               'b'          : 0.08 } 	        # = 0.0 for LTS, TC; =0.08 for RE  [nA]   
V_INIT = -60

def neo_spiketrains_to_numpy_gdf(segment):
    """
    takes a neo segment and converts the spiketrains contained therein into a
    numpy array in gdf format (1st column neuron ID, 2nd column time of spike).
    """
    spiketrains = segment.spiketrains
    allspikes_list = []
    for s in spiketrains:
        a = np.zeros((len(s), 2))
        a[:,0] = s.annotations["source_id"]
        a[:,1] = s.times
        allspikes_list.append(a)
    allspikes = np.concatenate(tuple(allspikes_list), axis=0)
    seq = np.argsort(allspikes[:,1])
    allspikes = allspikes[seq]
    return allspikes

def neo_analogsignals_to_numpy(segment):
    """
    takes a neo segment and converts the analog signals contained therein into
    numpy arrays.
    Returns dictionary with the following keys:
    voltages: NxM numpy array with N voltage samples for M neurons.
    times: Nx1 array with the timepoints of each sample.
    ids: Mx1 array with the IDs of the M recorded neurons.
    """
    analogsignals = segment.analogsignals
    if len(analogsignals) > 1:
        raise(Exception("more than 1 analog signal found, don't know what to do."))
    voltages = np.array(analogsignals[0])
    times = analogsignals[0].times
    ids = analogsignals[0].annotations['source_ids']
    return {"voltages": voltages,
            "times": times,
            "ids": ids}

cell = Population(1, EIF_cond_alpha_isfa_ista, cellparams=cell_params, initial_values={'v': V_INIT})

################################################################
#IClamp
step_current = DCSource(amplitude= 0.25, start=200.0, stop=700.0)
step_current.inject_into(cell)
################################################################

#Recording
cell.record('spikes')
cell.record('v')

#Running
sim.run(1000)

#Getting data
data_sp = cell.get_data().segments[0].spiketrains
data_vt = cell.get_data().segments[0].analogsignals
#data_vt.times = cell.get.data().segments[0].times
sim.end()

import numpy as np
block = cell.get_data()
segment = block.segments[0]
print "Segment::::::::::::::::::::", segment
analogsignals = segment.analogsignals
voltages = np.array(analogsignals[0])
times = analogsignals[0].times
ids = analogsignals[0].annotations['source_ids']



print "SpikeTimes:::::::::::::::::::::::::::::",data_sp
print "Voltage::::::::::::::::::::::::::::::::",data_vt
#print "Voltage times::::::::::::::::::::::::::",data_vt.times #AttributeError: 'list' object has no attribute 'times'

block = cell.get_data()
segment = block.segments[0]
allspikes = neo_spiketrains_to_numpy_gdf(segment)
np.savetxt('spikes.dat', allspikes, delimiter=",", fmt='%d, %f')
analogsigs = neo_analogsignals_to_numpy(segment)
np.savetxt('voltages.dat', analogsigs['voltages'], delimiter=',')
np.savetxt('times.dat', analogsigs['times'])
np.savetxt('ids.dat', analogsigs['ids'], fmt='%d')


end()
