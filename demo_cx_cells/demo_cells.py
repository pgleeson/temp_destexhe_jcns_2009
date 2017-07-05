"""
 Network of IF cells - Sustained Activity


License: Modified BSD (see LICENSE.txt)
"""

import sys
from math import sqrt, pi
import numpy as np

SIMULATOR = sys.argv[-1]
exec("import pyNN.%s as pyNN" % SIMULATOR)
    
#-----------------------------------------------------------------
#  Parameters
#-----------------------------------------------------------------

# General parameters

SEED_LTS = 428577
SEED_CONN = 193566
SEED_GEN = 983651

DT = 0.1                                        # (ms) Time step
TSTART  = 0                                     
TSTOP   = 1200
V_INIT  = -60.0

# Cell parameters

LENGTH          = sqrt(20000/pi)                # in um
DIAMETER        = sqrt(20000/pi)                # in um
AREA            = 1e-8 * pi * LENGTH * DIAMETER # membrane area in cm2
TAU             = 20                            # time constant in ms
CAPACITANCE     = 1                             # capacitance in muF/cm2
G_L             = 1e-3 * CAPACITANCE / TAU      # leak conductance in S/cm2
V_REST          = -60                           # resting potential     

a_RS            = 0.001 
b_RS            = 0.1   # full adaptation
b_RS            = 0.005 # weaker adaptation


b_RSa            = 0.04 # Figure 1a
b_RSb            = 0.005 # Figure 1b


a_LTS           = 0.02
b_LTS           = 0.0
a_FS            = 0.001
b_FS            = 0.0

TAU_W           = 600
DELTA           = 2.5

# Spike parameters

VTR             = -50           # threshold in mV
VTOP            = 40            # top voltage during spike in mV
VBOT            = -60           # reset voltage in mV
REFRACTORY      = 5.0/2         # refractory period in ms (correction for a bug in IF_CG4)

# Synapse parameters

scale = 1               # scaling factor (>1 = more synapses)

TAU_E           = 5             # excitation time constant in ms
TAU_I           = 10            # inhibition time constant in ms
V_E             = 0             # excitatory reversal potential
V_I             = -80           # inhibitory reversal potential


N_GEN = 1           # total number of cells


MODEL_ID        = "cx05_cells_%s" % SIMULATOR

NEURONS_TO_RECORD = range(0,N_GEN)

# NEURON-specific parameters

NEURONS_TO_PLOT = NEURONS_TO_RECORD
USE_CVODE = False # }


#-----------------------------------------------------------------
#  Create cells
#-----------------------------------------------------------------

# we now use a standard cell model from PyNN, so there is nothing to do here

all_populations = []


def netCreate():
    global nLTS, neurons_RS, neurons_LTS, neurons_FS, all_populations
    RS_parameters = {
        'cm': 1000*AREA*CAPACITANCE, 'tau_m': TAU, 'v_rest': V_REST,
        'v_thresh': VTR, 'tau_refrac': REFRACTORY+DT,
        'v_reset': VBOT, 'v_spike': VTR+1e-6, 'a': 1000.0*a_RS, 'b': b_RS,
        'tau_w': TAU_W, 'delta_T': DELTA, 'tau_syn_E': TAU_E, 'e_rev_E': V_E,
        'tau_syn_I': TAU_I, 'e_rev_I': V_I
    }
    #RS_parameters['i_offset'] = 0.1
    
    RSa_parameters = RS_parameters.copy()
    RSa_parameters.update({'b': b_RSa}) 
    
    RSb_parameters = RS_parameters.copy()
    RSb_parameters.update({'b': b_RSb}) 
    
    
    LTS_parameters = RS_parameters.copy()
    LTS_parameters.update({'a': 1000.0*a_LTS, 'b': b_LTS}) # 1000 is for uS --> nS
    FS_parameters = RS_parameters.copy()
    FS_parameters.update({'a': 1000.0*a_FS, 'b': b_FS})

    neurons_RSa = pyNN.Population(N_GEN, pyNN.EIF_cond_exp_isfa_ista, RSa_parameters, label="RSa")
    neurons_RSb = pyNN.Population(N_GEN, pyNN.EIF_cond_exp_isfa_ista, RSb_parameters, label="RSb")
    
    neurons_LTS = pyNN.Population(N_GEN, pyNN.EIF_cond_exp_isfa_ista, LTS_parameters, label="LTS")
    neurons_FS = pyNN.Population(N_GEN, pyNN.EIF_cond_exp_isfa_ista, FS_parameters, label="FS")
    
    all_populations = [neurons_RSa,neurons_RSb,neurons_LTS,neurons_FS]


    ################################################################
    #IClamp
    step_current = pyNN.DCSource(amplitude= 0.25, start=200.0, stop=700.0)
    
    for pop in all_populations:
        step_current.inject_into(pop)
        pyNN.initialize(pop, v=VBOT)
        
    ################################################################


#-----------------------------------------------------------------
# Simulation settings
#-----------------------------------------------------------------

pyNN.setup(DT, min_delay=DT, use_cvode=USE_CVODE, rng_seeds_seed=SEED_GEN)

#-----------------------------------------------------------------
#  Add graphs
#-----------------------------------------------------------------

   

print ""
print "======================================================================="
print "            Network of ",N_GEN," * 3 IF neurons in an active state"
print "======================================================================="
print ""

#------------------------------------------------------------------------------
#  creating cells
#------------------------------------------------------------------------------
print "----[ CREATING CELLS ]----"
netCreate()




#-----------------------------------------------------------------
# Procedure to run simulation and menu
#-----------------------------------------------------------------

def run_sim():

        # record the Vm
        
        for pop in all_populations:
            pop.record('v')
        
        print "----[ RUNNING SIMULATION ]----"
        pyNN.run(TSTOP)
        
        pyNN.end()
        
        for pop in all_populations:
            data =  pop.get_data('v', gather=False)
            filename = "%s_v.dat"%(pop.label)
            print("Writing data for %s"%pop)
            for segment in data.segments:
                vm = segment.analogsignals[0].transpose()[0]
                tt = np.array([t*DT/1000. for t in range(len(vm))])
                times_vm = np.array([tt, vm/1000.]).transpose()
                np.savetxt(filename, times_vm , delimiter = '\t', fmt='%s')
        '''
        if SIMULATOR in ['neuron', 'nest', 'brian']:
            import matplotlib.pyplot as plt

            plt.figure(1)
            vmrs = neurons_RS.get_v()
            vmlts = neurons_LTS.get_v()
            vmfs = neurons_FS.get_v()

            times = []
            voltsrs = []
            voltslts = []
            voltsfs = []
            for i in range(N_GEN):
                times.append([])
                voltsrs.append([])
                voltslts.append([])
                voltsfs.append([])

            ids = vmrs[:,0]
            ts = vmrs[:,1]
            vsrs = vmrs[:,2]
            vslts = vmlts[:,2]
            vsfs = vmfs[:,2]
            for i in range(len(ids)):
                times[int(ids[i])].append(ts[i])
                voltsrs[int(ids[i])].append(vsrs[i])
                voltslts[int(ids[i])].append(vslts[i])
                voltsfs[int(ids[i])].append(vsfs[i])

            for i in range(N_GEN):
                plt.plot(times[i],voltsrs[i], '-')
                plt.plot(times[i],voltslts[i], '-')
                plt.plot(times[i],voltsfs[i], '-')

            plt.show()'''
            
        print "----[ FINISHED ]----"

run_sim()


