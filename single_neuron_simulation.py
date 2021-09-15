#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import time

import pylab

def remove_files():
    for f in os.listdir('.'):
        if '.gdf' in f or '.dat' in f:
            os.remove(f)

remove_files()


# In[2]:


import nest
import numpy as np
import random


# In[3]:


tic = time.perf_counter()

nest.Install("cerebmodule")
nest.set_verbosity("M_WARNING")

n_simulation = 0


# In[4]:


print(nest.GetDefaults('eglif_hippo_CA1'))


# In[4]:


while(n_simulation < 1):
    n_simulation += 1
    nest.ResetKernel()
    nest.SetKernelStatus({"overwrite_files": True,				# Parameters for writing on files
                        "data_path": "/home/nrp/workspace/cereb-nest",
                        "data_prefix": "eglif_CA1PC_"+str(n_simulation)+'_',
                        "resolution": 1.0})

    random.seed()
    seed = random.randint(10, 10000)
    print(seed)
    nest.SetKernelStatus({'grng_seed' : seed})


# In[6]:


#param_all = [3.4178, 0.3265, 283.0757, 1.2063, 62.5273, 70.3318] # median
#param_all = [0.4000, 0.0920, 141.5567, 1.1414, 262.4390, 100.0000]      #median optim2
#param_all = [0.623, 0.064184852, 315.456, 0.254, 4416.384, 0.0]              # Test fitting Annalisa
#param_CA1PC = {'t_ref': 2.15, 'C_m': 90.0, 'tau_m': 15.0, 'V_th': -48.0, 'V_reset': -55.0,'E_L': -65.0}
#param_CA1PC = {'t_ref': 2.15, 'C_m': 189.79, 'tau_m': 15.58, 'V_th': 26.5, 'V_reset': -58.5,'E_L': -65.0}     # Test fitting Annalisa


# In[31]:


eglif_hippo_CA1 = {
 #                   "t_ref": param_CA1PC['t_ref'],
  #                  "C_m": param_CA1PC['C_m'],
   #                 "V_th": param_CA1PC['V_th'],
    #                "V_reset": param_CA1PC['V_reset'],
     #               "E_L": param_CA1PC['E_L'],
      #              "Vinit": -65.0,
       #             "lambda_0": 1.0,
        #            "tau_V":0.0000000000001,
         #           "tau_m": param_CA1PC['tau_m'],
          #          "I_e": param_all[5],
           #         "kadap": param_all[0],
            #        "k1": param_all[3],

             #       "k2": param_all[1],
               #     "A1": param_all[4],
                #    "A2":param_all[2],
                 #   "tau_syn1": 1.0,
                  #  "tau_syn2": 0.7,
                   # "E_rev1": 0.0,
                    #"E_rev2": -80.0,
                    #"E_rev3": 0.0
    }


# In[32]:


single_neuron_CA1PC = nest.Create("eglif_hippo_CA1")


# In[33]:


# External input current
num_dc = 6
num_freq = 6 # Number of different frequency considered - literature protocol
num_step = 10 # Number of current step in each square wave period - literature protocol
num_tot = num_freq*num_step

current_dc_CA1PC = []


# In[34]:


for i in range(num_dc):
        current_dc_CA1PC.append(1)


for i in range(0,num_dc):
    current_dc_CA1PC[i] = nest.Create("dc_generator")


# In[35]:


sd = nest.Create('spike_detector',
                params = {"to_file": True,
                "withgid": True,
                "label": "spikes"})


# In[36]:


m = nest.Create("multimeter",
                    params = {"interval": 1.0,
                            "record_from": ["V_m", "V_th", "I_dep", "I_adap", "I_gen"],
                            "withgid": True,
                            "to_file": True,
                            "label": "multimeter"})


# In[37]:


current_amplitude = [0.0,1800.0,180.0,1600.0,180.0,0.0]
durate = [0.03, 0.1, 0.12, 0.05, 0.2, 1]


# In[38]:


nest.SetStatus(sd,[{"withgid": True, "withtime": True}])


# In[39]:


cont = 1
for i in durate[:num_dc]:
        nest.SetStatus(current_dc_CA1PC[cont-1],{'amplitude' :current_amplitude[cont-1],'start' : (np.sum(durate[:cont-1]))*1000.0, 'stop' : (np.sum(durate[:cont]))*1000.0})
        cont += 1


# In[40]:


nest.SetStatus(single_neuron_CA1PC, eglif_cond_alpha_multisyn)


# In[18]:


for i in range(0,num_dc):
        nest.Connect(current_dc_CA1PC[i], single_neuron_CA1PC)


# In[18]:


nest.Connect(m, single_neuron_CA1PC)
nest.Connect(single_neuron_CA1PC, sd)


# In[19]:


nest.Simulate(500.0)
time.sleep(0.5)


# In[ ]:


print(n_simulation)


# In[20]:


nest.GetDefaults('eglif_hippo_CA1')

toc = time.perf_counter()
print(f"time: {toc - tic:0.4f} seconds")

# In[ ]:


# Plotting
events = nest.GetStatus(m)[0]["events"]
t = events["times"]

pylab.clf()


pylab.subplot(411)
pylab.plot(t,events["I_gen"])

pylab.subplot(412)
pylab.plot(t,events["V_m"])

pylab.subplot(413)
pylab.plot(t,events["I_adap"])

pylab.subplot(414)
pylab.plot(t,events["I_dep"])

pylab.show()
