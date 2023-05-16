#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 17:30:02 2023

@author: livillenave
"""

plt
for i in range(6):
    plt.title(code_gauge[i])
    plt.plot(smash_model.output.qsim[i]); 
    plt.plot(smash_model.input_data.qobs[i])
    plt.show()
    plt.close()
    
    
    
import matplotlib.pyplot as plt
fig, (ax1, ax3) = plt.subplots(2, 3)
fig.set_figheight(7.5)
fig.set_figwidth(15)
fig
i=0
ax1[0].set_title(code_gauge[i])
ax1[0].plot(smash_model.output.qsim[i], label = "sim"); 
ax1[0].plot(smash_model.input_data.qobs[i], label = "obs")
ax1[0].set_ylabel('Q (m3/s)')
ax1[0].legend()
i=1
ax1[1].set_title(code_gauge[i])
ax1[1].plot(smash_model.output.qsim[i], label = "sim"); 
ax1[1].plot(smash_model.input_data.qobs[i], label = "obs")

i=2
ax1[2].set_title(code_gauge[i])
ax1[2].plot(smash_model.output.qsim[i], label = "sim"); 
ax1[2].plot(smash_model.input_data.qobs[i], label = "obs")

i=3

ax3[0].set_title(code_gauge[i])
ax3[0].plot(smash_model.output.qsim[i], label = "sim"); 
ax3[0].plot(smash_model.input_data.qobs[i], label = "obs")
ax3[0].set_xlabel('Index calculation timestep')
ax3[0].set_ylabel('Q (m3/s)')

i=4
ax3[1].set_title(code_gauge[i])
ax3[1].plot(smash_model.output.qsim[i], label = "sim"); 
ax3[1].plot(smash_model.input_data.qobs[i], label = "obs")
ax3[1].set_xlabel('Index calculation timestep')


i=5
ax3[2].set_title(code_gauge[i])
ax3[2].plot(smash_model.output.qsim[i], label = "sim"); 
ax3[2].plot(smash_model.input_data.qobs[i], label = "obs")
ax3[2].set_xlabel('Index calculation timestep')


plt.show()
plt.close()





import matplotlib.pyplot as plt
fig, (ax1) = plt.subplots(1, 4)
fig.set_figheight(7.5)
fig.set_figwidth(15)
fig
i=0
ax1[0].set_title(code_gauge[i])
ax1[0].plot(smash_model.output.qsim[i], label = "sim"); 
ax1[0].plot(smash_model.input_data.qobs[i], label = "obs")
ax1[0].set_ylabel('Q (m3/s)')
ax1[0].legend()
ax1[0].set_xlabel('Index calculation timestep')
i=1
ax1[1].set_title(code_gauge[i])
ax1[1].plot(smash_model.output.qsim[i], label = "sim"); 
ax1[1].plot(smash_model.input_data.qobs[i], label = "obs")
ax1[1].legend()
ax1[1].set_xlabel('Index calculation timestep')

i=2
ax1[2].set_title(code_gauge[i])
ax1[2].plot(smash_model.output.qsim[i], label = "sim"); 
ax1[2].plot(smash_model.input_data.qobs[i], label = "obs")
ax1[2].legend()
ax1[2].set_xlabel('Index calculation timestep')

i=3
ax1[3].set_title(code_gauge[i])
ax1[3].plot(smash_model.output.qsim[i], label = "sim"); 
ax1[3].plot(smash_model.input_data.qobs[i], label = "obs")
ax1[3].legend()
ax1[3].set_xlabel('Index calculation timestep')


plt.show()
plt.close()




smash_model.input_data.prcp[:,:,0]


sumrain = np.sum(smash_model.input_data.prcp, axis = 2)# * np.shape(smash_model.input_data.prcp)[2]
plt.imshow( sumrain  )
plt.colorbar(label = "Cumul de pluies [mm]" )




sumrain = np.sum(smash_model.input_data.pet, axis = 2)# * np.shape(smash_model.input_data.prcp)[2]
plt.imshow( sumrain  )
plt.colorbar(label = "Cumul des ETP imposées à Smash [mm]" )
