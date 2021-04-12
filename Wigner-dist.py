#This code creates the initial conditions based on the wigner distribution 
import numpy as np

numTrajs=500
Rc=2.6
m=918.0
omega=0.05

sigma=np.sqrt(1/(omega*m))

R0 = np.zeros(numTrajs)
for i in np.arange(numTrajs):
    R0[i] = np.random.normal(Rc,sigma)
   
P0 = np.zeros(numTrajs)
for i in np.arange(numTrajs):
    P0[i] = (np.random.normal(0,0.5*sigma**(-1)))/m


for i in np.arange(numTrajs):
	print(R0[i],P0[i])
