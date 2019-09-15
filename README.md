# Mosfet-Poison-NEGF-solver

## Abstract : The 1D mosfet device N++-N+-N++ geometry are simulated by using the non-equilibrium Green function method. The numerical results are compared to two different numerical approaches, finite element and finite difference methods. The result shows that these two methods relative error are small, which means that the two methods are consistent in modeling the device physics.



## Numerical method references:

*    Two different approach using for simulating same problems. 
https://www.sciencedirect.com/science/article/pii/S0749603600909200


### The nonlinear poisson equation is solved by Newton–Raphson method.
https://en.wikipedia.org/wiki/Newton%27s_method



### Finite element method.
https://iopscience.iop.org/article/10.1088/0022-3727/42/10/105109/pdf

------------------------------------------------------------------------------------------------------------------------------

### Device Geometry :

![kk](https://github.com/Kuan-Ru-Chiou/Pic/blob/master/%E7%B0%A1%E5%A0%B11.jpg) 

x : electron transport direction

y and z are treated as infinite dimension approximation for dreducing the 3D problem to 1D.

--------------------------------------------------------------------------------------------------------------------------------
### Simulation Flowchart:

![kk](https://github.com/Kuan-Ru-Chiou/Pic/blob/master/4.jpg) 

-------------------------------------------------------------------------------------------------------------------------------
### Finite Difference Results :  

![kk](https://github.com/Kuan-Ru-Chiou/Pic/blob/master/1.png) 

Fig. 1 (a) electron density in a N++-N+-N++ structure. (b) Potential profile for electron in a N+-N-N+ structure. (c) Current voltage characteristics for the N++-N+-N++ structure. 


![kk](https://github.com/Kuan-Ru-Chiou/Pic/blob/master/2.png) 

Fig. 2 Energy position resolved energy/current profile for applied bias voltages V = 0.0625 V, 0.125 V, 0.1875 V, and 0.25 V in a N+-N-N+ structure. The quantum tunneling current can be seen in this simulation, which is different from thermal current.

### Finite Difference Compare To Finite Element Results:


![kk](https://github.com/Kuan-Ru-Chiou/Pic/blob/master/3.jpg) 

Fig. 3 Numerical results for finite difference method and finite element methods. The electron potential in the device and the source to drain current are shown in the left panels. Right panels are the relative errors for these two different methods.
