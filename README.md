# 1DGPE
Gross-Pitajevski equations solver in 1D through the Crank-Nicolson method with predictor-corrector. Realized in Fortran. It solves the system of two-coupled Gross-Pitajevskii equations describing the temporal dynamics of two interacting Bose-Einstein Condensates (BEC) in one dimension.

- Input Parameter file Aparameters.dat in the form (Example):
```
4.08775963342254       !g1
3.55457359428047       !g2
0.0                    !g12/g1	
3000                   !N1   
500                    !N2
-14.6027692307829      !z0 in units of lho=1.36960323647652e-06 m
```
- Interactions and Parameters:
	- N1 : number of particles in BEC of species 1
	- N2 : number of particles in BEC of species 2
	- g1 : intraspecies interactrion of BEC 1
	- g2 : intraspecies interactrion of BEC 2
	- g12: interspecies interaction between BEC 1 and BEC 2 (in units of g1)
  - z0 : initial displacement of the BEC 2 gaussian Wave-Packet 
