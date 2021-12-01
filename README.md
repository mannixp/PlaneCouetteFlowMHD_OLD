# PlaneCouetteFlowMHD
DNS script, initial condition and instructions required to simulate and reproduce the quasi-Keplarian plane Couette flow dynamo (Mannix et al. 2021/2). The script is written for the Dedalus code (Burns et al. 2020), which must be installed see [Dedalus](https://dedalus-project.org)

Having installed Dedalus and activated its conda environment, execute the script by running

`mpiexec -np procs python3 FWD_Solve_PCF_MHD.py`

where procs is to be replaced with the desired/availabe number. The time-series and snapshots produced by the simulation, can be viewed by executing

`python3 Plot_Paper_Data.py`

Figure 1. time-series of (Mannix et al. 2021/2) can be reproduced by specifying 

`Nx=64,Ny=128,Ny=64, dt=0.0125, N_ITERS = 10*((Pm*Re)/dt)`

along with the parameters specified therin, and running the script and plotting the data as above. 

Typically this will require a compute cluster or large workstation, for reference we find that using 32 cores on

`
Dual Intel Xeon Gold 6248 @ 2.5GHz – 20 cores processors - ??? cpu-hr
2 AMD Epyc 7302 @ 3 GHz – 16 cores processors            - ??? cpu-hr
`
Figure 1. snapshots of (Mannix et al. 2021/2) can be reproduced in Paraview. To prepare the data for Paraview: add the pyevtk to your conda environment

`conda install -c conda-forge pyevtk `

and execute

`python3 ViewHDF5_SLAB.py`


