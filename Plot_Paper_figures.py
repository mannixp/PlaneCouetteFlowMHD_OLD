import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches, cm
import h5py, sys, os

from mpl_toolkits.axes_grid1.inset_locator import inset_axes
#import logging
#logger = logging.getLogger(__name__)

#plt.rcParams['pcolor.shading']
#plt.style.use('./methods_paper.mplstyle'); # Used to set font_size ticks size font type defaults etc

##########################################################################
# Fig 1. Time-Series in Semilogy
##########################################################################

def Plot_UB_scalar_data(file_name, marker_times):

	"""
	Plot the Velocity field & Magnetic fields integrated energies 

	Input Parameters:

	file_name with types .hdf5

	Returns: None

	"""

	file = h5py.File(file_name,"r")
	#print(file['scales/'].keys()); print(file['tasks/'].keys()) #useful commands	

	# Set the time interval we take
	index = 0; index_end = -1; 	
	time = file['scales/sim_time'][index:index_end];

	# All of these are <u,u> = (1/V)*(int_v u*u dV) where dV = rdr*ds*dz
	u2 = file['tasks/u_x total kinetic energy'][index:index_end,0,0,0];
	v2 = file['tasks/v kinetic energy'][index:index_end,0,0,0];
	w2 = file['tasks/w kinetic energy'][index:index_end,0,0,0];
	Kinetic = u2 + v2 + w2;
	print(len(Kinetic))
	print('int_t <U^2> dt = ',np.mean(Kinetic[2400:-1]))
	#ME = v2 + w2;

	A2 = file['tasks/B_x   magnetic energy'][index:index_end,0,0,0];
	B2 = file['tasks/B_y magnetic energy'][index:index_end,0,0,0];
	C2 = file['tasks/B_z magnetic energy'][index:index_end,0,0,0];
	Magnetic = A2 + B2 + C2;
	print(len(Magnetic))
	print('int_t <B^2> dt = ',np.mean(Magnetic[2400:-1]))
	#BME = B2 + C2;

	x = time;# - time[0]*np.ones(len(time)); # Modify time so that it's zero'd

	dpi = 1200;
	fig1, ax = plt.subplots(figsize=(12,6));

	ax.semilogy(x,Kinetic[0]*np.ones(len(x)) ,'b:',linewidth=1.5,label=r'$<U_0(z)^2>$');
	ax.semilogy(x,Kinetic ,'b-',linewidth=1.5,label=r'$<U^2>$');
	ax.semilogy(x,Magnetic,'k-',linewidth=1.5,label=r'$<B^2>$');

	for ii in range(len(marker_times)):
		for jj in range(len(x)):
			if x[jj] == marker_times[ii]:
				print("Times are equal")
				ax.semilogy(x[jj],Magnetic[jj],'ro',markersize=5.);

	ax.set_ylabel(r'$M$',color='k',fontsize=20)
	ax.set_xlabel(r'Time $t^* = t U^{-1}$',fontsize=20)

	'''
	inset_ax = inset_axes(a[0][1],width="60%",height="60%",loc=5)
	inset_ax.plot(x[0:100],BKE[0:100],'-',label=r'$k=%i$'%i);
	inset_ax.set_ylabel(r'$<B^2>$');
	inset_ax.set_xlim([0,np.max(x[0:100])])
	inset_ax.grid()
	`	'''

	#plt.grid()
	plt.legend(fontsize=20)
	plt.xlim([min(x),max(x)])
	plt.tight_layout(pad=1, w_pad=1.5)
	fig1.savefig('Kinetic_And_Magnetic_Energy.pdf', dpi=dpi);
	plt.show();

	return None;


##########################################################################
# Fig 2,3 Flux + max|Div.B|,max|Div.U|
##########################################################################

def Plot_MAX_DivUB(filename,dt):

	"""
	Plot the max( abs(div(B)) ) & max( abs(div(U)) ) during a Direct Adjoint Looping (DAL) routine in-terms of....
	divB = DAL_file['max|div(B)|'][()];
	divU = DAL_file['max|div(U)|'][()];

	Input Parameters:

	filename - hdf5 file with dict 'key':value structure, all 'keys' are defined as above

	Returns: None

	"""
	import matplotlib.pyplot as plt
	
	# Grab Data
	DAL_file = h5py.File(filename,"r")

	divB = DAL_file['max|div(B)|'][()];
	divU = DAL_file['max|div(U)|'][()];

	FLUX_B = DAL_file['FLUX_B'][()];
	
	DAL_file.close();

	# ~~~~~~~~~~~~~ # ~~~~~~~~~~~~~ # ~~~~~~~~~~~~~

	outfile = 'MaxDivU_MaxDivB_vs_Time.pdf'; 

	dpi = 1200
	fig, ax1 = plt.subplots(figsize=(8,6)); # 1,2,

	# Plot figures
	color = 'tab:red'
	x = np.arange(0,len(divB),1)*dt
	ax1.semilogy(x,divB,color=color, linestyle=':',linewidth=1.5, markersize=3);
	ax1.tick_params(axis='y', labelcolor=color)
	ax1.set_ylabel(r'Max|div(B)|',color=color,fontsize=18)
	ax1.set_xlabel(r'Time $t$',fontsize=18)

	ax1.set_xlim([0,np.max(x)])

	ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
	color = 'tab:blue'
	x = np.arange(0,len(divU),1)*dt
	ax2.semilogy(x,divU,color=color, linestyle='-',linewidth=1.5, markersize=3);
	
	ax2.tick_params(axis='y', labelcolor=color)
	ax2.set_ylabel(r'Max|div(U)|',color=color,fontsize=18)

	plt.grid()
	plt.tight_layout(pad=1, w_pad=1.5)
	fig.savefig(outfile, dpi=dpi);
	plt.show();

	# ~~~~~~~~~~~~~ # ~~~~~~~~~~~~~ # ~~~~~~~~~~~~~

	outfile1 = 'FluxB_vs_Time.pdf'; 
	fig1, ax = plt.subplots(figsize=(8,6));
	x = np.arange(0,len(FLUX_B),1)*dt

	ax.semilogy(x,FLUX_B,'k-',linewidth=1.5, markersize=3);
	ax.set_ylabel(r'|<B>|',color='k',fontsize=18)
	ax.set_xlabel(r'Time $t$',fontsize=18)

	plt.grid()
	plt.xlim([min(x),max(x)])
	plt.tight_layout(pad=1, w_pad=1.5)
	fig1.savefig(outfile1, dpi=dpi);
	plt.show();
	
	return None;


##########################################################################
# Fig 4 Slip-velocities
########################################################################## 

def Check_BCS(file_name):
	

	file = h5py.File(file_name,"r")
	print(file['scales/'].keys()); print(file['tasks/'].keys()) #useful commands	

	#sys.exit()

	#(time,x,y,z)
	time = file['scales/sim_time'][()];
	u    = file['tasks/u'];  v    = file['tasks/v'];  w   = file['tasks/w'];
	BZ_x = file['tasks/Az']; BZ_y = file['tasks/Bz']; B_z = file['tasks/C'];

	Times = len(u[:,0,0,0]);

	Slip_U_z0 = np.zeros(Times)
	Slip_U_z1 = np.zeros(Times)

	Slip_V_z0 = np.zeros(Times)
	Slip_V_z1 = np.zeros(Times)

	Slip_W_z0 = np.zeros(Times)
	Slip_W_z1 = np.zeros(Times)


	Slip_Az_z0 = np.zeros(Times)
	Slip_Az_z1 = np.zeros(Times)

	Slip_Bz_z0 = np.zeros(Times)
	Slip_Bz_z1 = np.zeros(Times)

	Slip_C_z0 = np.zeros(Times)
	Slip_C_z1 = np.zeros(Times);

	for index in range(Times):

		#print("\n Time t=%e"%time[index]);

		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # ~~~~~~~~~~~~~~~~~~~

		#print("\n")
		#print("Slip u-streamwise")
		maxU = np.amax(abs(u[index,:,:,:]));
		
		G1 = np.max(abs(u[index,:,:,0]));
		Slip_U_z0[index] = G1/maxU;
		
		G1 = np.max(abs(u[index,:,:,-1]));
		Slip_U_z1[index] = G1/maxU;
		

		#print("\n")
		#print("Slip v-spanwise")
		maxV = np.amax(abs(v[index,:,:,:]));

		G1 = np.max(abs(v[index,:,:,0]));
		Slip_V_z0[index] = G1/maxV;

		G1 = np.max(abs(v[index,:,:,-1]));
		Slip_V_z1[index] = G1/maxV;

		#print("\n")
		#print("Slip w-shearwise")
		maxW = np.amax(abs(w[index,:,:,:]));

		G1 = np.max(abs(w[index,:,:,0]));
		Slip_W_z0[index] = G1/maxW;
		
		G1 = np.max(abs(w[index,:,:,-1]));
		Slip_W_z1[index] = G1/maxW;

		#print("\n")

		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # ~~~~~~~~~~~~~~~~~~~

		#print("\n")
		#print("Slip dz(Bx)-streamwise")
		maxBZx = np.amax(abs(BZ_x[index,:,:,:]));
		
		G1 = np.max(abs(BZ_x[index,:,:,0]));
		Slip_Az_z0[index] = G1/maxBZx;

		G1 = np.max(abs(BZ_x[index,:,:,-1]));
		Slip_Az_z1[index] = G1/maxBZx;

		#print("\n")
		#print("Slip dz(By)-spanwise")
		maxBZy = np.amax(abs(BZ_y[index,:,:,:]));

		G1 = np.max(abs(BZ_y[index,:,:,0]));
		Slip_Bz_z0[index] = G1/maxBZy

		G1 = np.max(abs(BZ_y[index,:,:,-1]));
		Slip_Bz_z1[index] = G1/maxBZy

		#print("\n")
		#print("Slip Bz-shearwise")
		maxBz  = np.amax(abs(B_z[index,:,:,:]));

		G1 = np.max(abs(B_z[index,:,:,0]));
		Slip_C_z0[index] = G1/maxBz
		
		G1 = np.max(abs(B_z[index,:,:,-1]));
		Slip_C_z1[index] = G1/maxBz
		
		#print("\n")


	# ~~~~~~~~~~~~~~~~ # ~~~~~~~~~~~~~~~~~ 	# ~~~~~~~~~~~~~~~~ # ~~~~~~~~~~~~~~~~~
	print("\n")
	print("Velocity_Slip")
	dpi = 1200;
	fig1, ax = plt.subplots(figsize=(12,6));

	ax.semilogy(time,Slip_U_z0,'b:',linewidth=1.5,label=r'$u,z= 1$');
	ax.semilogy(time,Slip_U_z1,'b-',linewidth=1.5,label=r'$u,z=-1$');
	
	ax.semilogy(time,Slip_V_z0,'k:',linewidth=1.5,label=r'$v,z= 1$');
	ax.semilogy(time,Slip_V_z1,'k-',linewidth=1.5,label=r'$v,z=-1$');

	ax.semilogy(time,Slip_W_z0,'r:',linewidth=1.5,label=r'$w,z= 1$');
	ax.semilogy(time,Slip_W_z1,'r-',linewidth=1.5,label=r'$w,z=-1$');

	ax.set_ylabel(r'$u_i(z=\pm1)/max(abs(u_i))$',color='k',fontsize=20)
	ax.set_xlabel(r'Time $t^* = t U^{-1}$',fontsize=20)

	plt.grid()
	plt.legend(fontsize=20)
	plt.xlim([min(time),max(time)])
	plt.tight_layout(pad=1, w_pad=1.5)
	fig1.savefig('Velocity_Slip.pdf', dpi=dpi);
	plt.show();

	# ~~~~~~~~~~~~~~~~ # ~~~~~~~~~~~~~~~~~ # ~~~~~~~~~~~~~~~~ # ~~~~~~~~~~~~~~~~~
	print("\n")
	print("Magnetic_Slip")
	dpi = 1200;
	fig1, ax = plt.subplots(figsize=(12,6));

	ax.semilogy(time,Slip_Az_z0,'b:',linewidth=1.5,label=r'$dz(B_x),z= 1$');
	ax.semilogy(time,Slip_Az_z1,'b-',linewidth=1.5,label=r'$dz(B_x),z=-1$');
	
	ax.semilogy(time,Slip_Bz_z0,'k:',linewidth=1.5,label=r'$dz(B_y),z= 1$');
	ax.semilogy(time,Slip_Bz_z1,'k-',linewidth=1.5,label=r'$dz(B_y),z=-1$');

	ax.semilogy(time,Slip_C_z0,'r:',linewidth=1.5,label=r'$B_z,z= 1$');
	ax.semilogy(time,Slip_C_z1,'r-',linewidth=1.5,label=r'$B_z,z=-1$');

	ax.set_ylabel(r'$B_i(z=\pm1)/max(abs(B_i))$',color='k',fontsize=20)
	ax.set_xlabel(r'Time $t^* = t U^{-1}$',fontsize=20)

	plt.grid()
	plt.legend(fontsize=20)
	plt.xlim([min(time),max(time)])
	plt.tight_layout(pad=1, w_pad=1.5)
	fig1.savefig('Magnetic_Slip.pdf', dpi=dpi);
	plt.show();
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # ~~~~~~~~~~~~~~~~~~~

	return None;


##########################################################################
# Fig 5, Time-averaged Plot <u,u>(k,m) & <B,B>(k,m) reduced spectra
########################################################################## 

def Plot_KUKB_pair(file_names,st_index,Just_B):

	"""
	Plot Kinetc <u,u>(k,m) & Magentic <B,B>(k,m) as 3D surfaces

	The Integrals are indexed by the waveniumbers as follows
	# k - axial wavenumber, 
	# m - azim  wavenumber, 
	# Radial dependancy has been integrated out

	Input Parameters:

	file_names = ['one_file','second_file'] - array like, with file-types .hdf5
	times - integer: indicies of the CheckPoint(s) you wish to plot.... Depends on how the analysis_tasks are defined in FWD_Solve_TC_MHD
	LEN - integer: Length of the filenames array
	CAD - int: the cadence at which to plot of the details of files contained
	Just_B - bool: If True Plots only the magnetic field components

	Returns: None

	"""

	#from mpl_toolkits.mplot3d import Axes3D
	#from matplotlib import cm
	#from matplotlib.ticker import LinearLocator, FormatStrFormatter

	# Make data.
	file = h5py.File(file_names,"r")
	print(file['tasks/'].keys()); print(file['scales/'].keys()); 
	time = file['scales/sim_time'][()];

	# Get wavenumbers and create a Mesh-grid for plotting
	kx = file['scales/kx']; ky = file['scales/ky']
	print(len(kx[:]));
	print(len(ky[:]));

	print("kx = ",kx[:],"\n")

	Ny = int( 1 + (len(ky[:]) - 1)/2);
	inds_z = np.r_[0:Ny]; #Nz+1:len(kz[:])
	print(inds_z)
	print(ky[inds_z])
	ky = ky[inds_z];#[::-1];
	kx = kx[:];
	X, Y = np.meshgrid(kx, ky);
	print("X = ",X.shape)
	print("Y = ",Y.shape)	
	#sys.exit();

	##########################################################################
	# ~~~~~~~~~~~~~~~~~~~ plotting Magnetic B ~~~~~~~~~~~~~~~~~~~~~
	##########################################################################

	if Just_B == True:

		# BE(kx,ky) = FFT{ int_z B**2 dz }
		B = np.zeros(file['tasks/BE per k'][st_index,:,inds_z,0].shape)
		for ii in range(len(time) - st_index):
			B += abs( file['tasks/BE per k'][st_index + ii,:,inds_z,0] ); #(time,k_x,k_y,0)
		
		B = B/len(time);
		BE = np.log10(abs(B) + 1e-16);

		fig = plt.figure(figsize=(8,6)); 
		ax = fig.gca(projection='3d'); dpi = 400;
		plt.title("Energy int_T <B,B>(k,m) dt Field Iter=i%i Time=t%i"%(k,index) );

		# Plot the surface.
		surf = ax.plot_surface(X, Y, BE[:,:].T, cmap=cm.Greys,linewidth=0, antialiased=False);

		# Customize the z axis.
		ax.set_zlim(-15.,2.)
		#ax.zaxis.set_major_locator(LinearLocator(10))
		#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
		ax.set_xlabel(r'$k_x$ - Streamwise',fontsize=18);
		ax.set_ylabel(r'$k_y$ - Spanwise  ',fontsize=18);
		ax.set_zlabel(r'$log10(\hat{E}_B(k_x,k_y))$',fontsize=18);

		# Save figure
		plt.tight_layout(pad=1, w_pad=1.5)
		fig.savefig('Magnetic_Reduced_Spectra_time_avgd.pdf', dpi=dpi)
		plt.show()

	elif Just_B == False:	
		##########################################################################
		# ~~~~~~~~~~~~~~~~~~~ plotting Velocity B ~~~~~~~~~~~~~~~~~~~~~
		##########################################################################

		# BE(kx,ky) = FFT{ int_z B**2 dz }
		B = np.zeros(file['tasks/BE per k'][st_index,:,inds_z,0].shape)
		for ii in range(len(time) - st_index):
			B += abs( file['tasks/BE per k'][st_index + ii,:,inds_z,0] ); #(time,k_x,k_y,0)
		
		B = B/len(time);
		BE = np.log10(abs(B) + 1e-16);

		#fig = plt.figure(figsize=(8,6)); 
		fig = plt.figure(figsize=plt.figaspect(0.5))
		ax = fig.add_subplot(1, 2, 1, projection='3d'); dpi = 1200;
		#ax = fig.gca(projection='3d'); 
		plt.title(r'Energy $\hat{E}_B = int_T <B,B>(k_x,k_y) dt$',fontsize=12)# Field Iter=i%i Time=t%i'%(k,index) );

		#cmaps['Perceptually Uniform Sequential'] = ['viridis', 'plasma', 'inferno', 'magma', 'cividis'];
		# Plot the surface.
		surf = ax.plot_surface(X, Y, BE[:,:].T, cmap=cm.Greys,linewidth=0, antialiased=False);

		# Customize the z axis.
		ax.set_zlim(np.min(BE),np.max(BE))
		ax.set_xlim(0,np.max(kx))
		ax.set_ylim(0,np.max(ky));

		#ax.zaxis.set_major_locator(LinearLocator(10))
		#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
		ax.set_xlabel(r'$k_x$ - Streamwise',fontsize=12);
		ax.set_ylabel(r'$k_y$ - Spanwise  ',fontsize=12);
		#ax.set_zlabel(r'$log_{10}(\hat{E}_B(m_{\phi},k_{z}))$',fontsize=12);

		ax.view_init(30, 30)

		##########################################################################
		# ~~~~~~~~~~~~~~~~~~~ plotting Magnetic U ~~~~~~~~~~~~~~~~~~~~~
		##########################################################################

		# KE(k,m) = FFT{ int_r U**2 dr }
		U = np.zeros(file['tasks/KE per k'][st_index,:,inds_z,0].shape)
		for ii in range(len(time) - st_index):
			U += abs( file['tasks/KE per k'][st_index + ii,:,inds_z,0] ); #(time,k_x,k_y,0)
		
		U = U/len(time);
		KE = np.log10(abs(U) + 1e-16);

		ax = fig.add_subplot(1, 2, 2, projection='3d')
		plt.title(r'Energy $\hat{E}_U = int_T <U,U>(k_x,k_y) dt$',fontsize=12)# Field Iter=i%i Time=t%i'%(k,index) );

		# Plot the surface.
		surf = ax.plot_surface(X, Y, KE[:,:].T, cmap=cm.Greys,linewidth=0, antialiased=False);
		#ax.plot_wireframe(X, Y, KE[:,::-1].T, rstride=10, cstride=10)

		# Customize the z axis.
		ax.set_zlim(np.min(KE),np.max(KE))
		ax.set_xlim(0,np.max(kx))
		ax.set_ylim(0,np.max(ky));

		ax.set_xlabel(r'$k_x$ - Streamwise',fontsize=12);
		ax.set_ylabel(r'$k_y$ - Spanwise  ',fontsize=12);
		#ax.set_zlabel(r'$log_{10}(\hat{E}_U(m_{\phi},k_{z}))$',fontsize=12);

		ax.view_init(30,30)

		# Save figure
		plt.tight_layout(pad=1, w_pad=1.5)
		fig.savefig('Kinetic_AND_Magnetic_Reduced_Spectra_time_avgd.pdf', dpi=dpi)
		plt.show()

	return None;


##########################################################################
# Plot a vec=U,B pair
##########################################################################

# CHECK
def Plot_UB_pair(file_names,times,LEN,CAD,Just_B):

	"""
	Plot the Magnetic & Velocity (Optional) fields, as determined by:
	1) the Magnetic Iduction equation 
	2) Full MHD equations

	Input Parameters:

	file_names = ['one_file','second_file'] - array like, with file-types .hdf5
	times - integer: indicies of the CheckPoint(s) you wish to plot.... Depends on how the analysis_tasks are defined in FWD_Solve_TC_MHD
	LEN - integer: Length of the filenames array
	CAD - int: the cadence at which to plot of the details of files contained
	Just_B - bool: If True Plots only the magnetic field components
	
	Returns: None
	
	"""

	for k in range(0,LEN,CAD):

		file = h5py.File(file_names[k],"r")
		#print(file['scales/'].keys()); print(file['tasks/'].keys()) #useful commands

		x = file['scales/x/1.5']; y = file['scales/y/1.5']; z = file['scales/z/1.5']
		#x = file['scales/x/1.0']; y = file['scales/y/1.0']; z = file['scales/z/1.0']
		#(time,x,y,z)
		
		if Just_B == False:
			u   = file['tasks/u']; v   = file['tasks/v']; w   = file['tasks/w'];
		
		B_x = file['tasks/A']; B_y = file['tasks/B']; B_z = file['tasks/C'];
		'''
		if Just_B == False:
			u   = file['tasks/nu_u']; v   = file['tasks/nu_v']; w   = file['tasks/nu_w'];
		
		B_x = file['tasks/G_A']; B_y = file['tasks/G_B']; B_z = file['tasks/G_C'];
		'''
		SLICE = 12; # This needs some modification

		for i in range(len(times)):
			
			index = times[i]; # The instant at which we plot out the vector

			outfile_U = "".join(['U_PLOTS_Iter_i%i_Time_t%i.pdf'%(k,index) ]);	
			outfile_B = "".join(['B_PLOTS_Iter_i%i_Time_t%i.pdf'%(k,index) ]);	
			
			##########################################################################
			# ~~~~~~~~~~~~~~~~~~~ plotting Magnetic B ~~~~~~~~~~~~~~~~~~~~~
			##########################################################################

			fig = plt.figure(figsize=(8,6))
			plt.suptitle("B Field Iter=i%i Time=t%i"%(k,index) ); dpi = 400;

			print("B_x shape =",B_x[index,SLICE,:,:].shape);
			#print("y =",y.shape);
			#print("z =",z.shape);
			maxBx = np.amax(B_x[i,:,:,:]);
			minBx = np.amin(B_x[i,:,:,:]);
			
			maxBy = np.amax(B_y[i,:,:,:]);
			minBy = np.amin(B_y[i,:,:,:]);

			maxBz = np.amax(B_z[i,:,:,:]);
			minBz = np.amin(B_z[i,:,:,:]);
			
			print("max(Bx) =",maxBx)
			print("min(Bx) =",minBx)
			
			print("max(Bz) =",maxBz)
			print("min(Bz) =",minBz)
			
			print("min(By) =",minBy)
			print("max(By) =",maxBy,"\n")

			#------------------------------------ #------------------------------------
			ax1 = plt.subplot(221)
			Y,Z = np.meshgrid(y,z);
			cs = ax1.contourf(Z,Y,B_x[index,SLICE,:,:].T,cmap='PuOr',levels=30)
			
			skip=(slice(None,None,2),slice(None,None,2))
			ax1.quiver(Z[skip],Y[skip],B_z[index,SLICE,:,:][skip].T,B_y[index,SLICE,:,:][skip].T,width=0.005);
			fig.colorbar(cs,ax=ax1);

			ax1.set_title(r'$B_x$, vecs - $(B_z,B_y)$');
			ax1.set_xlabel(r'Shearwise - $z$')#, fontsize=18)
			ax1.set_ylabel(r'Spanwise  - $y$')#, fontsize=18)

			#------------------------------------ #------------------------------------
			ax2 = plt.subplot(222)
			X,Z = np.meshgrid(x,z);
			SLICE_y = 12
			cs = ax2.contourf(Z,X,B_y[index,:,SLICE_y,:].T,cmap='PuOr',levels=30)
			
			skip=(slice(None,None,4),slice(None,None,4))
			ax2.quiver(Z[skip],X[skip], B_z[index,:,SLICE_y,:][skip].T,B_x[index,:,SLICE_y,:][skip].T,width=0.005);
			fig.colorbar(cs,ax=ax2);

			ax2.set_title(r'$B_y$, vecs - ($B_z,B_x$)');
			ax2.set_xlabel(r'Shearwise  - $z$')
			ax2.set_ylabel(r'Streamwise - $x$')

			#------------------------------------ #------------------------------------
			ax3 = plt.subplot(212)
			Y,X = np.meshgrid(y,x);
			#print("X =",X.shape);
			#print("Y =",Y.shape);
			#print("Bz =",B_z[index,:,:,SLICE].shape)
			cs = ax3.contourf(X,Y,B_y[index,:,:,SLICE],cmap='PuOr',levels=30)
			
			skip=(slice(None,None,4),slice(None,None,4))
			ax3.quiver(X[skip],Y[skip], B_x[index,:,:,SLICE][skip],B_y[index,:,:,SLICE][skip],width=0.005);
			fig.colorbar(cs,ax=ax3)#,loc='right');

			ax3.set_title(r'$B_y$, vecs - ($B_x,B_y$)');
			ax3.set_xlabel(r'Streamwise - $x$');#, fontsize=18)
			ax3.set_ylabel(r'Spanwise   - $y$')

			#------------------------------------ #------------------------------------
			# Save figure
			plt.tight_layout(pad=1, w_pad=1.5)
			fig.savefig(outfile_B, dpi=dpi)
			#plt.show()


			##########################################################################
			# ~~~~~~~~~~~~~~~~~~~ plotting Velocity U ~~~~~~~~~~~~~~~~~~~~~
			##########################################################################
			if Just_B == False:

				fig = plt.figure(figsize=(8,6))
				plt.suptitle("U Field Iter=i%i Time=t%i"%(k,index) ); dpi = 400;

				#------------------------------------ #------------------------------------
				ax1 = plt.subplot(221)
				Y,Z = np.meshgrid(y,z);
				cs = ax1.contourf(Z,Y,u[index,SLICE,:,:].T,cmap='PuOr',levels=10)

				skip=(slice(None,None,2),slice(None,None,2))
				ax1.quiver(Z[skip],Y[skip],w[index,SLICE,:,:][skip].T,v[index,SLICE,:,:][skip].T,width=0.005);
				fig.colorbar(cs,ax=ax1);

				ax1.set_title(r'$u$, vecs - $(w,v)$');
				ax1.set_xlabel(r'Shearwise - $z$')#, fontsize=18)
				ax1.set_ylabel(r'Spanwise  - $y$')#, fontsize=18)

				#------------------------------------ #------------------------------------
				ax2 = plt.subplot(222)
				X,Z = np.meshgrid(x,z);
				cs = ax2.contourf(Z,X,w[index,:,SLICE,:].T,cmap='PuOr',levels=10)
				
				skip=(slice(None,None,4),slice(None,None,4))
				ax2.quiver(Z[skip],X[skip], w[index,:,SLICE,:][skip].T,u[index,:,SLICE,:][skip].T,width=0.005);
				fig.colorbar(cs,ax=ax2);

				ax2.set_title(r'$w$, vecs - ($w,u$)');
				ax2.set_xlabel(r'Shearwise  - $z$')
				ax2.set_ylabel(r'Streamwise - $x$')

				#------------------------------------ #------------------------------------
				ax3 = plt.subplot(212)
				Y,X = np.meshgrid(y,x);
				#print("X =",X.shape);
				#print("Y =",Y.shape);
				#print("Bz =",B_z[index,:,:,SLICE].shape)
				cs = ax3.contourf(X,Y,w[index,:,:,SLICE],cmap='PuOr',levels=10)
				
				skip=(slice(None,None,4),slice(None,None,4))
				ax3.quiver(X[skip],Y[skip], u[index,:,:,SLICE][skip],v[index,:,:,SLICE][skip],width=0.005);
				fig.colorbar(cs,ax=ax3)#,loc='right');

				ax3.set_title(r'$w$, vecs - ($u,v$)');
				ax3.set_xlabel(r'Streamwise - $x$');#, fontsize=18)
				ax3.set_ylabel(r'Spanwise   - $y$')

				#------------------------------------ #------------------------------------
				# Save figure
				plt.tight_layout(pad=1, w_pad=1.5)
				fig.savefig(outfile_U, dpi=dpi)
				#plt.show()

	return None;


#####################################

if __name__ == "__main__":

	
	'''
	# Incase files unmerged
	from dedalus.tools  import post
	from mpi4py import MPI
	comm = MPI.COMM_WORLD
	#Checkpoints_filenames = ['CheckPoints/CheckPoints_s1.h5','CheckPoints/CheckPoints_s2.h5','CheckPoints/CheckPoints_s3.h5','CheckPoints/CheckPoints_s4.h5']
	#post.merge_sets("CheckPoints_Merged", Checkpoints_filenames,cleanup=False, comm=MPI.COMM_WORLD)
	post.merge_process_files("CheckPoints", cleanup=True, comm=MPI.COMM_WORLD);
	post.merge_process_files("scalar_data", cleanup=True, comm=MPI.COMM_WORLD);
	comm.Barrier();
	#sys.exit()
	'''

	try:
		atmosphere_file = h5py.File('./Params.h5', 'r+')
		print(atmosphere_file.keys())

		Re = atmosphere_file['Re'][()];
		Pm = atmosphere_file['Pm'][()];
		mu = atmosphere_file['mu'][()];

		alpha = atmosphere_file['alpha'][()];
		beta = atmosphere_file['beta'][()];

		# Numerical Params
		Nx = atmosphere_file['Nx'][()];
		Ny = atmosphere_file['Ny'][()];
		Nz = atmosphere_file['Nz'][()];
		dt = atmosphere_file['dt'][()];

		print("Nx = ",Nx);
		print("Ny = ",Ny);
		print("Nz = ",Nz);
		print("dt = ",dt,"\n");

		M_0 = atmosphere_file['M_0'][()];
		E_0 = atmosphere_file['M_0'][()];
		T   = atmosphere_file['T'][()]

		atmosphere_file.close()
	except:

		pass;		
		
	##########################################################################
	# Plot scalar data figures 1 -> 4.
	##########################################################################	
	file = h5py.File('CheckPoints/CheckPoints_s1.h5', 'r+')
	time = file['scales/sim_time'][()];
	Marker_times = [time[1],time[10],time[20]]
	Plot_UB_scalar_data('scalar_data/scalar_data_s1.h5',Marker_times) # Fig 1

	sys.exit()
	#Plot_MAX_DivUB('FWD_Solve_IVP_DIV_UB.h5',dt) # Fig 2 + 3 # Forgot to save the data for these as a readable array

	Check_BCS('CheckPoints/CheckPoints_s1.h5') # Fig 4

	print("\n ----> Scalar Data Plots Complete <------- \n")	
	
	##########################################################################
	# Plot Spectra figure 5
	##########################################################################

	Plot_KUKB_pair('CheckPoints/CheckPoints_s1.h5',0,False);

	print("\n ----> Vector Field Plots Complete <------- \n")

	'''
	Checkpoints_filenames = ['CheckPoints/CheckPoints_s1.h5']
	#Checkpoints_filenames = ['CheckPoints_iter_9.h5']

	LEN = len(Checkpoints_filenames);
	Plot_Cadence = 1;
	times = [0,1,2,3,4,5,6,7,8,9]; #First and Last Checkpoints
	times += [10,11,12,13,14,15,16,17,18,19]
	times += [20,21,22,23,24,25,26,27,28,29]
	times += [40,41,42,43,44,45,46,47,48,49,-1];
	times += [30,31,32,33,34,35,36,37,38,39]; #First and Last Checkpoints
	#times = [0,-1]; #First and Last Checkpoints
	Just_B = False; # If True only plots B-field

	Plot_UB_pair(Checkpoints_filenames,times,LEN,Plot_Cadence,Just_B)
	'''

	# Only do this if you want to remove the file
	#import shutil
	#shutil.rmtree('snapshots');
