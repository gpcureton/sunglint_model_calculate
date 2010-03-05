;------------------------------------------------------------------------------
; The purpose of this program is to plot the modelled glint moments and 
; moment functions, using simulated slope cumulants and cumulant functions as 
; input. This is so that the modelled quantities can be compared to the
; simulated quantities calculated in cumulantFunctionSimulate.pro
;
; Inputs: fileName of simulation results
; Outputs: HDF4 file of modelled results
;
; Author: Geoff Cureton, 29th April, 2009
;------------------------------------------------------------------------------

PRO main,infileName,outFileName,windowIndex

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;; Read in input HDF4 file containing simulation results, and initialise   ;;;
	;;; various data structures.                                                ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	inFileName = STRING(inFileName)
	inFileName = STRCOMPRESS(inFileName,/REMOVE_ALL)
	outFileName = STRCOMPRESS(outFileName,/REMOVE_ALL)

	fileID = HDF_SD_START(inFileName, /READ)

	PRINT,"Input File:  ",inFileName
	PRINT,"Output File: ",outFileName

	PRINT, "fileID = ",fileID

	;END
	;;;;;;;;;

	PRINT, 'Reading some scale values...'

	sds_index = hdf_sd_nametoindex(fileID,"Solar Zenith Angles")
	sds_id = hdf_sd_select(fileID,sds_index)
	hdf_sd_getdata,sds_id,Solar_Zenith_Angles
	HDF_SD_ENDACCESS, sds_id

	sds_index = hdf_sd_nametoindex(fileID,"Detector Zenith Angles")
	sds_id = hdf_sd_select(fileID,sds_index)
	hdf_sd_getdata,sds_id,Detector_Zenith_Angles
	HDF_SD_ENDACCESS, sds_id

	sds_index = hdf_sd_nametoindex(fileID,"Specular Slopes")
	sds_id = hdf_sd_select(fileID,sds_index)
	hdf_sd_getdata,sds_id,Specular_Slopes
	HDF_SD_ENDACCESS, sds_id

	sds_index = hdf_sd_nametoindex(fileID,"Min Slopes")
	sds_id = hdf_sd_select(fileID,sds_index)
	hdf_sd_getdata,sds_id,Min_Slopes
	HDF_SD_ENDACCESS, sds_id

	sds_index = hdf_sd_nametoindex(fileID,"Max Slopes")
	sds_id = hdf_sd_select(fileID,sds_index)
	hdf_sd_getdata,sds_id,Max_Slopes
	HDF_SD_ENDACCESS, sds_id

	;;;;;;;;;

	PRINT, 'Reading some geometry values...'

	sds_index = hdf_sd_nametoindex(fileID,"Cumulant Function 1D length scale")
	sds_id = hdf_sd_select(fileID,sds_index)
	hdf_sd_getdata,sds_id,lag_1D
	HDF_SD_ENDACCESS, sds_id

	sds_index = hdf_sd_nametoindex(fileID,"Cumulant Function 2D length scale")
	sds_id = hdf_sd_select(fileID,sds_index)
	hdf_sd_getdata,sds_id,lag_2D
	HDF_SD_ENDACCESS, sds_id

	;;;;;;;;;

	PRINT, 'Reading the elevation, slope and glint moments and cumulants...'

	sds_index = hdf_sd_nametoindex(fileID,"Slope Moments")
	sds_id = hdf_sd_select(fileID,sds_index)
	hdf_sd_getdata,sds_id,slopeMoments
	dim_id = hdf_sd_dimgetid(sds_id,0)
	hdf_sd_dimget, dim_id, COUNT=N_moments 
	HDF_SD_ENDACCESS, sds_id

	sds_index = hdf_sd_nametoindex(fileID,"Glint Moments")
	sds_id = hdf_sd_select(fileID,sds_index)
	hdf_sd_getdata,sds_id,glintFirstMoments
	HDF_SD_ENDACCESS, sds_id

	sds_index = hdf_sd_nametoindex(fileID,"Slope Cumulants")
	sds_id = hdf_sd_select(fileID,sds_index)
	hdf_sd_getdata,sds_id,slopeCumulants
	HDF_SD_ENDACCESS, sds_id

	sds_index = hdf_sd_nametoindex(fileID,"Glint Cumulants")
	sds_id = hdf_sd_select(fileID,sds_index)
	hdf_sd_getdata,sds_id,glintCumulants
	HDF_SD_ENDACCESS, sds_id

	;;;;;;;;;

	PRINT, 'Reading the elevation, slope and glint second moment functions...'

	sds_index = hdf_sd_nametoindex(fileID,"Average Slope Second Moment Function")
	sds_id = hdf_sd_select(fileID,sds_index)
	hdf_sd_getdata,sds_id,slopeSecondMomentFunction
	dim_id = hdf_sd_dimgetid(sds_id,0)
	hdf_sd_dimget, dim_id, COUNT=N 
	HDF_SD_ENDACCESS, sds_id
																		  
	sds_index = hdf_sd_nametoindex(fileID,"Average Glint Second Moment Function")
	sds_id = hdf_sd_select(fileID,sds_index)
	hdf_sd_getdata,sds_id,glintSecondMomentFunction
	HDF_SD_ENDACCESS, sds_id

	;;;;;;;;;

	PRINT, 'Reading the elevation, slope and glint second cumulant functions...'

	sds_index = hdf_sd_nametoindex(fileID,"Average Slope Second Cumulant Function")
	sds_id = hdf_sd_select(fileID,sds_index)
	hdf_sd_getdata,sds_id,slopeSecondCumulantFunction
	HDF_SD_ENDACCESS, sds_id
																		  
	sds_index = hdf_sd_nametoindex(fileID,"Average Glint Second Cumulant Function")
	sds_id = hdf_sd_select(fileID,sds_index)
	hdf_sd_getdata,sds_id,glintSecondCumulantFunction
	HDF_SD_ENDACCESS, sds_id

	;;;;;;;;;

	PRINT, 'Reading the slope and glint third moment functions...'

	sds_index = hdf_sd_nametoindex(fileID,"Average Slope Third Moment Function")
	sds_id = hdf_sd_select(fileID,sds_index)
	hdf_sd_getdata,sds_id,slopeThirdMomentFunction
	dim_id = hdf_sd_dimgetid(sds_id,0)
	hdf_sd_dimget, dim_id, COUNT=NN 
	HDF_SD_ENDACCESS, sds_id
																		  
	sds_index = hdf_sd_nametoindex(fileID,"Average Glint Third Moment Function")
	sds_id = hdf_sd_select(fileID,sds_index)
	hdf_sd_getdata,sds_id,glintThirdMomentFunction
	dim_id = hdf_sd_dimgetid(sds_id,2)
	hdf_sd_dimget, dim_id, COUNT=N_angles
	HDF_SD_ENDACCESS, sds_id

	;;;;;;;;;

	PRINT, 'Reading the slope and glint third cumulant functions...'
																		  
	sds_index = hdf_sd_nametoindex(fileID,"Average Slope Third Cumulant Function")
	sds_id = hdf_sd_select(fileID,sds_index)
	hdf_sd_getdata,sds_id,slopeThirdCumulantFunction
	HDF_SD_ENDACCESS, sds_id
																		  
	sds_index = hdf_sd_nametoindex(fileID,"Average Glint Third Cumulant Function")
	sds_id = hdf_sd_select(fileID,sds_index)
	hdf_sd_getdata,sds_id,glintThirdCumulantFunction
	HDF_SD_ENDACCESS, sds_id

	;;;;;;;;;

	PRINT, "Closing HDF access to file ",inFileName
	HDF_SD_END, fileID

	d2r = !DPI/180.D
	r2d = 180.D/!DPI

	PRINT,""
	PRINT,"N = ",N
	PRINT,"NN = ",NN
	PRINT,"N_angles = ",N_angles
	PRINT,"N_moments = ",N_moments
	PRINT,'delta_x (1D) = ',lag_1D[1],' meters', Format='(A,F20.10,A)'
	PRINT,'delta_x (2D) = ',lag_2D[1],' meters', Format='(A,F20.10,A)'
	delta_x = lag_1D[1]
	PRINT,'Solar Zenith Angles (radians) = ',Solar_Zenith_Angles
	PRINT,'Solar Zenith Angles (degrees) = ',Solar_Zenith_Angles*r2d
	PRINT,'Detector Zenith Angles (radians) = ',Detector_Zenith_Angles
	PRINT,'Detector Zenith Angles (degrees) = ',Detector_Zenith_Angles*r2d
	PRINT,'Minimum Slopes = ',Min_Slopes
	PRINT,'Specular Slopes = ',Specular_Slopes
	PRINT,'Maximum Slopes = ',Max_Slopes

	PRINT,""
	;PRINT,"Simulated Elevation Moments(N_moments) =  ",elevMoments
	PRINT,"Simulated Slope Moments(N_moments) =      ",slopeMoments
	PRINT,"Simulated Glint First Moments(N_angles) = "
	PRINT,glintFirstMoments
	PRINT,""
	;PRINT,"Simulated Elevation Cumulants(N_angles) = ",elevCumulants
	PRINT,"Simulated Slope Cumulants(N_angles) =     ",slopeCumulants
	PRINT,"Simulated Glint Cumulants(N_angles,N_moments) = "
	PRINT,glintCumulants
	PRINT,"Simulated Glint Cumulants(N_moments,N_angles) = "
	PRINT,TRANSPOSE(glintCumulants)
	PRINT,""

	x_max = FLOAT(N-1L)*delta_x ; meters
	k_max = 2.D*!DPI/delta_x
	delta_k = k_max/DOUBLE(N-1L)
	k_N = DOUBLE(N/2L)*delta_k

	PRINT,'x_max =   ',x_max,' meters', Format='(A,F20.10,A)'
	PRINT,'k_max =   ',k_max,' meters^{-1}', Format='(A,F20.10,A)'
	PRINT,'delta_k = ',delta_k,' meters^{-1}', Format='(A,F20.10,A)'
	PRINT,'Nyquist Wavenumber = ',k_N,' meters^{-1}',Format='(A,F20.12,A)'

	;;; Turn on error reporting
	!EXCEPT=2

	xwinsize=800
	ywinsize=450

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;; Initialise the COMMON blocks, for passing data to other functions       ;;;
	;;; and subroutines which cannot be easily passed as parameters             ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	;;; Common block for the integration limits
	COMMON GEOMCOMMON, xiMin, xiMax

	;;; Common Block to pass slope cumulants to procedures and functions
	;;; which cannot passed as parameters to those procedures and functions.
	COMMON SLOPECUMULANTS,kappa_xi_1,kappa_xi_2,kappa_xi_3,kappa_xi_11,kappa_xi_21,kappa_xi_12, $
		kappa_xi_110,kappa_xi_101,kappa_xi_011, $
		kappa_xi_210,kappa_xi_120, $
		kappa_xi_201,kappa_xi_102, $
		kappa_xi_012,kappa_xi_021, $
		kappa_xi_111

	;;; Common Block containing the glint moments
	COMMON GLINTMOMENTS,mu_L_1,mu_L_2,mu_L_3,mu_L_11,mu_L_21,mu_L_12,mu_L_111

	;;; Common Block to pass lambda variables to procedures and functions
	;;; which cannot passed as parameters to those procedures and functions.
	COMMON LAMBDA,lambda2_xi_Matrix,lambda3_xi_Matrix,lambda_xi_Det,lambda2_xi_invMatrix,lambda3_xi_invMatrix

	;;; COMMON block for debugging
	;COMMON DEBUGGING, debugVar,windowIndex
	;windowIndex = -1

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;; Initialise the various structures, which are used to pass data          ;;;
	;;; in an encapsulated fashion.                                             ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	;;; Initialise the GEOM structure, which will be passed to the glint moment routine(s)

	GEOM = {N_angles:0L,source_angle:DINDGEN(N_angles),xi_min:DINDGEN(N_angles),$
		xi_0:DINDGEN(N_angles), xi_max:DINDGEN(N_angles)}

	beta = 0.68D*d2r

	GEOM.N_angles = N_angles
	GEOM.source_angle = Solar_Zenith_Angles
	gamma = (GEOM.source_angle - Detector_Zenith_Angles)/2.D
	GEOM.xi_0 = Specular_Slopes
	GEOM.xi_min = Min_Slopes
	GEOM.xi_max = Max_Slopes
	dxi = GEOM.xi_max - GEOM.xi_min

	xiMin = GEOM.xi_min
	xiMax = GEOM.xi_max

	;;; Initialise the SCALE structure, which will contain the various values of 
	;;; spatial and wavenumber increments

	SCALE = {N:123L,delta_x:0.D,x_max:0.D,k_max:0.D,delta_k:0.D,NN:64L}

	x = DBLARR(N/2L+1L)
	;x = DINDGEN(N/2L+1L)*delta_x
	x = lag_1D[0L:N/2L]

	SCALE.N = N
	SCALE.NN = NN
	SCALE.delta_x = delta_x
	SCALE.x_max = x_max
	SCALE.k_max = k_max
	SCALE.delta_k = delta_k

	;;; Initialise the SLOPE structure, which contains the slope moments and 
	;;; second and third slope moment functions
	SLOPE = {moments:DBLARR(N_moments),secondMomentFunction:DBLARR(N/2L+1L),thirdMomentFunction:DBLARR(NN/2L+1L,NN/2L+1L)}

	SLOPE.moments = slopeMoments
	SLOPE.secondMomentFunction = slopeSecondMomentFunction[0:N/2L]
	SLOPE.thirdMomentFunction = slopeThirdMomentFunction[0:NN/2L,0:NN/2L]

	;;; Initialise the GLINT structure, which contains the first glint moment, 
	;;; second glint moment function and third glint moment function for each
	;;; geometry defined in GEOM.

	thirdMomentFunc = {thirdMomentFunction:DBLARR(NN/2L+1L,NN/2L+1L)}
	thirdMomentFunctionGeom = REPLICATE(thirdMomentFunc, N_angles)

	GLINT = {geomIndex:0L,firstMoment:DBLARR(N_angles),secondMomentFunction:DBLARR(N_angles,N/2L+1L),$
		thirdMomentFunctionGeom:thirdMomentFunctionGeom}
	GLINT.firstMoment = glintFirstMoments
	GLINT.secondMomentFunction = TRANSPOSE(glintSecondMomentFunction[0:N/2L,*])
	;HELP,glintSecondMomentFunction
	;HELP,GLINT.secondMomentFunction
	;HELP,glintSecondCumulantFunction
	;HELP,glintThirdMomentFunction
	;HELP,glintThirdCumulantFunction

	FOR geometry=0L,N_angles-1L DO BEGIN
		GLINT.thirdMomentFunctionGeom[geometry].thirdMomentFunction = $
			glintThirdMomentFunction[0:NN/2L,0:NN/2L,geometry]
		;HELP,GLINT.thirdMomentFunctionGeom[geometry].thirdMomentFunction
	ENDFOR

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;; Initialise the simulated slope moments and moment functions, and the    ;;;
	;;; simulated slope cumulants and cumulant functions from quantities read   ;;;
	;;; from the input file.                                                    ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	;;; The simulated slope moments and cumulants
	mu_xi_sim = DBLARR(N_moments)
	mu_xi_sim = slopeMoments
	kappa_xi_sim = DBLARR(N_moments)
	kappa_xi_sim = slopeCumulants

	;;; The second order slope moment and cumulant functions
	M2_xi_sim = DBLARR(N/2L+1L)
	C2_xi_sim = DBLARR(N/2L+1L)
	M2_xi_sim = SLOPE.secondMomentFunction
	C2_xi_sim = slopeSecondCumulantFunction[0L:N/2L]

	;;; The main slope third moment and cumulant functions
	M3_xi_sim = DBLARR(NN/2L+1L,NN/2L+1L)
	C3_xi_sim = DBLARR(NN/2L+1L,NN/2L+1L)
	M3_xi_sim = SLOPE.thirdMomentFunction
	C3_xi_sim = slopeThirdCumulantFunction[0:NN/2L,0:NN/2L]

	;;; The special slope third moment and cumulant functions, along
	;;; the tau_1 and the (tau_1 = tau_2) axes
	M21_xi_sim = DBLARR(NN/2L+1L)
	C21_xi_sim = DBLARR(NN/2L+1L)
	M12_xi_sim = DBLARR(NN/2L+1L)
	C12_xi_sim = DBLARR(NN/2L+1L)

	FOR j=0L,NN/2L DO BEGIN
		;;; The third slope moment and cumulant functions along the tau_1 axis
		M21_xi_sim[j] = M3_xi_sim[0L,j]
		C21_xi_sim[j] = C3_xi_sim[0L,j]

		;;; The third slope moment and cumulant functions along the (tau_1 = tau_2) axis
		M12_xi_sim[j] = M3_xi_sim[j,j]
		C12_xi_sim[j] = C3_xi_sim[j,j]
	ENDFOR

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;; Compute the simulated glint cumulants and cumulant functions from the 	;;;
	;;; simulated glint moments and moment functions. 							;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	;;; The simulated first glint moments and cumulants
	mu_L1_sim = DBLARR(N_angles)
	mu_L1_sim = GLINT.firstMoment
	kappa_L_sim = DBLARR(N_angles,N_moments)
	kappa_L_sim = TRANSPOSE(glintCumulants)

	;;; The second order glint moment and cumulant functions
	M2_L_sim = DBLARR(N_angles,N/2L+1L)
	C2_L_sim = DBLARR(N_angles,N/2L+1L)
	M2_L_sim = GLINT.secondMomentFunction
	C2_L_sim = TRANSPOSE(glintSecondCumulantFunction[0L:N/2L,*])

	;;; The main glint third moment and cumulant functions
	M3_L_sim = DBLARR(N_angles,NN/2L+1L,NN/2L+1L)
	C3_L_sim = DBLARR(N_angles,NN/2L+1L,NN/2L+1L)
	FOR geometry=0L,N_angles-1L DO BEGIN
		M3_L_sim[geometry,*,*] = GLINT.thirdMomentFunctionGeom[geometry].thirdMomentFunction
	ENDFOR
	C3_L_sim = TRANSPOSE(glintThirdCumulantFunction[0L:NN/2L,0L:NN/2L,*])
	
	;;; The special glint third moment and cumulant functions, along
	;;; the tau_1 and the (tau_1 = tau_2) axes
	M21_L_sim = DBLARR(N_angles,NN/2L+1L)
	C21_L_sim = DBLARR(N_angles,NN/2L+1L)
	M12_L_sim = DBLARR(N_angles,NN/2L+1L)
	C12_L_sim = DBLARR(N_angles,NN/2L+1L)

	FOR geometry=0L,N_angles-1L DO BEGIN
		FOR j=0L,NN/2L DO BEGIN
			;;; The third glint moment function along the tau_1 axis
			M21_L_sim[geometry,j] = M3_L_sim[geometry,0L,j]
			C21_L_sim[geometry,j] = C3_L_sim[geometry,0L,j]

			;;; The third glint moment function along the (tau_1 = tau_2) axis
			M12_L_sim[geometry,j] = M3_L_sim[geometry,j,j]
			C12_L_sim[geometry,j] = C3_L_sim[geometry,j,j]
		ENDFOR
	ENDFOR
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;;		Model the glint moments and cumulants from   ;;;
	;;;     the simulated slope cumulants                ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	PRINT,'#####################################################',Format='(/A/)'
	PRINT,'Modelling the glint moments and cumulants...',Format='(A/)'

	mu_L1_modelled = DBLARR(N_angles)
	kappa_L_modelled = DBLARR(N_angles,N_moments)
	
	firstGlintMoment,GEOM,kappa_xi_sim,mu_L1_modelled,pderiv
	glintCumulantsFromMoments,mu_L1_modelled,kappa_L_modelled
	
	PRINT,"Modelled Glint First Moments = "
	PRINT,mu_L1_modelled
	PRINT,"Modelled Glint Cumulants ="
	PRINT,kappa_L_modelled

	;phillips_elev_spectrum,SCALE,POWER,NLCOUPLING
	;END

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;; Model the glint second moment functions for the various geometries,
	;;; using the simulated slope cumulants, and the modelled glint moments associated
	;;; with those cumulants, for comparison with the simulated glint second moment 
	;;; functions for the various geometries
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	PRINT,'#####################################################',Format='(/A/)'
	PRINT,'Modelling the glint second moment function...',Format='(A/)'

	;;; Define the lag variable using the basic spatial increment
	tau = DBLARR(N/2L+1L)
	tau = DINDGEN(N/2L+1L)*delta_x

	;;; Define the lag values at which the glint second moment function will be modelled
	;;; (not including the case of zero lag)
	N_tau_samples = 2L*30L   ;;; Must be an even number!!! We have this many samples 
	                         ;;; in each scale region.
	tauLo = tau[1]
	tauMid = 10.D
	tauHi = tau[N/2];50.D

	tau_index = LONARR(2L*N_tau_samples+1L)
	tau_reduced = DBLARR(2L*N_tau_samples+1L)
	tau_index[0L] = 0L
	tau_reduced[0L] = 0.D

	createTauAbcissa,N,N_tau_samples,delta_x,tauLo,tauMid,tauHi,tau,tau_index,tau_reduced

	;;; Initialise the SLOPECUMULANTS COMMON block with the simulated slope cumulants,
	;;; so they can be used in the modelling of the second glint moment and cumulant functions

	kappa_xi_1 = kappa_xi_sim[0]
	kappa_xi_2 = kappa_xi_sim[1]
	kappa_xi_3 = kappa_xi_sim[2]

	;;; Initialise the GLINTMOMENTS COMMON block with the modelled glint moments,
	;;; so they can be used in the modelling of the second glint moment and cumulant functions
	mu_L_1 = mu_L1_modelled

	M2_L_modelled = DBLARR(N_angles,N_tau_samples + 1L)
	M2_L_temp = DBLARR(N_angles)

	;;; Compute the modelled second glint functions using the simulated
	;;; slope cumulants and the simulated slope second cumulant function,
	;;; to check against the simulated linear second glint functions

	M2_L_modelled[*,0L] = 1.D

	FOR tau_sample=1L,N_tau_samples DO BEGIN

		PRINT,'tau_sample:',tau_sample,'   tau_index:',tau_index[tau_sample],$
			'   tau_reduced:',tau_reduced[tau_sample],'   delta_tau:',tau_reduced[tau_sample]-tau_reduced[tau_sample-1]

		;;; Redefine kappa_xi_guess to include slope cumulant functions 
		;;; kappa_xi_11, kappa_xi_21 and kappa_xi_12. The initial values
		;;; of kappa_xi_21 and kappa_xi_12 are set to zero for the second
		;;; order moment function

		kappa_xi_11 = kappa_xi_sim[1]*C2_xi_sim[tau_index[tau_sample]]

		kappa_xi_21 = 0.D
		kappa_xi_12 = 0.D
		;IF(tau_index[tau_sample] LT NN/2) THEN BEGIN
			;kappa_xi_21 = kappa_xi_sim[2]*C21_xi_sim[tau_index[tau_sample]]
			;kappa_xi_12 = kappa_xi_sim[2]*C12_xi_sim[tau_index[tau_sample]]
			;PRINT,"Using full coefficients",kappa_xi_21,kappa_xi_12
			;PRINT,"Using full coefficients",C21_xi_sim[tau_index[tau_sample]],C12_xi_sim[tau_index[tau_sample]]
		;ENDIF ELSE BEGIN
			;kappa_xi_21 = 0.D
			;kappa_xi_12 = 0.D
		;ENDELSE

		kappa_xi_guess = [kappa_xi_11,kappa_xi_21,kappa_xi_12]

		;;; Compute the second glint Moment function, with the slope cumulant 
		;;; functions kappa_xi_11, kappa_xi_21 and kappa_xi_12 as inputs. The 
		;;; simulated slope cumulants are passed through the COMMON block

		M2_L_temp[*] = 0.D
		secondGlintMomentFunction2,GEOM,kappa_xi_guess,M2_L_temp
		M2_L_modelled[*,tau_sample] = M2_L_temp
	ENDFOR

	;;; Interpolate the modelled second glint moment function to the delta_x spatial increment
	N_interp = N/2+1
	PRINT,"N_interp = ",N_interp
	M2_L_modelled_interp = DBLARR(N_angles,N_interp)
	FOR geometry=0L,N_angles-1L DO BEGIN
		HELP,M2_L_modelled_interp[geometry,*]
		M2_L_modelled_interp[geometry,*] = INTERPOL(M2_L_modelled[geometry,*],tau_reduced,tau[0L:N_interp-1L],/SPLINE)
	ENDFOR

	M2_L_modelled = DBLARR(N_angles,N)
	M2_L_modelled[*,0L:N/2L] = M2_L_modelled_interp[*,0L:N/2L]
	M2_L_modelled[*,N/2L+1L:N-1L] = REVERSE(M2_L_modelled_interp[*,1L:N/2L-1L],2)
	M2_L_modelled[*,0L] = 1L

	PRINT,"M2_L_modelled[*,0]       = ",M2_L_modelled[*,0]
	PRINT,"M2_L_modelled[*,1]       = ",M2_L_modelled[*,1]
	PRINT,"M2_L_modelled[*,N/2L-1L] = ",M2_L_modelled[*,N/2L-1L]
	PRINT,"M2_L_modelled[*,N/2L]    = ",M2_L_modelled[*,N/2L]
	PRINT,"M2_L_modelled[*,N/2L+1L] = ",M2_L_modelled[*,N/2L+1L]
	PRINT,"M2_L_modelled[*,N-1L]    = ",M2_L_modelled[*,N-1L]


	;;; Compute the modelled linear second glint cumulant function using the simulated slope cumulants
	;;; and the modelled glint cumulants

	C2_L_modelled = DBLARR(N_angles,N)
	FOR geometry=0L,N_angles-1L DO BEGIN
		C2_L_modelled[geometry,*] = $
			mu_L1_modelled[geometry]*(M2_L_modelled[geometry,*]-mu_L1_modelled[geometry])/$
				kappa_L_modelled[geometry,1]
	ENDFOR
	C2_L_modelled[*,0L] = 1L

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;;		Plot the simulated and modelled second glint functions	;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	
	SET_PLOT,"X"
	!P.FONT = -1
	!P.MULTI = 0
	chsize = 1.5

	abcissaTex="\tau (m)" 	;;; tau_{1} (m)
	ytitleTex="M^{(2)}_{L}(\tau)"		;;; (tau_{1},0)

	SET_PLOT,'X'
	!P.FONT = -1
	windowIndex +=1
	WINDOW,windowIndex,xsize = xwinsize,ysize = ywinsize,title='Simulated second glint moment function',RETAIN=2
	abcissaStr=TeXtoIDL(abcissaTex,FONT=-1)
	ytitleStr=TeXtoIDL(ytitleTex,FONT=-1)

	PLOT,tau,M2_L_sim[0,*],xtitle=abcissaStr,xrange=[0.,tauHi],$
		ytitle=ytitleStr,yrange=[-0.005,0.03],linestyle=0,/ystyle,charsize=chsize
	FOR angleIndex=1L,N_angles-1L DO BEGIN
		OPLOT,tau,M2_L_sim[angleIndex,*],linestyle=angleIndex
	ENDFOR

;	;;;;;;;;

	abcissaTex="\tau (m)" 	;;; tau_{1} (m)
	ytitleTex="C^{(2)}_{L}(\tau)"		;;; (tau_{1},0)

	SET_PLOT,'X'
	!P.FONT = -1
	windowIndex +=1
	WINDOW,windowIndex,xsize = xwinsize,ysize = ywinsize,title='Simulated second glint cumulant function',RETAIN=2
	abcissaStr=TeXtoIDL(abcissaTex,FONT=-1)
	ytitleStr=TeXtoIDL(ytitleTex,FONT=-1)

	PLOT,tau,C2_L_sim[0,*],xtitle=abcissaStr,$;xrange=[0.,tauHi],$
		ytitle=ytitleStr,yrange=[-0.002,0.002],linestyle=0,/ystyle,charsize=chsize
	FOR angleIndex=1L,N_angles-1L DO BEGIN
		OPLOT,tau,C2_L_sim[angleIndex,*],linestyle=angleIndex
	ENDFOR

 	;;;;;;;;
 
 	abcissaTex="\tau (m)" 	;;; tau_{1} (m)
 	ytitleTex="M^{(2)}_{L}(\tau)"		;;; (tau_{1},0)
 
	SET_PLOT,'X'
	!P.FONT = -1
	windowIndex +=1
	WINDOW,windowIndex,xsize = xwinsize,ysize = ywinsize,title='Modelled second glint moment function',RETAIN=2
	abcissaStr=TeXtoIDL(abcissaTex,FONT=-1)
	ytitleStr=TeXtoIDL(ytitleTex,FONT=-1)
 
    ;PLOT,tau,M2_L_modelled[0,*],xtitle=abcissaStr,$;xrange=[0.,tauHi],$
 	PLOT,M2_L_modelled[0,*],xtitle=abcissaStr,$;xrange=[0.,tauHi],$
 		ytitle=ytitleStr,yrange=[-0.005,0.03],linestyle=0,/ystyle,charsize=chsize
 	FOR angleIndex=1L,N_angles-1L DO BEGIN
         ;OPLOT,tau,M2_L_modelled[angleIndex,*],linestyle=angleIndex
 		OPLOT,M2_L_modelled[angleIndex,*],linestyle=angleIndex
 	ENDFOR
 
 	;;;;;;;;
 
 	abcissaTex="\tau (m)" 	;;; tau_{1} (m)
 	ytitleTex="C^{(2)}_{L}(\tau)"		;;; (tau_{1},0)
 
	SET_PLOT,'X'
	!P.FONT = -1
	windowIndex +=1
	WINDOW,windowIndex,xsize = xwinsize,ysize = ywinsize,title='Modelled second glint cumulant function',RETAIN=2
	abcissaStr=TeXtoIDL(abcissaTex,FONT=-1)
	ytitleStr=TeXtoIDL(ytitleTex,FONT=-1)
 
     ;PLOT,tau,C2_L_modelled[0,*],xtitle=abcissaStr,$;xrange=[0.,tauHi],$
 	PLOT,C2_L_modelled[0,*],xtitle=abcissaStr,$;xrange=[0.,tauHi],$
 		ytitle=ytitleStr,yrange=[-0.002,0.002],linestyle=0,/ystyle,charsize=chsize
 	FOR angleIndex=1L,N_angles-1L DO BEGIN
         ;OPLOT,tau,C2_L_modelled[angleIndex,*],linestyle=angleIndex
 		OPLOT,C2_L_modelled[angleIndex,*],linestyle=angleIndex
 	ENDFOR
 
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;; Model the glint third moment functions along a couple of slices, for the 
	;;; various geometries, using the simulated slope cumulants, and the modelled 
	;;; glint moments associated with those cumulants, for comparison with the 
	;;; simulated glint second moment 
	;;; functions for the various geometries
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	PRINT,'#####################################################',Format='(/A/)'
	PRINT,'Modelling the third glint moment function along the tau_1 axis (21)...',Format='(/A/)'

	;;; Define the complete lag variable using the basic spatial increment
	tau = DBLARR(N/2L+1L)
	tau = DINDGEN(N/2L+1L)*delta_x

	M3_L_maxIndex = NN/2L

	;;; Initialise the SLOPECUMULANTS COMMON block with the simulated slope cumulants,
	;;; so they can be used in the modelling of the third glint moment and cumulant functions

	kappa_xi_1 = kappa_xi_sim[0]
	kappa_xi_2 = kappa_xi_sim[1]
	kappa_xi_3 = kappa_xi_sim[2]

	;;; Initialise the GLINTMOMENTS COMMON block with the modelled glint moments
	mu_L_1 = mu_L1_modelled

	M21_L_modelled = DBLARR(N_angles,NN/2L+1L)
	M21_L_temp = DBLARR(N_angles)

	;;; Compute the modelled linear second glint functions using the simulated
	;;; slope cumulants and the simulated slope second cumulant function,
	;;; to check against the simulated linear second glint functions

	M21_L_modelled[*,0L] = 1.D

	PRINT,"tau range --> ",tau[0],tau[M3_L_maxIndex] ;     kappa_xi_11       kappa_xi_21         kappa_xi_12"
	FOR tau1=1L,M3_L_maxIndex DO BEGIN

		;;; Redefine kappa_xi_guess to include slope cumulant functions 
		;;; kappa_xi_11, kappa_xi_21 and kappa_xi_12. The initial values
		;;; of kappa_xi_21 and kappa_xi_12 are tentative

		kappa_xi_11 = kappa_xi_sim[1]*C2_xi_sim[tau1]
		kappa_xi_21 = kappa_xi_sim[2]*C21_xi_sim[tau1]
		kappa_xi_12 = kappa_xi_sim[2]*C12_xi_sim[tau1]
		kappa_xi_guess = [kappa_xi_11,kappa_xi_21,kappa_xi_12]


		M21_L_temp[*] = 0.D
		thirdGlintMomentFunction21,GEOM,kappa_xi_guess,M2_L_modelled[*,tau1],M21_L_temp
		M21_L_modelled[*,tau1] = M21_L_temp

		;PRINT,kappa_xi_guess,M21_L_modelled[*,tau1]
	ENDFOR

	;;; Compute the modelled linear second glint cumulant function using the retrieved slope cumulants.

	C21_L_modelled = DBLARR(N_angles,NN/2L+1L)
	FOR geometry=0L,N_angles-1L DO BEGIN
		;;; The third glint cumulant function along the tau_1 axis
		C21_L_modelled[geometry,*] = (mu_L1_modelled[geometry]*M21_L_modelled[geometry,*] $
			- mu_L1_modelled[geometry]*mu_L1_modelled[geometry]*(1.D + 2.D*M2_L_sim[geometry,0L:NN/2L]) $
			+ 2.D*mu_L1_modelled[geometry]^3.D)/kappa_L_modelled[geometry,2]
	ENDFOR

 	;;;;;;;;
 
 	postScriptOut=0L
 	abcissaTex="\tau_{1} (m)" 	;;; tau_{1} (m)
 	ytitleTex="M^{(3)}_{L}(\tau_{1},0)"		;;; (tau_{1},0)
 
	SET_PLOT,'X'
	!P.FONT = -1
	windowIndex +=1
	WINDOW,windowIndex,xsize = xwinsize,ysize = ywinsize,title='Simulated Third glint moment function - tau_x',RETAIN=2
	abcissaStr=TeXtoIDL(abcissaTex,FONT=-1)
	ytitleStr=TeXtoIDL(ytitleTex,FONT=-1)
 
 	PLOT,tau[0L:M3_L_maxIndex],M21_L_sim[0,0L:M3_L_maxIndex],xtitle=abcissaStr, $
 		xrange=[0.,M3_L_maxIndex*delta_x],$
 		ytitle=ytitleStr, $
 		yrange=[-0.01,0.04], $
 		linestyle=0,/ystyle,/xstyle,charsize=chsize
 	FOR angleIndex=1L,N_angles-1L DO BEGIN
 		OPLOT,tau[0L:M3_L_maxIndex],M21_L_sim[angleIndex,0L:M3_L_maxIndex],linestyle=angleIndex
 	ENDFOR
 
 	;;;;;;;;
 
 	postScriptOut=0L
 	abcissaTex="\tau_{1} (m)" 	;;; tau_{1} (m)
 	ytitleTex="C^{(3)}_{L}(\tau_{1},0)"		;;; (tau_{1},0)
 
	SET_PLOT,'X'
	!P.FONT = -1
	windowIndex +=1
	WINDOW,windowIndex,xsize = xwinsize,ysize = ywinsize,title='Simulated Third glint cumulant function - tau_x',RETAIN=2
	abcissaStr=TeXtoIDL(abcissaTex,FONT=-1)
	ytitleStr=TeXtoIDL(ytitleTex,FONT=-1)
 
 	PLOT,tau[0L:M3_L_maxIndex],C21_L_sim[0,0L:M3_L_maxIndex],xtitle=abcissaStr, $
 		xrange=[0.,M3_L_maxIndex*delta_x],$
 		ytitle=ytitleStr, $
 		yrange=[-0.01,0.01], $
 		linestyle=0,/ystyle,/xstyle,charsize=chsize
 	FOR angleIndex=1L,N_angles-1L DO BEGIN
 		OPLOT,tau[0L:M3_L_maxIndex],C21_L_sim[angleIndex,0L:M3_L_maxIndex],linestyle=angleIndex
 	ENDFOR
 
 	;;;;;;;;
 
 	postScriptOut=0L
 	abcissaTex="\tau_{1} (m)" 	;;; tau_{1} (m)
 	ytitleTex="M^{(3)}_{L}(\tau_{1},0)"		;;; (tau_{1},0)
 
	SET_PLOT,'X'
	!P.FONT = -1
	windowIndex +=1
	WINDOW,windowIndex,xsize = xwinsize,ysize = ywinsize,title='Modelled Third glint moment function - tau_x',RETAIN=2
	abcissaStr=TeXtoIDL(abcissaTex,FONT=-1)
	ytitleStr=TeXtoIDL(ytitleTex,FONT=-1)
 
 	PLOT,tau[0L:M3_L_maxIndex],M21_L_modelled[0,0L:M3_L_maxIndex],xtitle=abcissaStr, $
 		xrange=[0.,M3_L_maxIndex*delta_x],$
 		ytitle=ytitleStr, $
 		yrange=[-0.01,0.04], $
 		linestyle=0,/ystyle,/xstyle,charsize=chsize
 	FOR angleIndex=1L,N_angles-1L DO BEGIN
 		OPLOT,tau[0L:M3_L_maxIndex],M21_L_modelled[angleIndex,0L:M3_L_maxIndex],linestyle=angleIndex
 	ENDFOR
 
 	;;;;;;;;
 
 	postScriptOut=0L
 	abcissaTex="\tau_{1} (m)" 	;;; tau_{1} (m)
 	ytitleTex="C^{(3)}_{L}(\tau_{1},0)"		;;; (tau_{1},0)
 
	SET_PLOT,'X'
	!P.FONT = -1
	windowIndex +=1
	WINDOW,windowIndex,xsize = xwinsize,ysize = ywinsize,title='Modelled Third glint cumulant function - tau_x',RETAIN=2
	abcissaStr=TeXtoIDL(abcissaTex,FONT=-1)
	ytitleStr=TeXtoIDL(ytitleTex,FONT=-1)
 
 	PLOT,tau[0L:M3_L_maxIndex],C21_L_modelled[0,0L:M3_L_maxIndex],xtitle=abcissaStr, $
 		xrange=[0.,M3_L_maxIndex*delta_x],$
 		ytitle=ytitleStr, $
 		yrange=[-0.01,0.01], $
 		linestyle=0,/ystyle,/xstyle,charsize=chsize
 	FOR angleIndex=1L,N_angles-1L DO BEGIN
 		OPLOT,tau[0L:M3_L_maxIndex],C21_L_modelled[angleIndex,0L:M3_L_maxIndex],linestyle=angleIndex
 	ENDFOR
 
	PRINT,'#####################################################',Format='(/A/)'
	PRINT,'Modelling the third glint moment function along the tau_1=tau_2 axis (12)...',Format='(/A/)'

	;;; Initialise the SLOPECUMULANTS COMMON block with the simulated slope cumulants,
	;;; so they can be used in the modelling of the second glint moment and cumulant functions

	kappa_xi_1 = kappa_xi_sim[0]
	kappa_xi_2 = kappa_xi_sim[1]
	kappa_xi_3 = kappa_xi_sim[2]

	;;; Initialise the GLINTMOMENTS COMMON block with the modelled glint moments
	mu_L_1 = mu_L1_modelled

	M12_L_modelled = DBLARR(N_angles,NN/2L+1L)
	M12_L_temp = DBLARR(N_angles)

	;;; Compute the modelled linear second glint functions using the retrieved
	;;; slope cumulants and the simulated slope second cumulant function,
	;;; to check against the simulated linear second glint functions

	M12_L_modelled[*,0L] = 1.D

	PRINT,"tau range --> ",tau[0],tau[M3_L_maxIndex] ;     kappa_xi_11       kappa_xi_21         kappa_xi_12"
	FOR tau1=1L,M3_L_maxIndex DO BEGIN

		;;; Redefine kappa_xi_guess to include slope cumulant functions 
		;;; kappa_xi_11, kappa_xi_21 and kappa_xi_12. The initial values
		;;; of kappa_xi_21 and kappa_xi_12 are tentative

		kappa_xi_11 = kappa_xi_sim[1]*C2_xi_sim[tau1]
		kappa_xi_21 = kappa_xi_sim[2]*C21_xi_sim[tau1]
		kappa_xi_12 = kappa_xi_sim[2]*C12_xi_sim[tau1]
		kappa_xi_guess = [kappa_xi_11,kappa_xi_21,kappa_xi_12]

		M12_L_temp[*] = 0.D
		thirdGlintMomentFunction12,GEOM,kappa_xi_guess,M2_L_modelled[*,tau1],M12_L_temp
		M12_L_modelled[*,tau1] = M12_L_temp
	ENDFOR

	;;; Compute the modelled linear second glint cumulant function using the retrieved slope cumulants.

	C12_L_modelled = DBLARR(N_angles,NN/2L+1L)
	FOR geometry=0L,N_angles-1L DO BEGIN
		;;; The third glint cumulant function along the tau_1 axis
		C12_L_modelled[geometry,*] = (mu_L1_modelled[geometry]*M12_L_modelled[geometry,*] $
			- mu_L1_modelled[geometry]*mu_L1_modelled[geometry]*(1.D + 2.D*M2_L_sim[geometry,0L:NN/2L]) $
			+ 2.D*mu_L1_modelled[geometry]^3.D)/kappa_L_modelled[geometry,2]
	ENDFOR

 	;;;;;;;;
 
 	postScriptOut=0L
 	abcissaTex="\tau (m)" 	;;; tau_{1} (m)
 	ytitleTex="M^{(3)}_{L}(\tau,\tau)"		;;; (tau_{1},0)
 
	SET_PLOT,'X'
	!P.FONT = -1
	windowIndex +=1
	WINDOW,windowIndex,xsize = xwinsize,ysize = ywinsize,title='Simulated Third glint moment function - tau_x = tau_y',RETAIN=2
	abcissaStr=TeXtoIDL(abcissaTex,FONT=-1)
	ytitleStr=TeXtoIDL(ytitleTex,FONT=-1)
 
 	PLOT,tau[0L:M3_L_maxIndex],M12_L_sim[0,0L:M3_L_maxIndex],xtitle=abcissaStr, $
 		xrange=[0.,M3_L_maxIndex*delta_x],$
 		ytitle=ytitleStr, $
 		yrange=[-0.01,0.04], $
 		linestyle=0,/ystyle,/xstyle,charsize=chsize
 	FOR angleIndex=1L,N_angles-1L DO BEGIN
 		OPLOT,tau[0L:M3_L_maxIndex],M12_L_sim[angleIndex,0L:M3_L_maxIndex],linestyle=angleIndex
 	ENDFOR
 
 	;;;;;;;;
 
 	postScriptOut=0L
 	abcissaTex="\tau (m)" 	;;; tau_{1} (m)
 	ytitleTex="C^{(3)}_{L}(\tau,\tau)"		;;; (tau_{1},0)
 
	SET_PLOT,'X'
	!P.FONT = -1
	windowIndex +=1
	WINDOW,windowIndex,xsize = xwinsize,ysize = ywinsize,title='Simulated Third glint cumulant function - tau_x = tau_y',RETAIN=2
	abcissaStr=TeXtoIDL(abcissaTex,FONT=-1)
	ytitleStr=TeXtoIDL(ytitleTex,FONT=-1)
 
 	PLOT,tau[0L:M3_L_maxIndex],C12_L_sim[0,0L:M3_L_maxIndex],xtitle=abcissaStr, $
 		xrange=[0.,M3_L_maxIndex*delta_x],$
 		ytitle=ytitleStr, $
 		yrange=[-0.01,0.01], $
 		linestyle=0,/ystyle,/xstyle,charsize=chsize
 	FOR angleIndex=1L,N_angles-1L DO BEGIN
 		OPLOT,tau[0L:M3_L_maxIndex],C12_L_sim[angleIndex,0L:M3_L_maxIndex],linestyle=angleIndex
 	ENDFOR
 
 	;;;;;;;;
 
 	postScriptOut=0L
 	abcissaTex="\tau (m)" 	;;; tau_{1} (m)
 	ytitleTex="M^{(3)}_{L}(\tau,\tau)"		;;; (tau_{1},0)
 
	SET_PLOT,'X'
	!P.FONT = -1
	windowIndex +=1
	WINDOW,windowIndex,xsize = xwinsize,ysize = ywinsize,title='Modelled Third glint moment function - tau_x = tau_y',RETAIN=2
	abcissaStr=TeXtoIDL(abcissaTex,FONT=-1)
	ytitleStr=TeXtoIDL(ytitleTex,FONT=-1)
  
 	PLOT,tau[0L:M3_L_maxIndex],M12_L_modelled[0,0L:M3_L_maxIndex],xtitle=abcissaStr, $
 		xrange=[0.,M3_L_maxIndex*delta_x],$
 		ytitle=ytitleStr, $
 		yrange=[-0.01,0.04], $
 		linestyle=0,/ystyle,/xstyle,charsize=chsize
 	FOR angleIndex=1L,N_angles-1L DO BEGIN
 		OPLOT,tau[0L:M3_L_maxIndex],M12_L_modelled[angleIndex,0L:M3_L_maxIndex],linestyle=angleIndex
 	ENDFOR
 
 	;;;;;;;;
 
 	postScriptOut=0L
 	abcissaTex="\tau (m)" 	;;; tau_{1} (m)
 	ytitleTex="C^{(3)}_{L}(\tau,\tau)"		;;; (tau_{1},0)
 
	SET_PLOT,'X'
	!P.FONT = -1
	windowIndex +=1
	WINDOW,windowIndex,xsize = xwinsize,ysize = ywinsize,title='Modelled Third glint cumulant function - tau_x = tau_y',RETAIN=2
	abcissaStr=TeXtoIDL(abcissaTex,FONT=-1)
	ytitleStr=TeXtoIDL(ytitleTex,FONT=-1)
 
 	PLOT,tau[0L:M3_L_maxIndex],C12_L_modelled[0,0L:M3_L_maxIndex],xtitle=abcissaStr, $
 		xrange=[0.,M3_L_maxIndex*delta_x],$
 		ytitle=ytitleStr, $
		yrange=[-0.01,0.01], $
 		linestyle=0,/ystyle,/xstyle,charsize=chsize
 	FOR angleIndex=1L,N_angles-1L DO BEGIN
 		OPLOT,tau[0L:M3_L_maxIndex],C12_L_modelled[angleIndex,0L:M3_L_maxIndex],linestyle=angleIndex
 	ENDFOR
 
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;; Model the glint third moment functions for the various geometries,
	;;; using the simulated slope cumulant function, and the modelled glint 
	;;; moments associated with those cumulants, for comparison with the 
	;;; simulated glint second moment functions for the various geometries.
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	PRINT,'#####################################################',Format='(/A/)'
	PRINT,'Modelling the glint third moment function...',Format='(/A/)'

	;;; How much of the matrix to we compute (for debugging purposes)
	;M3_L_maxIndex = 10L 
	;M3_L_maxIndex = 20L
	M3_L_maxIndex = NN/2L-1L

	;;; Initialise the SLOPECUMULANTS COMMON block with the simulated slope cumulants,
	;;; so they can be used in the modelling of the second glint moment and cumulant functions

	kappa_xi_1 = kappa_xi_sim[0]
	kappa_xi_2 = kappa_xi_sim[1]
	kappa_xi_3 = kappa_xi_sim[2]

	;;; Initialise the GLINTMOMENTS COMMON block with the modelled glint moments,
	;;; so they can be used in the modelling of the second glint moment and cumulant functions
	mu_L_1 = mu_L1_modelled

	;M3_L_modelled = DBLARR(N_angles,NN/2L+1L,NN/2L+1L)
	M3_L_modelled = DBLARR(NN,NN,N_angles)
	M3_L_temp = DBLARR(N_angles)

	M3_L_modelled[0L,0L,*] = 1.D

	;;; Compute the modelled glint third moment function,
	;;; avoiding the Nyquist line at 
	FOR tau1=1L,M3_L_maxIndex DO BEGIN
		FOR tau2=1L,tau1-1L DO BEGIN

			PRINT,'tau1:',tau1,'   tau2:',tau2
			
			;;; The second order slope cumulants
			kappa_xi_110 = kappa_xi_2*C2_xi_sim[tau1]
			kappa_xi_101 = kappa_xi_2*C2_xi_sim[tau2]
			kappa_xi_011 = kappa_xi_2*C2_xi_sim[ABS(tau1-tau2)]

			;;; The third order slope cumulants. Since we are calculating
			;;; in the region (tau1>tau2), and the  lags for (021) and (012)
			;;; are positive by defn., 
			kappa_xi_210 =      kappa_xi_3*C3_xi_sim[0L,tau1]
			kappa_xi_120 =      kappa_xi_3*C3_xi_sim[tau1,tau1]
			kappa_xi_201 =      kappa_xi_3*C3_xi_sim[0L,tau2]
			kappa_xi_102 =      kappa_xi_3*C3_xi_sim[tau2,tau2]
			kappa_xi_021 = -1.D*kappa_xi_3*C3_xi_sim[0L,ABS(tau2-tau1)]
			kappa_xi_012 = -1.D*kappa_xi_3*C3_xi_sim[ABS(tau2-tau1),ABS(tau2-tau1)]
			kappa_xi_111 =      kappa_xi_3*C3_xi_sim[tau1,tau2]

			kappa_xi_guess = [kappa_xi_210,kappa_xi_120,kappa_xi_111]

			;;; Compute the third glint Moment function, with the slope cumulant 
			;;; functions kappa_xi_111, kappa_xi_210 and kappa_xi_120 as inputs. The 
			;;; simulated slope cumulants are passed through the COMMON block

			M3_L_temp[*] = 0.D
			thirdGlintMomentFunction2,GEOM,kappa_xi_guess,M3_L_temp
			M3_L_modelled[tau1,tau2,*] = M3_L_temp
			M3_L_modelled[tau2,tau1,*] = M3_L_temp
		ENDFOR
	ENDFOR

	;;; Plug in the values of M3_L_modelled along the tau1 and tau1=tau2 axes...
	FOR tau1=1L,M3_L_maxIndex DO BEGIN
		M3_L_modelled[tau1,0L,*] = M21_L_modelled[*,tau1]
		M3_L_modelled[tau1,tau1,*] = M12_L_modelled[*,tau1]
		M3_L_modelled[0L,tau1,*] = M12_L_modelled[*,tau1]
	ENDFOR


	;;; Compute the modelled third glint cumulant function
	;C3_L_modelled = DBLARR(N_angles,NN/2L+1L,NN/2L+1L)
	C3_L_modelled = DBLARR(NN,NN,N_angles)
	FOR geometry=0L,N_angles-1L DO BEGIN
		FOR tau1=0L,M3_L_maxIndex DO BEGIN
			FOR tau2=0L,tau1 DO BEGIN
				C3_L_modelled[tau1,tau2,geometry] = (mu_L1_modelled[geometry]*M3_L_modelled[tau1,tau2,geometry] $
					- (mu_L1_modelled[geometry]^2.D)*(M2_L_modelled[geometry,tau1] $
					+ M2_L_modelled[geometry,tau2] $
					+ M2_L_modelled[geometry,ABS(tau2-tau1)]) $
					+ 2.D*(mu_L1_modelled[geometry]^3.D))/kappa_L_modelled[geometry,2L]

				C3_L_modelled[tau2,tau1,geometry] = C3_L_modelled[tau1,tau2,geometry]
			ENDFOR
		ENDFOR
		C3_L_modelled[0L,0L,geometry] = 1.D
	ENDFOR

	HELP,M3_L_modelled
	HELP,C3_L_modelled
	temp_C3 = DBLARR(NN,NN)

	;;; Fill the rest of the third moment function array
	FOR geometry=0L,N_angles-1L DO BEGIN
		temp_C3 = M3_L_modelled[*,*,geometry]
		biCovarianceSymmetry,temp_C3,NN
		M3_L_modelled[*,*,geometry] = temp_C3
	ENDFOR

	;;; Fill the rest of the third cumulant function array
	FOR geometry=0L,N_angles-1L DO BEGIN
		temp_C3 = C3_L_modelled[*,*,geometry]
		biCovarianceSymmetry,temp_C3,NN
		C3_L_modelled[*,*,geometry] = temp_C3
	ENDFOR

 	;;;;;;;;

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;;	Glint Third Moment Function                 ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	;;--- 3D plot setup
	xwinsize=800
	ywinsize=450
	
	rotate_z=-20.
	rotate_x=10.
	z_axis=2
	ysize=1.*ywinsize
	xsize=ROUND(1.2*ysize)

	;plotLim = 16L
	;plotLim = 32L
	;plotLim = 64L
	;plotLim = NN/8-1L
	plotLim = NN/16-1L
	;plotLim = NN/2-1L

	bispecMin = MIN(DOUBLE(glintThirdMomentFunction[*,*,0]))
	bispecMax = MAX(DOUBLE(glintThirdMomentFunction[*,*,0]))
	PRINT,"glint third moment function minimum: Angle "+STRING(0)+" = ",bispecMin
	PRINT,"glint third moment function maximum: Angle "+STRING(0)+" = ",bispecMax

	FOR angle=0,N_angles-1L DO BEGIN
		!P.MULTI = 0
		windowIndex++
		window,windowIndex,xsize = 2*xsize,ysize = 2*ysize,title='Modelled Glint Third Moment Function: Angle:'+STRING(angle),RETAIN=2

		;bispecMin = MIN(DOUBLE(M3_L_modelled[0:plotLim,0:plotLim,angle]))
		;bispecMax = MAX(DOUBLE(M3_L_modelled[0:plotLim,0:plotLim,angle]))
		;PRINT,"glint third moment function minimum: Angle "+STRING(angle)+" = ",bispecMin
		;PRINT,"glint third moment function maximum: Angle "+STRING(angle)+" = ",bispecMax

		shiftBispectra = SHIFT(M3_L_modelled[*,*,angle],NN/2-1L,NN/2-1L)
		
		SURFACE,DOUBLE(shiftBispectra[NN/2-1L-plotLim:NN/2-1L+plotLim,NN/2-1L-plotLim:NN/2-1L+plotLim]),AZ=rotate_z,$
			;zrange=[bispecMin,bispecMax], $
			zrange=[-0.01,0.04], $
			/xstyle,/ystyle,/zstyle,$
			charsize=3.5,xtitle='i',ytitle='j',ZAXIS=z_axis
	ENDFOR
 
	;END
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;;   Open the output HDF file   ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	PRINT,"Open output filename: ",outFileName

	outFile = fileinfo(outFileName)
	help, fileinfo(outFileName), /structure
	PRINT, outFile

	IF (NOT outFile.EXIST) THEN BEGIN
		;;; Create and open file using SD interface
		fileID = HDF_SD_START(outFileName, /CREATE)
		;fileID = HDF_OPEN(fileName, /CREATE,/WRITE)
		PRINT, 'Created new HDF file: ',outFileName
	ENDIF ELSE BEGIN
		;;; Create and open file using SD interface
		PRINT, 'HDF file ',outFileName,' exists, opening...'
		fileID = HDF_SD_START(outFileName, /RdWr)
		;fileID = HDF_OPEN(fileName, /WRITE)
		PRINT, 'Opened HDF file ',outFileName,' for reading and writing'
	ENDELSE

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;;   Add some attributes to the file   ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	
	IF (NOT outFile.EXIST) THEN BEGIN
		PRINT,"Writing global attributes to ",outFileName
		HDF_SD_ATTRSET, fileID, 'DATE', SYSTIME()
		HDF_SD_ATTRSET, fileID, 'EXPERIMENT', 'sunglintModelCalculate.pro'
		HDF_SD_ATTRSET, fileID, 'NAME', 'Geoff Cureton'
		HDF_SD_ATTRSET, fileID, 'EMAIL ADDRESS', 'geoff.cureton@physics.org'
	ENDIF ELSE BEGIN
		HDF_SD_ATTRSET, fileID, 'DATE', SYSTIME()+" Zulu"
	ENDELSE

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;;   Add some datasets to the file   ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	   
	IF (NOT outFile.EXIST) THEN BEGIN

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;;   Save the geometry information to global attributes, and variables   ;;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

		PRINT, 'Writing geometry angles and slopes...'

		sourceAngleID  = HDF_SD_CREATE(fileID, "Solar Zenith Angles", [N_angles], /FLOAT)
		HDF_SD_ADDDATA, sourceAngleID, Solar_Zenith_Angles
		numAnglesDimID = HDF_SD_DIMGETID(sourceAngleID, 0)
		HDF_SD_DIMSET, numAnglesDimID, LABEL='Number of Angles', NAME='N_angles', UNIT='radians'
		HDF_SD_ENDACCESS, sourceAngleID 

		detectorAngleID  = HDF_SD_CREATE(fileID, "Detector Zenith Angles", [N_angles], /FLOAT)
		HDF_SD_ADDDATA, detectorAngleID,Detector_Zenith_Angles
		numAnglesDimID = HDF_SD_DIMGETID(detectorAngleID, 0)
		HDF_SD_DIMSET, numAnglesDimID, LABEL='Number of Angles', NAME='N_angles', UNIT='radians'
		HDF_SD_ENDACCESS, detectorAngleID 

		specularSlopeID  = HDF_SD_CREATE(fileID, "Specular Slopes", [N_angles], /FLOAT)
		HDF_SD_ADDDATA, specularSlopeID,Specular_Slopes
		numAnglesDimID = HDF_SD_DIMGETID(specularSlopeID, 0)
		HDF_SD_DIMSET, numAnglesDimID, LABEL='Number of Angles', NAME='N_angles', UNIT='radians'
		HDF_SD_ENDACCESS, specularSlopeID 

		minSlopeID  = HDF_SD_CREATE(fileID, "Min Slopes", [N_angles], /FLOAT)
		HDF_SD_ADDDATA, minSlopeID,Min_Slopes
		numAnglesDimID = HDF_SD_DIMGETID(minSlopeID, 0)
		HDF_SD_DIMSET, numAnglesDimID, LABEL='Number of Angles', NAME='N_angles', UNIT='radians'
		HDF_SD_ENDACCESS, minSlopeID 

		maxSlopeID  = HDF_SD_CREATE(fileID, "Max Slopes", [N_angles], /FLOAT)
		HDF_SD_ADDDATA, maxSlopeID,Max_Slopes
		numAnglesDimID = HDF_SD_DIMGETID(maxSlopeID, 0)
		HDF_SD_DIMSET, numAnglesDimID, LABEL='Number of Angles', NAME='N_angles', UNIT='radians'
		HDF_SD_ENDACCESS, maxSlopeID

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;;   Save the elevation, slope and glint moments   ;;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		
		PRINT, 'Writing the simulated slope and modelled glint moments...'

		slopeMomentID = HDF_SD_CREATE(fileID, 'Simulated Slope Moments',     [N_Moments],    /FLOAT) 
		HDF_SD_ADDDATA, slopeMomentID, slopeMoments
		numMomentsDimID = HDF_SD_DIMGETID(slopeMomentID, 0)
		HDF_SD_DIMSET, numMomentsDimID, LABEL='Number of Moments', NAME='N_moments'
		HDF_SD_ENDACCESS, slopeMomentID
		
		glintMomentID = HDF_SD_CREATE(fileID, 'Modelled Glint First Moments',     [N_angles], /FLOAT)
		HDF_SD_ADDDATA, glintMomentID, mu_L1_modelled
		numAnglesDimID = HDF_SD_DIMGETID(glintMomentID, 0)
		HDF_SD_DIMSET, numAnglesDimID, LABEL='Number of Angles', NAME='N_angles', UNIT='radians'
		HDF_SD_ENDACCESS, glintMomentID

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;;   Save the elevation, slope and glint cumulants ;;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		
		PRINT, 'Writing the simulated slope and modelled glint cumulants...'
		
		slopeCumulantID = HDF_SD_CREATE(fileID, 'Simulated Slope Cumulants',[N_Moments],    /FLOAT) 
		HDF_SD_ADDDATA, slopeCumulantID, slopeCumulants
		numCumulantsDimID = HDF_SD_DIMGETID(slopeCumulantID, 0)
		HDF_SD_DIMSET, numCumulantsDimID, LABEL='Number of Cumulants', NAME='N_moments'
		HDF_SD_ENDACCESS, slopeCumulantID
		
		glintCumulantID = HDF_SD_CREATE(fileID, 'Modelled Glint Cumulants',     [N_angles,N_Moments], /FLOAT)
		HDF_SD_ADDDATA, glintCumulantID, kappa_L_modelled
		numAnglesDimID = HDF_SD_DIMGETID(glintCumulantID, 0)
		HDF_SD_DIMSET, numAnglesDimID, LABEL='Number of Angles', NAME='N_angles', UNIT='radians'
		numCumulantsDimID = HDF_SD_DIMGETID(glintCumulantID, 1)
		HDF_SD_DIMSET, numCumulantsDimID, LABEL='Number of Cumulants', NAME='N_moments'
		HDF_SD_ENDACCESS, glintCumulantID

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;;   Set 1D wavenumber scale   ;;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		
		wavenumber = FINDGEN(N)*delta_k
		wavenumberID  = HDF_SD_CREATE(fileID, "Power Spectrum wavenumber scale", [N], /FLOAT)
		HDF_SD_ADDDATA, wavenumberID,  wavenumber
		HDF_SD_ATTRSET, wavenumberID,  'units', 'meters^{-1}'
		HDF_SD_ATTRSET, wavenumberID,  'increment', delta_k
		wavenumberDimID = HDF_SD_DIMGETID(wavenumberID, 0)
		HDF_SD_DIMSET, wavenumberDimID, LABEL='Data length', NAME='N', SCALE=wavenumber, UNIT='meters^{-1}'
		HDF_SD_ENDACCESS, wavenumberID

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;;   Set 1D lag scale   ;;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		
		length = FINDGEN(N)*delta_x
		lengthID  = HDF_SD_CREATE(fileID, "Cumulant Function 1D length scale", [N], /FLOAT)
		HDF_SD_ADDDATA, lengthID,  length
		HDF_SD_ATTRSET, lengthID,  'units', 'meters'
		HDF_SD_ATTRSET, lengthID,  'increment', delta_x
		lengthDimID = HDF_SD_DIMGETID(lengthID, 0)
		HDF_SD_DIMSET, lengthDimID, LABEL='Data length', NAME='N', SCALE=length, UNIT='meters'
		HDF_SD_ENDACCESS, lengthID

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;;   Save the Second Moment Functions   ;;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		
		PRINT, 'Writing simulated slope and modelled glint second moment functions ...'
		
		;;; Slope Second Moment Function
		
		slopeSecondMomentFunctionID = HDF_SD_CREATE(fileID, "Simulated Slope Second Moment Function"    , [N], /FLOAT)
		HDF_SD_ADDDATA, slopeSecondMomentFunctionID, slopeSecondMomentFunction
		HDF_SD_ATTRSET, slopeSecondMomentFunctionID, 'long_name', $
			'Realisation averaged slope Second Moment Function, secondMomentFunction[0 ... N-1]'
		secondMomentFunctionDimID = HDF_SD_DIMGETID(slopeSecondMomentFunctionID, 0)
		HDF_SD_DIMSET, secondMomentFunctionDimID, LABEL='Second Moment Function Data length', NAME='N', SCALE=length, UNIT='meters'
		HDF_SD_ENDACCESS, slopeSecondMomentFunctionID
		
		;;; Glint Second Moment Function 
		
		glintSecondMomentFunctionID = HDF_SD_CREATE(fileID, "Modelled Glint Second Moment Function"    , [N_angles, N], /FLOAT)
		HDF_SD_ADDDATA, glintSecondMomentFunctionID, M2_L_modelled
		HDF_SD_ATTRSET, glintSecondMomentFunctionID, 'long_name', $
			'Modelled glint Second Moment Function, secondMomentFunction[0 ... N_angles-1][0 ... N-1]'
		numAnglesDimID = HDF_SD_DIMGETID(glintSecondMomentFunctionID, 0)
		HDF_SD_DIMSET, numAnglesDimID, LABEL='Number of Angles', NAME='N_angles', UNIT='radians'
		secondMomentFunctionDimID = HDF_SD_DIMGETID(glintSecondMomentFunctionID, 1)
		HDF_SD_DIMSET, secondMomentFunctionDimID, LABEL='Second Moment Function Data length', NAME='N', SCALE=length, UNIT='meters'
		HDF_SD_ENDACCESS, glintSecondMomentFunctionID

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;;   Save the Second Cumulant Functions   ;;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		
		PRINT, 'Writing simulated slope and modelled glint second cumulant functions ...'
		
		;;; Slope Second Cumulant Function
		
		slopeSecondCumulantFunctionID = HDF_SD_CREATE(fileID, "Simulated Slope Second Cumulant Function"    , [N], /FLOAT)
		HDF_SD_ADDDATA, slopeSecondCumulantFunctionID, slopeSecondCumulantFunction
		HDF_SD_ATTRSET, slopeSecondCumulantFunctionID, 'long_name', $
			'Realisation averaged slope Second Cumulant Function, secondCumulantFunction[0 ... N-1]'
		secondCumulantFunctionDimID = HDF_SD_DIMGETID(slopeSecondCumulantFunctionID, 0)
		HDF_SD_DIMSET, secondCumulantFunctionDimID, LABEL='Second Cumulant Function Data length', NAME='N', SCALE=length, UNIT='meters'
		HDF_SD_ENDACCESS, slopeSecondCumulantFunctionID
		
		;;; Glint Second Cumulant Function
		
		glintSecondCumulantFunctionID = HDF_SD_CREATE(fileID, "Modelled Glint Second Cumulant Function"    , [N_angles, N], /FLOAT)
		HDF_SD_ADDDATA, glintSecondCumulantFunctionID, C2_L_modelled
		HDF_SD_ATTRSET, glintSecondCumulantFunctionID, 'long_name', $
			'Modelled glint Second Cumulant Function, secondCumulantFunction[0 ... N_angles-1][0 ... N-1]'
		numAnglesDimID = HDF_SD_DIMGETID(glintSecondCumulantFunctionID, 0)
		HDF_SD_DIMSET, numAnglesDimID, LABEL='Number of Angles', NAME='N_angles', UNIT='radians'
		secondCumulantFunctionDimID = HDF_SD_DIMGETID(glintSecondCumulantFunctionID, 1)
		HDF_SD_DIMSET, secondCumulantFunctionDimID, LABEL='Second Cumulant Function Data length', NAME='N', SCALE=length, UNIT='meters'
		HDF_SD_ENDACCESS, glintSecondCumulantFunctionID

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;;   Set 2D lag scale   ;;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		
		length2 = FINDGEN(NN)*delta_x
		length2ID  = HDF_SD_CREATE(fileID, "Cumulant Function 2D length scale", [NN], /FLOAT)
		HDF_SD_ADDDATA, length2ID,  length2
		HDF_SD_ATTRSET, length2ID,  'units', 'meters'
		HDF_SD_ATTRSET, length2ID,  'increment', delta_x
		length2DimID = HDF_SD_DIMGETID(length2ID, 0)
		HDF_SD_DIMSET, length2DimID, LABEL='Data length', NAME='NN', SCALE=length2, UNIT='meters'
		HDF_SD_ENDACCESS, length2ID

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;;   Save the Third Moment Functions   ;;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		
		PRINT, 'Writing simulated slope and modelled glint third moment functions ...'
		
		;;; Slope Third Moment Function
		
		slopeThirdMomentFunctionID = HDF_SD_CREATE(fileID, "Simulated Slope Third Moment Function"    , [NN,NN], /FLOAT)
		HDF_SD_ADDDATA, slopeThirdMomentFunctionID, slopeThirdMomentFunction
		HDF_SD_ATTRSET, slopeThirdMomentFunctionID, 'long_name', $
			'Realisation averaged slope Third Moment Function, ThirdMomentFunction[0 ... NN-1][0 ... NN-1]'
		thirdMomentFunctionDimID = HDF_SD_DIMGETID(slopeThirdMomentFunctionID, 0)
		HDF_SD_DIMSET, thirdMomentFunctionDimID, LABEL='Third Moment Function Data length', NAME='NN', SCALE=length2, UNIT='meters'
		thirdMomentFunctionDimID = HDF_SD_DIMGETID(slopeThirdMomentFunctionID, 1)
		HDF_SD_DIMSET, thirdMomentFunctionDimID, LABEL='Third Moment Function Data length', NAME='NN', SCALE=length2, UNIT='meters'
		HDF_SD_ENDACCESS, slopeThirdMomentFunctionID
		
		;;; Glint Third Moment Function
		
		glintThirdMomentFunctionID = HDF_SD_CREATE(fileID, "Modelled Glint Third Moment Function"    , [NN,NN, N_angles], /FLOAT)
		HDF_SD_ADDDATA, glintThirdMomentFunctionID, M3_L_modelled
		HDF_SD_ATTRSET, glintThirdMomentFunctionID, 'long_name', $
			'Modelled glint Third Moment Function, ThirdMomentFunction[0 ... NN-1][0 ... NN-1][0 ... N_angles-1]'
		thirdMomentFunctionDimID = HDF_SD_DIMGETID(glintThirdMomentFunctionID, 0)
		HDF_SD_DIMSET, thirdMomentFunctionDimID, LABEL='Third Moment Function Data length', NAME='NN', SCALE=length2, UNIT='meters'
		thirdMomentFunctionDimID = HDF_SD_DIMGETID(glintThirdMomentFunctionID, 1)
		HDF_SD_DIMSET, thirdMomentFunctionDimID, LABEL='Third Moment Function Data length', NAME='NN', SCALE=length2, UNIT='meters'
		numAnglesDimID = HDF_SD_DIMGETID(glintThirdMomentFunctionID, 2)
		HDF_SD_DIMSET, numAnglesDimID, LABEL='Number of Angles', NAME='N_angles', UNIT='radians'
		HDF_SD_ENDACCESS, glintThirdMomentFunctionID

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;;   Save the Third Cumulant Functions   ;;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		
		PRINT, 'Writing simulated slope and modelled third Cumulant functions ...'
		
		;;; Slope Third Cumulant Function
		
		slopeThirdCumulantFunctionID = HDF_SD_CREATE(fileID, "Simulated Slope Third Cumulant Function"    , [NN,NN], /FLOAT)
		HDF_SD_ADDDATA, slopeThirdCumulantFunctionID, slopeThirdCumulantFunction
		HDF_SD_ATTRSET, slopeThirdCumulantFunctionID, 'long_name', $
			'Realisation averaged slope Third Cumulant Function, ThirdCumulantFunction[0 ... NN-1][0 ... NN-1]'
		thirdCumulantFunctionDimID = HDF_SD_DIMGETID(slopeThirdCumulantFunctionID, 0)
		HDF_SD_DIMSET, thirdCumulantFunctionDimID, LABEL='Third Cumulant Function Data length', NAME='NN', SCALE=length2, UNIT='meters'
		thirdCumulantFunctionDimID = HDF_SD_DIMGETID(slopeThirdCumulantFunctionID, 1)
		HDF_SD_DIMSET, thirdCumulantFunctionDimID, LABEL='Third Cumulant Function Data length', NAME='NN', SCALE=length2, UNIT='meters'
		HDF_SD_ENDACCESS, slopeThirdCumulantFunctionID
		
		;;; Glint Third Cumulant Function 
		
		glintThirdCumulantFunctionID = HDF_SD_CREATE(fileID, "Modelled Glint Third Cumulant Function"    , [NN,NN, N_angles], /FLOAT)
		HDF_SD_ADDDATA, glintThirdCumulantFunctionID, C3_L_modelled
		HDF_SD_ATTRSET, glintThirdCumulantFunctionID, 'long_name', $
			'Modelled glint Third Cumulant Function, ThirdCumulantFunction[0 ... NN-1][0 ... NN-1][0 ... N_angles-1]'
		thirdCumulantFunctionDimID = HDF_SD_DIMGETID(glintThirdCumulantFunctionID, 0)
		HDF_SD_DIMSET, thirdCumulantFunctionDimID, LABEL='Third Cumulant Function Data length', NAME='NN', SCALE=length2, UNIT='meters'
		thirdCumulantFunctionDimID = HDF_SD_DIMGETID(glintThirdCumulantFunctionID, 1)
		HDF_SD_DIMSET, thirdCumulantFunctionDimID, LABEL='Third Cumulant Function Data length', NAME='NN', SCALE=length2, UNIT='meters'
		numAnglesDimID = HDF_SD_DIMGETID(glintThirdCumulantFunctionID, 2)
		HDF_SD_DIMSET, numAnglesDimID, LABEL='Number of Angles', NAME='N_angles', UNIT='radians'
		HDF_SD_ENDACCESS, glintThirdCumulantFunctionID

	ENDIF ELSE BEGIN

	ENDELSE

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;;   Close the output HDF file   ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	PRINT, 'Write Operation Completed'
	HDF_SD_END, fileID
	PRINT, '*** File Closed ***'
	PRINT, ''

END
