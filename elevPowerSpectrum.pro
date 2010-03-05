;;; Function to calculate the elevation power spectrum (nonlinear)

PRO gaussian_elev_spectrum,SCALE,POWER,NLCOUPLING
	print,'Calculating the Gaussian power spectrum...', Format='(/A)'

	elevStDev = 0.13D ; meters
	elevVar = elevStDev*elevStDev
	slopeStDev = 0.2121D
	slopeVar = slopeStDev*slopeStDev
	elevCorr = SQRT(2.0D)*elevStDev/slopeStDev ; meters

	PRINT,''
	PRINT,'Theoretical elevation stdev:   ',elevStDev,' meters',Format='(A,F12.6,A)'
	PRINT,'Theoretical elevation variance: ',elevVar,' meters^{2}',Format='(A,F12.6,A)'
	PRINT,'Theoretical slope stdev:       ',slopeStDev,Format='(A,F12.6,A)'
	PRINT,'Theoretical slope variance:       ',slopeVar,Format='(A,F12.6,A)'
	PRINT,'Elevation correlation length:  ',elevCorr,' meters',Format='(A,F12.6,A)'
	PRINT,''

	N = SCALE.N
	dk = SCALE.delta_k

	k = DBLARR(LONG(N))
    k = DINDGEN(LONG(N))*dk;
	exp_arg = -(elevCorr*elevCorr*k*k)/4.0D;
	;;; Causes arithmetic underflow! OK though
	power = elevCorr*elevVar*exp(exp_arg)/SQRT(!DPI)
	power[0L]=0.D

;	xwinsize=800
;	ywinsize=450
;	WINDOW,15,xsize = xwinsize,ysize = ywinsize,title='Elevation Power Spectrum',RETAIN=2
;	PLOT,k,0.5*power,xrange=[0.01D,10.0],xtitle='k',/xstyle,$
;		ytitle='power',linestyle=0,/ystyle,charsize=chsize

END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; Function to calculate a Phillips elevation power spectrum,
;;; and the desired phase relationships

PRO phillips_elev_spectrum,SCALE,POWER,NLCOUPLING
	
	COMMON DEBUGGING, debugVar,windowIndex
	
	PRINT,'Calculating the Phillips power spectrum...', Format='(/A)'

	N = SCALE.N
	dk = SCALE.delta_k
	k_max = SCALE.k_max

;	PRINT,'k_max = ',k_max,' meters^{-1}', Format='(A25,F20.10,A)'
;	PRINT,'dk = ',dk,' meters^{-1}', Format='(A25,F20.10,A)'
;	PRINT,'Nyquist Wavenumber = ',DOUBLE(N/2L)*dk,' meters^{-1}',Format='(A25,F20.12,A)'

	B = 0.005D		;;; Dimensionless constant
	g = 9.81D		;;; Acceleration due to gravity in m * s^{-2}
	v = 5.D			;;; Wind velocity in m * s^{-1}
	gammaT = 0.073D	;;; Water surface tension in N * m^{-1}
	rho = 1000.D	;;; Water density in kg * m^{3}

	kc0 = g/(v*v)	;;; Lower wavenumber cutoff for wind waves
	k_gamma = MIN(SQRT(rho*g/gammaT),k_max)	;;; Upper wavenumber cutoff for the capilliary wave region

	elevStDev = SQRT(0.25D*B*(kc0^(-2.D)-k_gamma^(-2.D)))
	elevVar = 0.25D*B*(kc0^(-2.D)-k_gamma^(-2.D))
	slopeStDev = SQRT(0.5D*B*ALOG(k_gamma/kc0))
	slopeVar = 0.5D*B*ALOG(k_gamma/kc0)

	k = DBLARR(LONG(N))
    k = DINDGEN(LONG(N))*dk

	PRINT,''
	PRINT,'Theoretical elevation stdev:   ',elevStDev,' meters',Format='(A,F12.6,A)'
	PRINT,'Theoretical elevation variance: ',elevVar,' meters^{2}',Format='(A,F12.6,A)'
	PRINT,'Theoretical slope stdev:       ',slopeStDev,Format='(A,F12.6,A)'
	PRINT,'Theoretical slope variance:       ',slopeVar,Format='(A,F12.6,A)'
	PRINT,''

	PRINT,'kc0 cutoff is: ',kc0,' meters^{-1}',Format='(A,F12.6,A)'
	PRINT,'1.5 times the kc0 cutoff is: ',1.5D*kc0,' meters^{-1}',Format='(A,F12.6,A)'
	PRINT,'k_gamma cutoff is: ',k_gamma,' meters^{-1}',Format='(A,F12.6,A/)'

	RETURN

;	PRINT,"k= ",k[0L:9L],FORMAT='(A,9F10.6/)'
	
	;;; Determine the indices of components in the wavenumber range for the free waves
;	PRINT,"Indices between kc0 and 1.5*kc0 (free1 and free2, the primary power indicies)..."
	sourceIndex=WHERE((k[0:N-1] GT kc0)*(k[0:N-1] LT 1.5D*kc0))
;	PRINT,"sourceindex...",sourceIndex
;	PRINT,"sourceindex...",sourceIndex," at wavenumber ",k[sourceIndex]

	;;; Define the structure containing the phase relationships between the free
	;;; and bound waves
	NLCOUPLING = {Nbound:N_ELEMENTS(sourceIndex),bound:LONARR(N_ELEMENTS(sourceIndex)),$
		free1:LONARR(N_ELEMENTS(sourceIndex)),free2:LONARR(N_ELEMENTS(sourceIndex))}

;	PRINT,"There are ",NLCOUPLING.Nbound," bound-free pairs",Format='(A,I3,A/)'

	;;; Determine the indices of the bound waves
;	PRINT,"Indices between 2*kc0 and 3*kc0 (bound, the coupled power indicies)..."
	coupleIndex = 2L*WHERE((k[0:N-1] GT kc0)*(k[0:N-1] LT 1.5D*kc0))
;	PRINT,"coupleIndex...",coupleIndex," at wavenumber ",k[coupleIndex]

	NLCOUPLING.free1 = sourceIndex
	NLCOUPLING.free2 = sourceIndex
	NLCOUPLING.bound = coupleIndex

	;; Compute the total power S(k)
	totalPower = DBLARR(N)
	totalPower[0] = 0.D
	;totalPower[1:N-1] = (k[1:N-1] GT kc0)*(k[1:N-1] LT k_gamma)*B/(2.D*(k[1:N-1]*k[1:N-1]*k[1:N-1]))
	totalPower[1:N-1] = (k[1:N-1] GT kc0)*(k[1:N-1] LT k_gamma)*B/(2.D*(k[1:N-1]*k[1:N-1]*k[1:N-1]*k[1:N-1]))

	;;; Set the bound power at the bound wavenumbers to 35% of the total power
	;;; at those wavenumbers
	nlPower = DBLARR(N)
	nlPower[0L] = 0.D
	;nlPower[coupleIndex] = 0.35D*B/(2.D*(k[coupleIndex]*k[coupleIndex]*k[coupleIndex]))
	nlPower[coupleIndex] = 0.35D*B/(2.D*(k[coupleIndex]*k[coupleIndex]*k[coupleIndex]*k[coupleIndex]))
	POWER.nlPower = nlPower

	;;; Define the primary power spectrum
	primaryPower = DBLARR(N)
	primaryPower = totalPower - nlPower
	POWER.primaryPower = primaryPower

END
