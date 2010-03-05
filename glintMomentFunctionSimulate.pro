;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; Subroutine to simulate a number of glint realisations, generating a
;;; simulated glint first moment, and glint second moment function

	;;; SCALE: Structure containing the various values of 
	;;;		spatial and wavenumber increments
	;;; GEOM: Structure containing angular dependent parameters
	;;;		for a range of geometries.
	;;; SLOPE : Structure containing the slope moments and 
	;;; 		second slope moment function
	;;; GLINT: Structure containing the first glint moment and 
	;;;		second glint moment function for each 
	;;;		geometry defined in GEOM
	;;; POWER: Structure containing the primary and nonlinear 
	;;; 		elevation power spectra, and the components 
	;;; 		connected by the three way coupling
	;;; Nruns: The number of runs in the simulation

PRO glintMomentFunctionSimulate,SCALE,GEOM,SLOPE,GLINT,POWER,NLCOUPLING,Nruns,nlSwitch

	COMMON DEBUGGING, debugVar,windowIndex

	IF(windowIndex GT 10) THEN windowIndex -= 21

	PRINT,'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>',FORMAT='(/A)'
	PRINT,"Inside glintMomentFunctionSimulate..."
	PRINT,'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>',FORMAT='(A/)'

	d2r = !DPI/180.D
	r2d = 180.D/!DPI

	;;; The scale parameters
	N = SCALE.N
	NN = SCALE.NN
	delta_x = SCALE.delta_x
	x_max = SCALE.x_max
	k_max = SCALE.k_max
	delta_k = SCALE.delta_k

	x = FINDGEN(N)*delta_x	

	k_N = DOUBLE(N/2L)*delta_k
	k = DBLARR(N)
	k = DINDGEN(N)*delta_k

	tau = DBLARR(N)
	tau = DINDGEN(N)*delta_x

	;;; The geometry parameters
	N_angles = GEOM.N_angles
	slopeMin = GEOM.xi_min[GLINT.geomIndex]
	slopeMax = GEOM.xi_max[GLINT.geomIndex]

	;;; Set some plot limits based on the spectrum used (in terms of
    ;;;	the Nyquist wavenumber k_N
	CASE POWER.spectrumType OF
	'gaussian': BEGIN
			k_PlotMin = 0.D
			k_PlotMax = 0.1D
		END
	'phillips': BEGIN
			k_PlotMin = 0.D
			k_PlotMax = 0.02D
		END
	ENDCASE

	;;; Define the total elevation power based on the primary and nonlinear components
	totalElevPower = DBLARR(N)
	totalElevPower = POWER.primaryPower + POWER.nlPower
	PRINT,"Total Elevation stdev from power vector:    ",SQRT(TOTAL(totalElevPower)*delta_k)," meters",FORMAT='(A,F10.6,A)'
	PRINT,"Total Elevation variance from power vector: ",TOTAL(totalElevPower)*delta_k," meters^{2}",FORMAT='(A,F10.6,A)'

	;;; Define the slope primary, nonlinear and total power from the elevation equivalents
	totalSlopePower = DBLARR(N)
	totalSlopePower = k*k*totalElevPower
	primarySlopePower = DBLARR(N)
	primarySlopePower = k*k*POWER.primaryPower
	nlSlopePower = DBLARR(N)
	nlSlopePower = k*k*POWER.nlPower
	PRINT,"Total Slope stdev from power vector: ",SQRT(TOTAL(totalSlopePower)*delta_k),FORMAT='(A,F10.6)'
	PRINT,"Total Slope variance from power vector: ",TOTAL(totalSlopePower)*delta_k,FORMAT='(A,F10.6/)'

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;; Compute the elevation amplitude, phase and spectrum       ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	xwinsize=800
	ywinsize=450
	!P.MULTI=0

	;windowIndex++
	;WINDOW,windowIndex,xsize = xwinsize,ysize = ywinsize,title='Elevation Power Spectrum',RETAIN=2
	;PLOT,k/k_N,totalElevPower,xrange=[k_PlotMin,k_PlotMax],yrange=[0.D,1.2D*MAX(totalElevPower)], $
		;xtitle='k/k_N',/xstyle,ytitle='Total Elevation Power',linestyle=0,/ystyle,charsize=chsize
	;OPLOT,k/k_N,totalElevPower,PSYM=2
	;OPLOT,k/k_N,POWER.primaryPower,linestyle=1
	;OPLOT,k/k_N,POWER.nlPower,linestyle=2
	;OPLOT,k[NLCOUPLING.free1]/k_N,totalElevPower[NLCOUPLING.free1],PSYM=6

	totalElevAmplitude = DBLARR(N)
	totalElevAmplitude = SQRT(0.5D*totalElevPower*delta_k)
	totalElevAmplitude[N/2L+1L:N-1L] = REVERSE(totalElevAmplitude[1L:N/2L-1L])

	primaryElevAmplitude = DBLARR(N)
	primaryElevAmplitude = SQRT(0.5D*POWER.primaryPower*delta_k)
	primaryElevAmplitude[N/2L+1L:N-1L] = REVERSE(primaryElevAmplitude[1L:N/2L-1L])

	nlElevAmplitude = DBLARR(N)
	nlElevAmplitude = SQRT(0.5D*POWER.nlPower*delta_k)
	nlElevAmplitude[N/2L+1L:N-1L] = REVERSE(nlElevAmplitude[1L:N/2L-1L])

	PRINT,"Total Elevation stdev from amplitude vector: ",SQRT(TOTAL(totalElevAmplitude^2.D))
	PRINT,"Total Elevation Variance from amplitude vector: ",TOTAL(totalElevAmplitude^2.D)

	seed = 30L
	totalElevPhase = DBLARR(N)
	totalElevPhase = RANDOMU(seed,N)*2.D*!DPI - !DPI

	totalElevSpectrum = totalElevAmplitude*DCOMPLEX(COS(totalElevPhase),SIN(totalElevPhase))
	totalElevSpectrum[N/2L+1L:N-1L] = CONJ(REVERSE(totalElevSpectrum[1L:N/2L-1L]))

	totalElevation = DCOMPLEXARR(N)
	totalElevation = FFT(totalElevSpectrum,/INVERSE,/DOUBLE)
	PRINT,"Total Elevation stdev from surface: ",SQRT(VARIANCE(DOUBLE(totalElevation),/DOUBLE))
	PRINT,"Total Elevation Variance from surface: ",VARIANCE(DOUBLE(totalElevation),/DOUBLE)

	elevSkewness = 0.D
	elevSkewness = SKEWNESS(DOUBLE(totalElevation),/DOUBLE)
	PRINT,"Elevation skewness from surface: ",elevSkewness*(VARIANCE(DOUBLE(totalElevation),/DOUBLE))^1.5D

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;; Compute the total slope amplitude, phase and spectrum     ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	xwinsize=800
	ywinsize=450
	!P.MULTI=0

	;windowIndex++
	;WINDOW,windowIndex,xsize = xwinsize,ysize = ywinsize,title='Slope Power Spectrum',RETAIN=2
	;PLOT,k/k_N,totalSlopePower,xrange=[k_PlotMin,k_PlotMax],yrange=[0.D,1.2D*MAX(totalSlopePower)], $
		;xtitle='k',/xstyle,ytitle='Total Slope Power',linestyle=0,/ystyle,charsize=chsize
	;OPLOT,k/k_N,totalSlopePower,PSYM=2
	;OPLOT,k/k_N,primarySlopePower,linestyle=1
	;OPLOT,k/k_N,nlSlopePower,linestyle=2
	;OPLOT,k[NLCOUPLING.free1]/k_N,totalSlopePower[NLCOUPLING.free1],PSYM=6
	;windowIndex--
	;windowIndex--

	totalSlopeAmplitude = DBLARR(N)
	totalSlopeAmplitude = SQRT(0.5D*totalSlopePower*delta_k)
	totalSlopeAmplitude[N/2L+1L:N-1L] = REVERSE(totalSlopeAmplitude[1L:N/2L-1L])
	totalSlopeSpectrum = totalSlopeAmplitude*DCOMPLEX(-SIN(totalElevPhase),COS(totalElevPhase))
	totalSlopeSpectrum[N/2L+1L:N-1L] = CONJ(REVERSE(totalSlopeSpectrum[1L:N/2L-1L]))
	PRINT,""
	PRINT,"Slope stdev from amplitude vector: ",SQRT(TOTAL(totalSlopeAmplitude^2.D))
	PRINT,"Slope Variance from amplitude vector: ",TOTAL(totalSlopeAmplitude^2.D)
	totalSlopeSurface = DCOMPLEXARR(N)
	totalSlopeSurface = FFT(totalSlopeSpectrum,/INVERSE,/DOUBLE)
	PRINT,"Slope stdev from total surface: ",SQRT(VARIANCE(DOUBLE(totalSlopeSurface),/DOUBLE))
	PRINT,"Slope Variance from total surface: ",VARIANCE(DOUBLE(totalSlopeSurface),/DOUBLE)
	slopeSkewness = 0.D
	slopeSkewness = SKEWNESS(DOUBLE(totalSlopeSurface),/DOUBLE)
	PRINT,"Slope skewness from surface: ",slopeSkewness*(VARIANCE(DOUBLE(totalSlopeSurface),/DOUBLE))^1.5D

	primarySlopeAmplitude = DBLARR(N)
	primarySlopeAmplitude = SQRT(0.5D*primarySlopePower*delta_k)
	primarySlopeAmplitude[N/2L+1L:N-1L] = REVERSE(primarySlopeAmplitude[1L:N/2L-1L])
	primarySlopeSpectrum = primarySlopeAmplitude*DCOMPLEX(-SIN(totalElevPhase),COS(totalElevPhase))
	primarySlopeSpectrum[N/2L+1L:N-1L] = CONJ(REVERSE(primarySlopeSpectrum[1L:N/2L-1L]))
	primarySlopeSurface = DCOMPLEXARR(N)
	primarySlopeSurface = FFT(primarySlopeSpectrum,/INVERSE,/DOUBLE)
	PRINT,"Slope stdev from primary surface: ",SQRT(VARIANCE(DOUBLE(primarySlopeSurface),/DOUBLE))
	PRINT,"Slope Variance from primary surface: ",VARIANCE(DOUBLE(primarySlopeSurface),/DOUBLE)

	nlSlopeAmplitude = DBLARR(N)
	nlSlopeAmplitude = SQRT(0.5D*nlSlopePower*delta_k)
	nlSlopeAmplitude[N/2L+1L:N-1L] = REVERSE(nlSlopeAmplitude[1L:N/2L-1L])
	nlSlopeSpectrum = nlSlopeAmplitude*DCOMPLEX(-SIN(totalElevPhase),COS(totalElevPhase))
	nlSlopeSpectrum[N/2L+1L:N-1L] = CONJ(REVERSE(nlSlopeSpectrum[1L:N/2L-1L]))
	nlSlopeSurface = DCOMPLEXARR(N)
	nlSlopeSurface = FFT(nlSlopeSpectrum,/INVERSE,/DOUBLE)
	PRINT,"Slope stdev from nl surface: ",SQRT(VARIANCE(DOUBLE(nlSlopeSurface),/DOUBLE))
	PRINT,"Slope Variance from nl surface: ",VARIANCE(DOUBLE(nlSlopeSurface),/DOUBLE)

;	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;	;;; Plot a realisation of the elevation and slope to check    ;;;
;	;;; that they seem to have the correct relationship           ;;;
;	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;	xsize=450
;	ysize=ROUND(2.4*xsize)
;	; 	!P.MULTI=0
;	!P.MULTI = [0,1,4]
;
;	windowIndex++
;	WINDOW,windowIndex,xsize = 2*xsize,ysize = ysize,title='Elevation and Slope Realisations',RETAIN=2
;	PLOT,x,elevation,xrange=[0.,5.],xtitle='x',/xstyle
;	OPLOT,x,IMAGINARY(elevation),LINESTYLE=2
;
;	PLOT,x,slope,xrange=[0.,5.],xtitle='x',/xstyle
;	OPLOT,x,IMAGINARY(slope),LINESTYLE=2
;	OPLOT,x,DERIV(x,elevation),LINESTYLE=3
;
;	PLOT,k/k_N,2.D*((ABS(FFT(elevation,/DOUBLE)))^2.D)/delta_k,$
;		xrange=[k_PlotMin,k_PlotMax],xtitle='x',yrange=[-0.025,0.025],/xstyle
;	OPLOT,k/k_N,DOUBLE(FFT(elevation,/DOUBLE)), LINESTYLE=2
;	OPLOT,k/k_N,IMAGINARY(FFT(slope,/DOUBLE)),LINESTYLE=3
;	OPLOT,k/k_N,2.D*((DOUBLE(FFT(elevation,/DOUBLE)))^2.D + (IMAGINARY(FFT(elevation,/DOUBLE)))^2.D)/delta_k, $
;		LINESTYLE=4
;
;	PLOT,k/k_N,2.D*((ABS(FFT(slope,/DOUBLE)))^2.D)/delta_k,$
;		xrange=[k_PlotMin,k_PlotMax],xtitle='x',yrange=[-0.025,0.025],/xstyle
;	OPLOT,k/k_N,DOUBLE(FFT(slope,/DOUBLE)), LINESTYLE=2
;	OPLOT,k/k_N,IMAGINARY(FFT(slope,/DOUBLE)),LINESTYLE=3
;	OPLOT,k/k_N,2.D*((DOUBLE(FFT(slope,/DOUBLE)))^2.D + (IMAGINARY(FFT(slope,/DOUBLE)))^2.D)/delta_k, $
;		LINESTYLE=4

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;; Define the various quantities for the calculation of the  ;;;
	;;; bispectrum and the component power spectra                ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	NN2 = LONG(FIX(NN/2L))
	NN4 = LONG(FIX(NN/4L))
	PRINT,"NN = ",NN,FORMAT='(/A,I8)'

	kx = DINDGEN(NN)*delta_k
	ky = DINDGEN(NN)*delta_k
	tau_x = DINDGEN(NN)*delta_x
	tau_y = DINDGEN(NN)*delta_x

;	elevBispectrum =     DCOMPLEXARR(NN,NN)
;	elevComponentPower = DBLARR(NN,NN)
;	elevSumPower =       DBLARR(NN,NN)
;	elevBicoherence =    DBLARR(NN,NN)

	slopeMean = 0.D
	slopeVariance = 0.D
	slopeSkewness = 0.D

	slopeFirstCumulant = 0.D
	slopeSecondCumulant = 0.D
	slopeThirdCumulant = 0.D

	slopeBispectrum =     DCOMPLEXARR(NN,NN)
	slopeComponentPower = DBLARR(NN,NN)
	slopeSumPower =       DBLARR(NN,NN)
	slopeBicoherence =    DBLARR(NN,NN)

	glintBispectrum =     DCOMPLEXARR(NN,NN)
	glintComponentPower = DBLARR(NN,NN)
	glintSumPower =       DBLARR(NN,NN)
	glintBicoherence =    DBLARR(NN,NN)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	temp = DCOMPLEXARR(N)
	sunGlint = DBLARR(N)
	sunglintPower = DBLARR(N)

	totalSlopeAvgPower = DBLARR(N)
	primarySlopeAvgPower = DBLARR(N)
	nlSlopeAvgPower = DBLARR(N)

	tempSecondMomentFunction = DBLARR(N/2L+1L)

	mySeed = 30L
	SLOPE.moments[*] = 0.D
	SLOPE.secondMomentFunction[*] = 0.D
	SLOPE.thirdMomentFunction[*,*] = 0.D
	GLINT.firstMoment[GLINT.geomIndex] = 0.D
	GLINT.secondMomentFunction[GLINT.geomIndex,*] = 0.D
	GLINT.thirdMomentFunctionGeom[GLINT.geomIndex].thirdMomentFunction[*,*] = 0.D

	PRINT,'Commencing',Nruns,' runs for geometry ',GLINT.geomIndex

	FOR t = 1L,Nruns DO BEGIN

;		mySeed += 1L

		;;; Compute the independent phases for this realisation
		primaryElevPhase = RANDOMU(mySeed,N)*2.D*!DPI - !DPI
		nlElevPhase = RANDOMU(mySeed,N)*2.D*!DPI - !DPI

		;;; Calculate the elevation spectrum for the independent phases
;		primaryElevSpectrum = primaryElevAmplitude*DCOMPLEX(COS(elevPhase),SIN(elevPhase))
;		primaryElevSpectrum += 0.00001*MAX(totalElevAmplitude)*RANDOMN(seed,N)
;		primaryElevSpectrum[N/2L+1L:N-1L] = CONJ(REVERSE(primaryElevSpectrum[1L:N/2L-1L]))

		;;; Calculate the slope spectrum for the independent phases
		primarySlopeSpectrum = primarySlopeAmplitude*DCOMPLEX(-SIN(primaryElevPhase),COS(primaryElevPhase))
;		primarySlopeSpectrum += 0.00001*MAX(totalSlopeAmplitude)*RANDOMN(seed,N)
		primarySlopeSpectrum[N/2L+1L:N-1L] = CONJ(REVERSE(primarySlopeSpectrum[1L:N/2L-1L]))

		;;; Calculate the elevation spectrum for the coupled phases
;		nlElevSpectrum = nlElevAmplitude*DCOMPLEX(COS(nlElevPhase),SIN(nlElevPhase))
;		nlElevSpectrum[N/2L+1L:N-1L] = CONJ(REVERSE(nlElevSpectrum[1L:N/2L-1L]))

		;;; Apply the phase correlations between the free and bound wavenumbers for the nonlinear
		;;; component
		IF (nlSwitch EQ 1) THEN $
			nlElevPhase[NLCOUPLING.bound] = primaryElevPhase[NLCOUPLING.free1] + primaryElevPhase[NLCOUPLING.free2]

		;;; Calculate the slope spectrum for the coupled phases
		nlSlopeSpectrum = nlSlopeAmplitude*DCOMPLEX(-SIN(nlElevPhase),COS(nlElevPhase))
		nlSlopeSpectrum[N/2L+1L:N-1L] = CONJ(REVERSE(nlSlopeSpectrum[1L:N/2L-1L]))

		;;; Add the primary and nonlinear spectra
;		totalElevSpectrum = primaryElevSpectrum + nlElevSpectrum
;		totalSlopeSpectrum = primarySlopeSpectrum + nlSlopeSpectrum

		;;; Compute specific realisation of the elevation and slope
;		elevation = FFT(totalElevSpectrum,/INVERSE,/DOUBLE)
		primarySlopeSurface = FFT(primarySlopeSpectrum,/INVERSE,/DOUBLE)
		nlSlopeSurface = FFT(nlSlopeSpectrum,/INVERSE,/DOUBLE)
		totalSlopeSurface = primarySlopeSurface + nlSlopeSurface

		;;; Compute the average slope power spectrum
		primarySlopeAvgPower += ABS(FFT(primarySlopeSurface,/DOUBLE))^2.D
		nlSlopeAvgPower += ABS(FFT(nlSlopeSurface,/DOUBLE))^2.D
		totalSlopeAvgPower += ABS(FFT(totalSlopeSurface,/DOUBLE))^2.D

		;;; Generate glint realisation
		sunGlint = (REAL_PART(totalSlopeSurface) GT slopeMin)*$
			(REAL_PART(totalSlopeSurface) LT slopeMax)

		;;; Count the number of elements that are equal to unity
;		PRINT,'Nonzero sunglint elements: ',N_ELEMENTS(WHERE(sunGlint EQ 1.0))

		;;; Check if all sunGlint elements vanish
		result=WHERE(sunGlint,NCOMPLEMENT=vanishCount)

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;; The Slope Stuff
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;			PRINT,'Successful glint realisation...',Format='(A/)'

		;;; Compute the slope moments
		slopeMean += MEAN(DOUBLE(totalSlopeSurface),/DOUBLE)
		slopeVariance += VARIANCE(DOUBLE(totalSlopeSurface),/DOUBLE)
		slopeSkewness += SKEWNESS(DOUBLE(totalSlopeSurface),/DOUBLE)

		SLOPE.moments[0L] += TOTAL(DOUBLE(totalSlopeSurface))/DOUBLE(N) 
		SLOPE.moments[1L] += TOTAL(DOUBLE(totalSlopeSurface)^2.D)/DOUBLE(N)
		SLOPE.moments[2L] += TOTAL(DOUBLE(totalSlopeSurface)^3.D)/DOUBLE(N)

		;;; Compute the slope second moment function
		temp = FFT(totalSlopeSurface,/DOUBLE)
		tempSecondMomentFunction[0L:N/2L] = (REAL_PART(FFT(ABS(temp)^2.D,/INVERSE)))[0L:N/2L]
		tempSecondMomentFunction /= tempSecondMomentFunction[0L]
		SLOPE.secondMomentFunction[0L:N/2L] += tempSecondMomentFunction[0L:N/2L]


		;;; Calculate the average bispectrum (for the reduced domain)
		FOR j=0L,NN4 DO BEGIN
			FOR i=j,NN2-j DO BEGIN
				slopeBispectrum[i,j] += temp[i]*temp[j]*CONJ(temp[i+j])
			ENDFOR
		ENDFOR

		;;; Calculate the average component power (for the reduced domain)
		FOR j=0L,NN4 DO BEGIN
			FOR i=j,NN2-j DO BEGIN
				slopeComponentPower[i,j] += (ABS(temp[i]*temp[j]))^2.D
			ENDFOR
		ENDFOR

		;;; Calculate the average sum power (for the reduced domain)
		FOR j=0L,NN4 DO BEGIN
			FOR i=j,NN2-j DO BEGIN
				slopeSumPower[i,j] += (ABS(temp[i+j]))^2.D
			ENDFOR
		ENDFOR

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;;; The Glint Stuff
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;

		IF(vanishCount EQ N) THEN BEGIN
			PRINT,'Glint realisation gave zero glints',Format='(A/)'
			;t--
		ENDIF ELSE BEGIN

;			sunGlint = FFT(((ABS(MAX(ABS(FFT(sunglint,/DOUBLE))^2.D))*EXP(-(k/(0.8*k_N))^2.D))/ $
;				(ABS(MAX(ABS(FFT(sunglint,/DOUBLE))^2.D))*EXP(-(k/(0.8*k_N))^2.D))[0L])*FFT(sunGlint,/DOUBLE),/INVERSE,DOUBLE)

			;;; Calculate the glint first moment
			PRINT,'Run ',t;,'  seed',mySeed
			GLINT.firstMoment[GLINT.geomIndex] += TOTAL(sunGlint)/DOUBLE(N)

			;;; Compute the sunglint power spectrum
			sunglintPower += ABS(FFT(sunGlint,/DOUBLE))^2.D

			;;; Compute the glint second moment function
			temp = FFT(sunGlint,/DOUBLE)
			tempSecondMomentFunction[0L:N/2L] = (REAL_PART(FFT(ABS(temp)^2.D,/INVERSE)))[0L:N/2L]
			tempSecondMomentFunction /= tempSecondMomentFunction[0L]
			GLINT.secondMomentFunction[GLINT.geomIndex,0L:N/2L] += tempSecondMomentFunction[0L:N/2L]

			;;; Calculate the average bispectrum (for the reduced domain)
			FOR j=0L,NN4 DO BEGIN
				FOR i=j,NN2-j DO BEGIN
					glintBispectrum[i,j] += temp[i]*temp[j]*CONJ(temp[i+j])
				ENDFOR
			ENDFOR

			;;; Calculate the average component power (for the reduced domain)
			FOR j=0L,NN4 DO BEGIN
				FOR i=j,NN2-j DO BEGIN
					glintComponentPower[i,j] += (ABS(temp[i]*temp[j]))^2.D
				ENDFOR
			ENDFOR

			;;; Calculate the average sum power (for the reduced domain)
			FOR j=0L,NN4 DO BEGIN
				FOR i=j,NN2-j DO BEGIN
					glintSumPower[i,j] += (ABS(temp[i+j]))^2.D
				ENDFOR
			ENDFOR

		ENDELSE

	ENDFOR

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;; Compute the slope estimators, including the average power ;;;
	;;; of the slope components, and the bispectrum 			  ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	
	;;; Compute the estimators of the IDL moments
	slopeMean /= DOUBLE(Nruns)
	slopeVariance /= DOUBLE(Nruns)
	slopeSkewness /= DOUBLE(Nruns)

	SLOPE.moments /= DOUBLE(Nruns)
	SLOPE.secondMomentFunction /= DOUBLE(Nruns)
	SLOPE.secondMomentFunction /= SLOPE.secondMomentFunction[0L]
	primarySlopeAvgPower /= DOUBLE(Nruns)
	nlSlopeAvgPower /= DOUBLE(Nruns)
	totalSlopeAvgPower /= DOUBLE(Nruns)

	slopeBispectrum /= DOUBLE(Nruns)
	slopeComponentPower /= DOUBLE(Nruns)
	slopeSumPower /= DOUBLE(Nruns)
	FOR j=0L,NN4 DO BEGIN
		FOR i=j,NN2-j DO BEGIN
			IF (SQRT(slopeComponentPower[i,j])*SQRT(slopeSumPower[i,j]) GT 10.D^(-12.D)) THEN BEGIN
				slopeBicoherence[i,j] = ABS(slopeBispectrum[i,j])/(SQRT(slopeComponentPower[i,j])*SQRT(slopeSumPower[i,j]))
			ENDIF ELSE BEGIN
				slopeBicoherence[i,j] = 0.D
			ENDELSE
		ENDFOR
	ENDFOR

	PRINT,'IDL Slope Moments:   ',slopeMean,slopeVariance,slopeSkewness*(slopeVariance^1.5D),FORMAT='(A,3E16.7)'
	PRINT,'    Slope Moments:   ',SLOPE.moments,FORMAT='(A,3E16.7)'
;	PRINT,'    Slope Cumulants: ',slopeFirstCumulant,slopeSecondCumulant,slopeThirdCumulant,FORMAT='(A,3E16.7)'
	PRINT,''

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;; Compute the glint estimators, including bispectrum 		  ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	GLINT.firstMoment[GLINT.geomIndex] /= DOUBLE(Nruns)
	GLINT.secondMomentFunction[GLINT.geomIndex,*] /= DOUBLE(Nruns)
	GLINT.secondMomentFunction[GLINT.geomIndex,*] /= GLINT.secondMomentFunction[GLINT.geomIndex,0L]
	sunglintPower /= DOUBLE(Nruns)

	glintBispectrum /= DOUBLE(Nruns)
	glintComponentPower /= DOUBLE(Nruns)
	glintSumPower /= DOUBLE(Nruns)
	FOR j=0L,NN4 DO BEGIN
		FOR i=j,NN2-j DO BEGIN
			IF (SQRT(glintComponentPower[i,j])*SQRT(glintSumPower[i,j]) GT 10.D^(-12.D)) THEN BEGIN
				glintBicoherence[i,j] = ABS(glintBispectrum[i,j])/(SQRT(glintComponentPower[i,j])*SQRT(glintSumPower[i,j]))
			ENDIF ELSE BEGIN
				glintBicoherence[i,j] = 0.D
			ENDELSE
		ENDFOR
	ENDFOR
	PRINT,'Glint first moment: ',GLINT.firstMoment[GLINT.geomIndex],FORMAT='(A,E16.7)'

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;; Compute the slope and glint triple cumulant functions     ;;;
	;;; from the slope and glint bispectra                        ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	;;; Apply the symmetry relationships to the slope and glint bispectra
	bispectrumSymmetry,slopeBispectrum,NN
	bispectrumSymmetry,glintBispectrum,NN

	;;; Compute the covariance functions
	thirdSlopeMomentFunction = DCOMPLEXARR(NN,NN)
	thirdGlintMomentFunction = DCOMPLEXARR(NN,NN)
	thirdSlopeMomentFunction = FFT(slopeBispectrum,/DOUBLE,/INVERSE)
	thirdGlintMomentFunction = FFT(glintBispectrum,/DOUBLE,/INVERSE)

	PRINT,''
	PRINT,"    Slope third moment from bicovariance: ",DOUBLE(thirdSlopeMomentFunction[0L,0L]),FORMAT='(A,E16.7)'
	PRINT,"      Slope third moment from bispectrum: ",TOTAL(DOUBLE(slopeBispectrum)),FORMAT='(A,E16.7)'
	PRINT,''
	PRINT,"    Glint third moment from bicovariance: ",DOUBLE(thirdGlintMomentFunction[0L,0L]),FORMAT='(A,E16.7)'
	PRINT,"      Glint third moment from bispectrum: ",TOTAL(DOUBLE(glintBispectrum)),FORMAT='(A,E16.7)'

	SLOPE.thirdMomentFunction[0L:NN/2L,0L:NN/2L] = $
		DOUBLE(thirdSlopeMomentFunction[0L:NN/2L,0L:NN/2L])/DOUBLE(thirdSlopeMomentFunction[0L,0L])

	GLINT.thirdMomentFunctionGeom[GLINT.geomIndex].thirdMomentFunction[0L:NN/2L,0L:NN/2L] = $
		DOUBLE(thirdGlintMomentFunction[0L:NN/2L,0L:NN/2L])/DOUBLE(thirdGlintMomentFunction[0L,0L])

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;                                                                     ;;;;;
;;;;;          Do the plotting for the bispectra and bicovariance         ;;;;;
;;;;;                                                                     ;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	;!P.FONT=-1
	;DEVICE, DECOMPOSED=0, RETAIN=2
	;LOADCT,39
	CONTOUR_LEVELS=20L
	N_plot = MIN([NN,N_ELEMENTS(WHERE(k LT k_PlotMax*k_N))])
;	N_plot = N_ELEMENTS(WHERE(k LT k_PlotMax*k_N))
	PRINT,"N_Plot = ",N_Plot

	rotate_z=-20.
	rotate_x=10.
	z_axis=2

	xwinsize=800
	ywinsize=450
	!P.MULTI=0

;	windowIndex++
;	WINDOW,windowIndex,xsize = xwinsize,ysize = ywinsize,title='Total Average Slope Power Spectrum',RETAIN=2
;	PLOT,k/k_N,primarySlopeAvgPower,xrange=[k_PlotMin,k_PlotMax],yrange=[0.D,1.2D*MAX(primarySlopeAvgPower)], $
;		xtitle='k',/xstyle,ytitle='Total Slope AvgPower',linestyle=0,/ystyle,charsize=chsize
;	OPLOT,k/k_N,totalSlopeAvgPower,PSYM=2
;;	OPLOT,k/k_N,primarySlopeAvgPower,linestyle=1
;	OPLOT,k/k_N,nlSlopeAvgPower,linestyle=2
;;	OPLOT,k[NLCOUPLING.free1]/k_N,totalSlopeAvgPower[NLCOUPLING.free1],PSYM=6
;
;	ysize=450
;	xsize=ROUND(1.2*ysize)
;; 	!P.MULTI=0
;	!P.MULTI = [0,2,1]
;
;	windowIndex++
;	window,windowIndex,xsize = 2*xsize,ysize = ysize,title='Slope real bispectrum',RETAIN=2
;	bispecMin = MIN([DOUBLE(slopeBispectrum),IMAGINARY(slopeBispectrum)])
;	bispecMax = MAX([DOUBLE(slopeBispectrum),IMAGINARY(slopeBispectrum)])
;	SURFACE,DOUBLE(slopeBispectrum[0L:N_plot-1L,0L:N_plot-1L]),kx[0L:N_plot-1L]/k_N,ky[0L:N_plot-1L]/k_N,AZ=rotate_z,$
;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax],$
;		zrange=[bispecMin,bispecMax], $
;		/xstyle,/ystyle,/zstyle,$
;		charsize=2.5,xtitle='k_x/k_N',ytitle='k_y/k_N',ZAXIS=z_axis
;	contourMin=MIN(DOUBLE(slopeBispectrum))
;	contourMax=MAX(DOUBLE(slopeBispectrum))
;	contourRange=contourMax-contourMin
;	CONTOUR,DOUBLE(slopeBispectrum[0L:N_plot-1L,0L:N_plot-1L]),kx[0L:N_plot-1L]/k_N,ky[0L:N_plot-1L]/k_N, $
;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax], $
;		xtitle='k_x/k_N',ytitle='k_y/k_N', $
;		/CELL_FILL,/ISOTROPIC,LEVELS=contourMin+contourRange*FINDGEN(CONTOUR_LEVELS)/DOUBLE(CONTOUR_LEVELS-1L)
	PRINT,'Maximum Value of Slope Real Bispectrum: ',MAX([DOUBLE(slopeBispectrum)])
	PRINT,'Minimum Value of Slope Real Bispectrum: ',MIN([DOUBLE(slopeBispectrum)])
;
;	ysize=450
;	xsize=ROUND(1.2*ysize)
;; 	!P.MULTI=0
;	!P.MULTI = [0,2,1]
;
;	windowIndex++
;	window,windowIndex,xsize = 2*xsize,ysize = ysize,title='Total Slope real bispectrum',RETAIN=2
;	SURFACE,DOUBLE(SHIFT(slopeBispectrum,NN/2L-1L,NN/2L-1L)),kx/k_N,ky/k_N,AZ=rotate_z,$
;;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax],$
;		zrange=[bispecMin,bispecMax], $
;		/xstyle,/ystyle,/zstyle,$
;		charsize=2.5,xtitle='k_x/k_N',ytitle='k_y/k_N',ZAXIS=z_axis
;	contourMin=MIN(DOUBLE(slopeBispectrum))
;	contourMax=MAX(DOUBLE(slopeBispectrum))
;	contourRange=contourMax-contourMin
;	CONTOUR,DOUBLE(SHIFT(slopeBispectrum,NN/2L-1L,NN/2L-1L)),kx/k_N,ky/k_N, $
;;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax], $
;		xtitle='k_x/k_N',ytitle='k_y/k_N',/xstyle,/ystyle,/zstyle, $
;		/CELL_FILL,/ISOTROPIC,LEVELS=contourMin+contourRange*FINDGEN(CONTOUR_LEVELS)/DOUBLE(CONTOUR_LEVELS-1L)
	PRINT,'Maximum Value of Total Slope Real Bispectrum: ',MAX([DOUBLE(slopeBispectrum)])
	PRINT,'Minimum Value of Total Slope Real Bispectrum: ',MIN([DOUBLE(slopeBispectrum)])
;
;	ysize=450
;	xsize=ROUND(1.2*ysize)
;; 	!P.MULTI=0
;	!P.MULTI = [0,2,1]
;
;	windowIndex++
;	window,windowIndex,xsize = 2*xsize,ysize = ysize,title='Slope imaginary bispectrum',RETAIN=2
;	SURFACE,IMAGINARY(slopeBispectrum[0L:N_plot-1L,0L:N_plot-1L]),kx[0L:N_plot-1L]/k_N,ky[0L:N_plot-1L]/k_N,AZ=rotate_z,$
;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax],$
;		zrange=[bispecMin,bispecMax], $
;		/xstyle,/ystyle,/zstyle,$
;		charsize=2.5,xtitle='k_x/k_N',ytitle='k_y/k_N',ZAXIS=z_axis
;	contourMin=MIN(IMAGINARY(slopeBispectrum))
;	contourMax=MAX(IMAGINARY(slopeBispectrum))
;	contourRange=contourMax-contourMin
;	CONTOUR,IMAGINARY(slopeBispectrum[0L:N_plot-1L,0L:N_plot-1L]),kx[0L:N_plot-1L]/k_N,ky[0L:N_plot-1L]/k_N, $
;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax], $
;		xtitle='k_x/k_N',ytitle='k_y/k_N', $
;		/CELL_FILL,/ISOTROPIC,LEVELS=contourMin+contourRange*FINDGEN(CONTOUR_LEVELS)/DOUBLE(CONTOUR_LEVELS-1L)
	PRINT,'Maximum Value of Slope Imaginary Bispectrum: ',MAX([IMAGINARY(slopeBispectrum)])
	PRINT,'Minimum Value of Slope Imaginary Bispectrum: ',MIN([IMAGINARY(slopeBispectrum)])
;
;	ysize=450
;	xsize=ROUND(1.2*ysize)
;; 	!P.MULTI=0
;	!P.MULTI = [0,2,1]
;
;	windowIndex++
;	window,windowIndex,xsize = 2*xsize,ysize = ysize,title='Total Slope imaginary bispectrum',RETAIN=2
;	SURFACE,IMAGINARY(SHIFT(slopeBispectrum,NN/2L-1L,NN/2L-1L)),kx/k_N,ky/k_N,AZ=rotate_z,$
;;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax],$
;		zrange=[bispecMin,bispecMax], $
;		/xstyle,/ystyle,/zstyle,$
;		charsize=2.5,xtitle='k_x/k_N',ytitle='k_y/k_N',ZAXIS=z_axis
;	contourMin=MIN(IMAGINARY(slopeBispectrum))
;	contourMax=MAX(IMAGINARY(slopeBispectrum))
;	contourRange=contourMax-contourMin
;	CONTOUR,IMAGINARY(SHIFT(slopeBispectrum,NN/2L-1L,NN/2L-1L)),kx/k_N,ky/k_N, $
;;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax], $
;		xtitle='k_x/k_N',ytitle='k_y/k_N',/xstyle,/ystyle,/zstyle, $
;		/CELL_FILL,/ISOTROPIC,LEVELS=contourMin+contourRange*FINDGEN(CONTOUR_LEVELS)/DOUBLE(CONTOUR_LEVELS-1L)
	PRINT,'Maximum Value of Total Slope Imaginary Bispectrum: ',MAX([IMAGINARY(slopeBispectrum)])
	PRINT,'Minimum Value of Total Slope Imaginary Bispectrum: ',MIN([IMAGINARY(slopeBispectrum)])
;
;	ysize=450
;	xsize=ROUND(1.2*ysize)
;; 	!P.MULTI=0
;	!P.MULTI = [0,2,1]
;
;	windowIndex++
;	window,windowIndex,xsize = 2*xsize,ysize = ysize,title='Total Uncentered Real Third Slope Moment Function',RETAIN=2
;	bispecMin = MIN([DOUBLE(thirdSlopeMomentFunction),IMAGINARY(thirdSlopeMomentFunction)])
;	bispecMax = MAX([DOUBLE(thirdSlopeMomentFunction),IMAGINARY(thirdSlopeMomentFunction)])
;	SURFACE,DOUBLE(thirdSlopeMomentFunction),AZ=rotate_z,$
;;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax],$
;		zrange=[5.*bispecMin,5.*bispecMax], $
;		/xstyle,/ystyle,/zstyle,$
;		charsize=2.5,xtitle='k_x/k_N',ytitle='k_y/k_N',ZAXIS=z_axis
;	contourMin=MIN(DOUBLE(thirdSlopeMomentFunction))
;	contourMax=MAX(DOUBLE(thirdSlopeMomentFunction))
;	contourRange=contourMax-contourMin
;	CONTOUR,DOUBLE(thirdSlopeMomentFunction), $
;;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax], $
;		xtitle='x',ytitle='y',/xstyle,/ystyle,/zstyle, $
;		/CELL_FILL,/ISOTROPIC,LEVELS=contourMin+contourRange*FINDGEN(CONTOUR_LEVELS)/DOUBLE(CONTOUR_LEVELS-1L)
	PRINT,'Maximum Value of Third Slope Moment Function (Real): ',MAX([DOUBLE(thirdSlopeMomentFunction)])
	PRINT,'Minimum Value of Third Slope Moment Function (Real): ',MIN([DOUBLE(thirdSlopeMomentFunction)])
;
;	ysize=450
;	xsize=ROUND(1.2*ysize)
;; 	!P.MULTI=0
;	!P.MULTI = [0,2,1]
;
;	windowIndex++
;	window,windowIndex,xsize = 2*xsize,ysize = ysize,title='Total Uncenetered Imaginary Third Slope Moment Function',RETAIN=2
;	SURFACE,IMAGINARY(thirdSlopeMomentFunction),tau_x,tau_y,AZ=rotate_z,$
;;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax],$
;		zrange=[bispecMin,bispecMax], $
;		/xstyle,/ystyle,/zstyle,$
;		charsize=2.5,xtitle='k_x/k_N',ytitle='k_y/k_N',ZAXIS=z_axis
;	contourMin=MIN(IMAGINARY(thirdSlopeMomentFunction))
;	contourMax=MAX(IMAGINARY(thirdSlopeMomentFunction))
;	contourRange=contourMax-contourMin
;	CONTOUR,IMAGINARY(thirdSlopeMomentFunction),tau_x,tau_y, $
;;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax], $
;		xtitle='x',ytitle='y',/xstyle,/ystyle,/zstyle, $
;		/CELL_FILL,/ISOTROPIC,LEVELS=contourMin+contourRange*FINDGEN(CONTOUR_LEVELS)/DOUBLE(CONTOUR_LEVELS-1L)
	PRINT,'Maximum Value of Third Slope Moment Function (Imag): ',MAX([IMAGINARY(thirdSlopeMomentFunction)])
	PRINT,'Minimum Value of Third Slope Moment Function (Imag): ',MIN([IMAGINARY(thirdSlopeMomentFunction)])

;	ysize=450
;	xsize=ROUND(1.2*ysize)
;; 	!P.MULTI=0
;	!P.MULTI = [0,2,1]
;	windowIndex++
;	window,windowIndex,xsize = 2*xsize,ysize = ysize,title='Slope Component Power Spectrum',RETAIN=2
;	SURFACE,slopeComponentPower[0L:N_plot-1L,0L:N_plot-1L],kx[0L:N_plot-1L]/k_N,ky[0L:N_plot-1L]/k_N,AZ=rotate_z,$
;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax], $
;		/xstyle,/ystyle,/zstyle,$
;		charsize=2.5,xtitle='k_x/k_N',ytitle='k_y/k_N',ZAXIS=z_axis
;	contourMin=MIN(DOUBLE(slopeComponentPower))
;	contourMax=MAX(DOUBLE(slopeComponentPower))
;	contourRange=contourMax-contourMin
;	CONTOUR,DOUBLE(slopeComponentPower[0L:N_plot-1L,0L:N_plot-1L]),kx[0L:N_plot-1L]/k_N,ky[0L:N_plot-1L]/k_N, $
;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax], $
;		xtitle='k_x/k_N',ytitle='k_y/k_N', $
;		/CELL_FILL,/ISOTROPIC,LEVELS=contourMin+contourRange*FINDGEN(CONTOUR_LEVELS)/DOUBLE(CONTOUR_LEVELS-1L)
;	PRINT,'Maximum Value of Slope Component Power Spectrum: ',contourMax
;	PRINT,'Minimum Value of Slope Component Power Spectrum: ',contourMin
;
;
;	ysize=450
;	xsize=ROUND(1.2*ysize)
;; 	!P.MULTI=0
;	!P.MULTI = [0,2,1]
;	windowIndex++
;	window,windowIndex,xsize = 2*xsize,ysize = ysize,title='Total Slope Component Power Spectrum',RETAIN=2
;	SURFACE,slopeComponentPower,kx/k_N,ky/k_N,AZ=rotate_z,$
;;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax], $
;		/xstyle,/ystyle,/zstyle,$
;		charsize=2.5,xtitle='k_x/k_N',ytitle='k_y/k_N',ZAXIS=z_axis
;	contourMin=MIN(DOUBLE(slopeComponentPower))
;	contourMax=MAX(DOUBLE(slopeComponentPower))
;	contourRange=contourMax-contourMin
;	CONTOUR,DOUBLE(slopeComponentPower),kx/k_N,ky/k_N, $
;;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax], $
;		xtitle='k_x/k_N',ytitle='k_y/k_N',/xstyle,/ystyle,/zstyle, $
;		/CELL_FILL,/ISOTROPIC,LEVELS=contourMin+contourRange*FINDGEN(CONTOUR_LEVELS)/DOUBLE(CONTOUR_LEVELS-1L)
;	PRINT,'Maximum Value of Total Slope Component Power Spectrum: ',contourMax
;	PRINT,'Minimum Value of Total Slope Component Power Spectrum: ',contourMin
;
;	ysize=450
;	xsize=ROUND(1.2*ysize)
;; 	!P.MULTI=0
;	!P.MULTI = [0,2,1]
;	windowIndex++
;	window,windowIndex,xsize = 2*xsize,ysize = ysize,title='Slope Sum Power Spectrum',RETAIN=2
;	SURFACE,slopeSumPower[0L:N_plot-1L,0L:N_plot-1L],kx[0L:N_plot-1L]/k_N,ky[0L:N_plot-1L]/k_N,AZ=rotate_z,$
;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax], $
;		/xstyle,/ystyle,/zstyle,$
;		charsize=2.5,xtitle='k_x/k_N',ytitle='k_y/k_N',ZAXIS=z_axis
;	contourMin=MIN(DOUBLE(slopeSumPower))
;	contourMax=MAX(DOUBLE(slopeSumPower))
;	contourRange=contourMax-contourMin
;	CONTOUR,DOUBLE(slopeSumPower[0L:N_plot-1L,0L:N_plot-1L]),kx[0L:N_plot-1L]/k_N,ky[0L:N_plot-1L]/k_N, $
;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax], $
;		xtitle='k_x/k_N',ytitle='k_y/k_N', $
;		/CELL_FILL,/ISOTROPIC,LEVELS=contourMin+contourRange*FINDGEN(CONTOUR_LEVELS)/DOUBLE(CONTOUR_LEVELS-1L)
;	PRINT,'Maximum Value of Slope Sum Power Spectrum: ',contourMax
;	PRINT,'Minimum Value of Slope Sum Power Spectrum: ',contourMin
;
;
;	ysize=450
;	xsize=ROUND(1.2*ysize)
;; 	!P.MULTI=0
;	!P.MULTI = [0,2,1]
;	windowIndex++
;	window,windowIndex,xsize = 2*xsize,ysize = ysize,title='Total Slope Sum Power Spectrum',RETAIN=2
;	SURFACE,slopeSumPower,kx/k_N,ky/k_N,AZ=rotate_z,$
;;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax], $
;		/xstyle,/ystyle,/zstyle,$
;		charsize=2.5,xtitle='k_x/k_N',ytitle='k_y/k_N',ZAXIS=z_axis
;	contourMin=MIN(DOUBLE(slopeSumPower))
;	contourMax=MAX(DOUBLE(slopeSumPower))
;	contourRange=contourMax-contourMin
;	CONTOUR,DOUBLE(slopeSumPower),kx/k_N,ky/k_N, $
;;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax], $
;		xtitle='k_x/k_N',ytitle='k_y/k_N',/xstyle,/ystyle,/zstyle, $
;		/CELL_FILL,/ISOTROPIC,LEVELS=contourMin+contourRange*FINDGEN(CONTOUR_LEVELS)/DOUBLE(CONTOUR_LEVELS-1L)
;	PRINT,'Maximum Value of Total Slope Sum Power Spectrum: ',contourMax
;	PRINT,'Minimum Value of Total Slope Sum Power Spectrum: ',contourMin
;
;	ysize=450
;	xsize=ROUND(1.2*ysize)
;; 	!P.MULTI=0
;	!P.MULTI = [0,2,1]
;	windowIndex++
;	window,windowIndex,xsize = 2*xsize,ysize = ysize,title='Slope bicoherence',RETAIN=2
;	SURFACE,slopeBicoherence[0L:N_plot-1L,0L:N_plot-1L]^2.,kx[0L:N_plot-1L]/k_N,ky[0L:N_plot-1L]/k_N,AZ=rotate_z,$
;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax], $
;		zrange=[0.,1.2D],/xstyle,/ystyle,/zstyle,$
;		charsize=2.5,xtitle='k_x/k_N',ytitle='k_y/k_N',ZAXIS=z_axis
;	contourMin=MIN(DOUBLE(slopeBicoherence^2.D))
;	contourMax=MAX(DOUBLE(slopeBicoherence^2.D))
;	contourRange=contourMax-contourMin
;	CONTOUR,DOUBLE(slopeBicoherence[0L:N_plot-1L,0L:N_plot-1L]^2.D),kx[0L:N_plot-1L]/k_N,ky[0L:N_plot-1L]/k_N, $
;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax], $
;		xtitle='k_x/k_N',ytitle='k_y/k_N', $
;		/CELL_FILL,/ISOTROPIC,LEVELS=contourMin+contourRange*FINDGEN(CONTOUR_LEVELS)/DOUBLE(CONTOUR_LEVELS-1L)
	PRINT,'Maximum Value of Slope Bicoherence: ',MIN(DOUBLE(slopeBicoherence^2.D))
	PRINT,'Minimum Value of Slope Bicoherence: ',MAX(DOUBLE(slopeBicoherence^2.D))
;
;
;	ysize=450
;	xsize=ROUND(1.2*ysize)
;; 	!P.MULTI=0
;	!P.MULTI = [0,2,1]
;	windowIndex++
;	window,windowIndex,xsize = 2*xsize,ysize = ysize,title='Total Slope bicoherence',RETAIN=2
;	SURFACE,slopeBicoherence^2.,kx/k_N,ky/k_N,AZ=rotate_z,$
;;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax], $
;		zrange=[0.,1.2D],/xstyle,/ystyle,/zstyle,$
;		charsize=2.5,xtitle='k_x/k_N',ytitle='k_y/k_N',ZAXIS=z_axis
;	contourMin=MIN(DOUBLE(slopeBicoherence^2.D))
;	contourMax=MAX(DOUBLE(slopeBicoherence^2.D))
;	contourRange=contourMax-contourMin
;	CONTOUR,DOUBLE(slopeBicoherence^2.D),kx/k_N,ky/k_N, $
;;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax], $
;		xtitle='k_x/k_N',ytitle='k_y/k_N',/xstyle,/ystyle,/zstyle, $
;		/CELL_FILL,/ISOTROPIC,LEVELS=contourMin+contourRange*FINDGEN(CONTOUR_LEVELS)/DOUBLE(CONTOUR_LEVELS-1L)
;	PRINT,'Maximum Value of Slope Bicoherence: ',MIN(DOUBLE(slopeBicoherence^2.D))
;	PRINT,'Minimum Value of Slope Bicoherence: ',MAX(DOUBLE(slopeBicoherence^2.D))
;
;	xwinsize=800
;	ywinsize=450
;	!P.MULTI=0
;
;	windowIndex++
;	WINDOW,windowIndex,xsize = xwinsize,ysize = ywinsize,title='Sunglint Power Spectrum',RETAIN=2
;	PLOT,k/k_N,sunglintPower,xrange=[0.D,2.0D],yrange=[0.D,1.2D*MAX(sunglintPower)], $
;		xtitle='k/k_N',/xstyle,ytitle='Total Sunglint Power',linestyle=0,/ystyle,charsize=chsize
;	OPLOT,k/k_N,ABS(MAX(sunglintPower))*EXP(-(k/(0.8*k_N))^2.D),linestyle=2
;
;	xwinsize=800
;	ywinsize=450
;	!P.MULTI=0
;
;	windowIndex++
;	WINDOW,windowIndex,xsize = xwinsize,ysize = ywinsize,title='Sunglint autocorrelation',RETAIN=2
;	PLOT,tau,FFT(sunglintPower,/INVERSE,/DOUBLE), xrange=[0.D,1.0D], $
;		yrange=[-0.0001,0.0001], $
;		xtitle='tau',/xstyle,ytitle='Glint autocorrelation',PSYM=2,/ystyle,charsize=chsize
;;	OPLOT,k/k_N,sunglintPower,PSYM=2
;;	PRINT,(DOUBLE(FFT(sunglintPower,/INVERSE,/DOUBLE)))[0L:50L]
;
;	ysize=450
;	xsize=ROUND(1.2*ysize)
;; 	!P.MULTI=0
;	!P.MULTI = [0,2,1]
;
;	windowIndex++
;	window,windowIndex,xsize = 2*xsize,ysize = ysize,title='glint real bispectrum',RETAIN=2
;	bispecMin = MIN([DOUBLE(glintBispectrum),IMAGINARY(glintBispectrum)])
;	bispecMax = MAX([DOUBLE(glintBispectrum),IMAGINARY(glintBispectrum)])
;	SURFACE,DOUBLE(glintBispectrum[0L:N_plot-1L,0L:N_plot-1L]),kx[0L:N_plot-1L]/k_N,ky[0L:N_plot-1L]/k_N,AZ=rotate_z,$
;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax],$
;		zrange=[bispecMin,bispecMax], $
;		/xstyle,/ystyle,/zstyle,$
;		charsize=2.5,xtitle='k_x/k_N',ytitle='k_y/k_N',ZAXIS=z_axis
;	contourMin=MIN(DOUBLE(glintBispectrum))
;	contourMax=MAX(DOUBLE(glintBispectrum))
;	contourRange=contourMax-contourMin
;	CONTOUR,DOUBLE(glintBispectrum[0L:N_plot-1L,0L:N_plot-1L]),kx[0L:N_plot-1L]/k_N,ky[0L:N_plot-1L]/k_N, $
;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax], $
;		xtitle='k_x/k_N',ytitle='k_y/k_N', $
;		/CELL_FILL,/ISOTROPIC,LEVELS=contourMin+contourRange*FINDGEN(CONTOUR_LEVELS)/DOUBLE(CONTOUR_LEVELS-1L)
	PRINT,'Maximum Value of glint Real Bispectrum: ',MAX([DOUBLE(glintBispectrum)])
	PRINT,'Minimum Value of glint Real Bispectrum: ',MIN([DOUBLE(glintBispectrum)])
;
;	ysize=450
;	xsize=ROUND(1.2*ysize)
;; 	!P.MULTI=0
;	!P.MULTI = [0,2,1]
;
;	windowIndex++
;	window,windowIndex,xsize = 2*xsize,ysize = ysize,title='Total glint real bispectrum',RETAIN=2
;	SURFACE,DOUBLE(SHIFT(glintBispectrum,NN/2L-1L,NN/2L-1L)),kx/k_N,ky/k_N,AZ=rotate_z,$
;;	SURFACE,DOUBLE(glintBispectrum),kx/k_N,ky/k_N,AZ=rotate_z,$
;;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax],$
;		zrange=[bispecMin,bispecMax], $
;		/xstyle,/ystyle,/zstyle,$
;		charsize=2.5,xtitle='k_x/k_N',ytitle='k_y/k_N',ZAXIS=z_axis
;	contourMin=MIN(DOUBLE(glintBispectrum))
;	contourMax=MAX(DOUBLE(glintBispectrum))
;	contourRange=contourMax-contourMin
;	CONTOUR,DOUBLE(SHIFT(glintBispectrum,NN/2L-1L,NN/2L-1L)),kx/k_N,ky/k_N, $
;;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax], $
;		xtitle='k_x/k_N',ytitle='k_y/k_N',/xstyle,/ystyle,/zstyle, $
;		/CELL_FILL,/ISOTROPIC,LEVELS=contourMin+contourRange*FINDGEN(CONTOUR_LEVELS)/DOUBLE(CONTOUR_LEVELS-1L)
;	PRINT,'Maximum Value of Total glint Real Bispectrum: ',MAX([DOUBLE(glintBispectrum)])
;	PRINT,'Minimum Value of Total glint Real Bispectrum: ',MIN([DOUBLE(glintBispectrum)])
;
;	ysize=450
;	xsize=ROUND(1.2*ysize)
;; 	!P.MULTI=0
;	!P.MULTI = [0,2,1]
;
;	windowIndex++
;	window,windowIndex,xsize = 2*xsize,ysize = ysize,title='glint imaginary bispectrum',RETAIN=2
;	SURFACE,IMAGINARY(glintBispectrum[0L:N_plot-1L,0L:N_plot-1L]),kx[0L:N_plot-1L]/k_N,ky[0L:N_plot-1L]/k_N,AZ=rotate_z,$
;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax],$
;		zrange=[bispecMin,bispecMax], $
;		/xstyle,/ystyle,/zstyle,$
;		charsize=2.5,xtitle='k_x/k_N',ytitle='k_y/k_N',ZAXIS=z_axis
;	contourMin=MIN(IMAGINARY(glintBispectrum))
;	contourMax=MAX(IMAGINARY(glintBispectrum))
;	contourRange=contourMax-contourMin
;	CONTOUR,IMAGINARY(glintBispectrum[0L:N_plot-1L,0L:N_plot-1L]),kx[0L:N_plot-1L]/k_N,ky[0L:N_plot-1L]/k_N, $
;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax], $
;		xtitle='k_x/k_N',ytitle='k_y/k_N', $
;		/CELL_FILL,/ISOTROPIC,LEVELS=contourMin+contourRange*FINDGEN(CONTOUR_LEVELS)/DOUBLE(CONTOUR_LEVELS-1L)
	PRINT,'Maximum Value of glint Imaginary Bispectrum: ',MAX([IMAGINARY(glintBispectrum)])
	PRINT,'Minimum Value of glint Imaginary Bispectrum: ',MIN([IMAGINARY(glintBispectrum)])
;
;	windowIndex++
;	window,windowIndex,xsize = 2*xsize,ysize = ysize,title='Total glint imaginary bispectrum',RETAIN=2
;	SURFACE,IMAGINARY(SHIFT(glintBispectrum,NN/2L-1L,NN/2L-1L)),kx/k_N,ky/k_N,AZ=rotate_z,$
;;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax],$
;		zrange=[bispecMin,bispecMax], $
;		/xstyle,/ystyle,/zstyle,$
;		charsize=2.5,xtitle='k_x/k_N',ytitle='k_y/k_N',ZAXIS=z_axis
;	contourMin=MIN(IMAGINARY(glintBispectrum))
;	contourMax=MAX(IMAGINARY(glintBispectrum))
;	contourRange=contourMax-contourMin
;	CONTOUR,IMAGINARY(SHIFT(glintBispectrum,NN/2L-1L,NN/2L-1L)),kx/k_N,ky/k_N, $
;;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax], $
;		xtitle='k_x/k_N',ytitle='k_y/k_N',/xstyle,/ystyle,/zstyle, $
;		/CELL_FILL,/ISOTROPIC,LEVELS=contourMin+contourRange*FINDGEN(CONTOUR_LEVELS)/DOUBLE(CONTOUR_LEVELS-1L)
;	PRINT,'Maximum Value of Total glint Imaginary Bispectrum: ',MAX([IMAGINARY(glintBispectrum)])
;	PRINT,'Minimum Value of Total glint Imaginary Bispectrum: ',MIN([IMAGINARY(glintBispectrum)])
;
;	ysize=450
;	xsize=ROUND(1.2*ysize)
;; 	!P.MULTI=0
;	!P.MULTI = [0,2,1]
;
;	windowIndex++
;	window,windowIndex,xsize = 2*xsize,ysize = ysize,title='Total Uncentered Real Third glint Moment Function',RETAIN=2
;	bispecMin = MIN([DOUBLE(thirdGlintMomentFunction),IMAGINARY(thirdGlintMomentFunction)])
;	bispecMax = MAX([DOUBLE(thirdGlintMomentFunction),IMAGINARY(thirdGlintMomentFunction)])
;	SURFACE,DOUBLE(thirdGlintMomentFunction),tau_x,tau_y,AZ=rotate_z,$
;;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax],$
;		zrange=[5.*bispecMin,5.*bispecMax], $
;		/xstyle,/ystyle,/zstyle,$
;		charsize=2.5,xtitle='k_x/k_N',ytitle='k_y/k_N',ZAXIS=z_axis
;	contourMin=MIN(DOUBLE(thirdGlintMomentFunction))
;	contourMax=MAX(DOUBLE(thirdGlintMomentFunction))
;	contourRange=contourMax-contourMin
;	CONTOUR,DOUBLE(thirdGlintMomentFunction),tau_x,tau_y, $
;;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax], $
;		xtitle='x',ytitle='y',/xstyle,/ystyle,/zstyle, $
;		/CELL_FILL,/ISOTROPIC,LEVELS=contourMin+contourRange*FINDGEN(CONTOUR_LEVELS)/DOUBLE(CONTOUR_LEVELS-1L)
	PRINT,'Maximum Value of Real Third glint Moment Function: ',MAX([DOUBLE(thirdGlintMomentFunction)])
	PRINT,'Minimum Value of Real Third glint Moment Function: ',MIN([DOUBLE(thirdGlintMomentFunction)])
;
;	ysize=450
;	xsize=ROUND(1.2*ysize)
;; 	!P.MULTI=0
;	!P.MULTI = [0,2,1]
;
;	windowIndex++
;	window,windowIndex,xsize = 2*xsize,ysize = ysize,title='Total Uncentered Imaginary Third glint Moment Function',RETAIN=2
;	bispecMin = MIN([DOUBLE(thirdGlintMomentFunction),IMAGINARY(thirdGlintMomentFunction)])
;	bispecMax = MAX([DOUBLE(thirdGlintMomentFunction),IMAGINARY(thirdGlintMomentFunction)])
;	SURFACE,IMAGINARY(thirdGlintMomentFunction),tau_x,tau_y,AZ=rotate_z,$
;;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax],$
;		zrange=[5.*bispecMin,5.*bispecMax], $
;		/xstyle,/ystyle,/zstyle,$
;		charsize=2.5,xtitle='k_x/k_N',ytitle='k_y/k_N',ZAXIS=z_axis
;	contourMin=MIN(IMAGINARY(thirdGlintMomentFunction))
;	contourMax=MAX(IMAGINARY(thirdGlintMomentFunction))
;	contourRange=contourMax-contourMin
;	CONTOUR,IMAGINARY(thirdGlintMomentFunction),tau_x,tau_y, $
;;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax], $
;		xtitle='x',ytitle='y',/xstyle,/ystyle,/zstyle, $
;		/CELL_FILL,/ISOTROPIC,LEVELS=contourMin+contourRange*FINDGEN(CONTOUR_LEVELS)/DOUBLE(CONTOUR_LEVELS-1L)
	PRINT,'Maximum Value of Imaginary Third glint Moment Function: ',MAX([IMAGINARY(thirdGlintMomentFunction)])
	PRINT,'Minimum Value of Imaginary Third glint Moment Function: ',MIN([IMAGINARY(thirdGlintMomentFunction)])
;
;	ysize=450
;	xsize=ROUND(1.2*ysize)
;; 	!P.MULTI=0
;	!P.MULTI = [0,2,1]
;	windowIndex++
;	window,windowIndex,xsize = 2*xsize,ysize = ysize,title='glint Component Power Spectrum',RETAIN=2
;	SURFACE,glintComponentPower[0L:N_plot-1L,0L:N_plot-1L],kx[0L:N_plot-1L]/k_N,ky[0L:N_plot-1L]/k_N,AZ=rotate_z,$
;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax], $
;		/xstyle,/ystyle,/zstyle,$
;		charsize=2.5,xtitle='k_x/k_N',ytitle='k_y/k_N',ZAXIS=z_axis
;	contourMin=MIN(DOUBLE(glintComponentPower))
;	contourMax=MAX(DOUBLE(glintComponentPower))
;	contourRange=contourMax-contourMin
;	CONTOUR,DOUBLE(glintComponentPower[0L:N_plot-1L,0L:N_plot-1L]),kx[0L:N_plot-1L]/k_N,ky[0L:N_plot-1L]/k_N, $
;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax], $
;		xtitle='k_x/k_N',ytitle='k_y/k_N', $
;		/CELL_FILL,/ISOTROPIC,LEVELS=contourMin+contourRange*FINDGEN(CONTOUR_LEVELS)/DOUBLE(CONTOUR_LEVELS-1L)
;	PRINT,'Maximum Value of glint Component Power Spectrum: ',contourMax
;	PRINT,'Minimum Value of glint Component Power Spectrum: ',contourMin
;
;
;	ysize=450
;	xsize=ROUND(1.2*ysize)
;; 	!P.MULTI=0
;	!P.MULTI = [0,2,1]
;	windowIndex++
;	window,windowIndex,xsize = 2*xsize,ysize = ysize,title='Total glint Component Power Spectrum',RETAIN=2
;	SURFACE,glintComponentPower,kx/k_N,ky/k_N,AZ=rotate_z,$
;;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax], $
;		/xstyle,/ystyle,/zstyle,$
;		charsize=2.5,xtitle='k_x/k_N',ytitle='k_y/k_N',ZAXIS=z_axis
;	contourMin=MIN(DOUBLE(glintComponentPower))
;	contourMax=MAX(DOUBLE(glintComponentPower))
;	contourRange=contourMax-contourMin
;	CONTOUR,DOUBLE(glintComponentPower),kx/k_N,ky/k_N, $
;;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax], $
;		xtitle='k_x/k_N',ytitle='k_y/k_N',/xstyle,/ystyle,/zstyle, $
;		/CELL_FILL,/ISOTROPIC,LEVELS=contourMin+contourRange*FINDGEN(CONTOUR_LEVELS)/DOUBLE(CONTOUR_LEVELS-1L)
;	PRINT,'Maximum Value of Total glint Component Power Spectrum: ',contourMax
;	PRINT,'Minimum Value of Total glint Component Power Spectrum: ',contourMin
;
;	ysize=450
;	xsize=ROUND(1.2*ysize)
;; 	!P.MULTI=0
;	!P.MULTI = [0,2,1]
;	windowIndex++
;	window,windowIndex,xsize = 2*xsize,ysize = ysize,title='glint Sum Power Spectrum',RETAIN=2
;	SURFACE,glintSumPower[0L:N_plot-1L,0L:N_plot-1L],kx[0L:N_plot-1L]/k_N,ky[0L:N_plot-1L]/k_N,AZ=rotate_z,$
;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax], $
;		/xstyle,/ystyle,/zstyle,$
;		charsize=2.5,xtitle='k_x/k_N',ytitle='k_y/k_N',ZAXIS=z_axis
;	contourMin=MIN(DOUBLE(glintSumPower))
;	contourMax=MAX(DOUBLE(glintSumPower))
;	contourRange=contourMax-contourMin
;	CONTOUR,DOUBLE(glintSumPower[0L:N_plot-1L,0L:N_plot-1L]),kx[0L:N_plot-1L]/k_N,ky[0L:N_plot-1L]/k_N, $
;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax], $
;		xtitle='k_x/k_N',ytitle='k_y/k_N', $
;		/CELL_FILL,/ISOTROPIC,LEVELS=contourMin+contourRange*FINDGEN(CONTOUR_LEVELS)/DOUBLE(CONTOUR_LEVELS-1L)
;	PRINT,'Maximum Value of glint Sum Power Spectrum: ',contourMax
;	PRINT,'Minimum Value of glint Sum Power Spectrum: ',contourMin
;
;
;	ysize=450
;	xsize=ROUND(1.2*ysize)
;; 	!P.MULTI=0
;	!P.MULTI = [0,2,1]
;	windowIndex++
;	window,windowIndex,xsize = 2*xsize,ysize = ysize,title='Total glint Sum Power Spectrum',RETAIN=2
;	SURFACE,glintSumPower,kx/k_N,ky/k_N,AZ=rotate_z,$
;;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax], $
;		/xstyle,/ystyle,/zstyle,$
;		charsize=2.5,xtitle='k_x/k_N',ytitle='k_y/k_N',ZAXIS=z_axis
;	contourMin=MIN(DOUBLE(glintSumPower))
;	contourMax=MAX(DOUBLE(glintSumPower))
;	contourRange=contourMax-contourMin
;	CONTOUR,DOUBLE(glintSumPower),kx/k_N,ky/k_N, $
;;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax], $
;		xtitle='k_x/k_N',ytitle='k_y/k_N',/xstyle,/ystyle,/zstyle, $
;		/CELL_FILL,/ISOTROPIC,LEVELS=contourMin+contourRange*FINDGEN(CONTOUR_LEVELS)/DOUBLE(CONTOUR_LEVELS-1L)
;	PRINT,'Maximum Value of Total glint Sum Power Spectrum: ',contourMax
;	PRINT,'Minimum Value of Total glint Sum Power Spectrum: ',contourMin
;
;	ysize=450
;	xsize=ROUND(1.2*ysize)
;; 	!P.MULTI=0
;	!P.MULTI = [0,2,1]
;	windowIndex++
;	window,windowIndex,xsize = 2*xsize,ysize = ysize,title='glint bicoherence',RETAIN=2
;	SURFACE,glintBicoherence[0L:N_plot-1L,0L:N_plot-1L]^2.,kx[0L:N_plot-1L]/k_N,ky[0L:N_plot-1L]/k_N,AZ=rotate_z,$
;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax], $
;;		zrange=[0.,1.2D], $
;		/xstyle,/ystyle,/zstyle,$
;		charsize=2.5,xtitle='k_x/k_N',ytitle='k_y/k_N',ZAXIS=z_axis
;	contourMin=MIN(DOUBLE(glintBicoherence^2.D))
;	contourMax=MAX(DOUBLE(glintBicoherence^2.D))
;	contourRange=contourMax-contourMin
;	CONTOUR,DOUBLE(glintBicoherence[0L:N_plot-1L,0L:N_plot-1L]^2.D),kx[0L:N_plot-1L]/k_N,ky[0L:N_plot-1L]/k_N, $
;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax], $
;		xtitle='k_x/k_N',ytitle='k_y/k_N', $
;		/CELL_FILL,/ISOTROPIC,LEVELS=contourMin+contourRange*FINDGEN(CONTOUR_LEVELS)/DOUBLE(CONTOUR_LEVELS-1L)
	PRINT,'Maximum Value of glint Bicoherence: ',MIN(DOUBLE(glintBicoherence^2.D))
	PRINT,'Minimum Value of glint Bicoherence: ',MAX(DOUBLE(glintBicoherence^2.D))
;
;
;	ysize=450
;	xsize=ROUND(1.2*ysize)
;; 	!P.MULTI=0
;	!P.MULTI = [0,2,1]
;	windowIndex++
;	window,windowIndex,xsize = 2*xsize,ysize = ysize,title='Total glint bicoherence',RETAIN=2
;	SURFACE,glintBicoherence^2.,kx/k_N,ky/k_N,AZ=rotate_z,$
;;		xrange=[k_PlotMin,k_PlotMax], $
;;		yrange=[k_PlotMin,k_PlotMax], $
;;		zrange=[0.,1.2D], $
;		/xstyle,/ystyle,/zstyle,$
;		charsize=2.5,xtitle='k_x/k_N',ytitle='k_y/k_N',ZAXIS=z_axis
;	contourMin=MIN(DOUBLE(glintBicoherence^2.D))
;	contourMax=MAX(DOUBLE(glintBicoherence^2.D))
;	contourRange=contourMax-contourMin
;	CONTOUR,DOUBLE(glintBicoherence^2.D),kx/k_N,ky/k_N, $
;;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax], $
;		xtitle='k_x/k_N',ytitle='k_y/k_N',/xstyle,/ystyle,/zstyle, $
;		/CELL_FILL,/ISOTROPIC,LEVELS=contourMin+contourRange*FINDGEN(CONTOUR_LEVELS)/DOUBLE(CONTOUR_LEVELS-1L)
;	PRINT,'Maximum Value of glint Bicoherence: ',MIN(DOUBLE(glintBicoherence^2.D))
;	PRINT,'Minimum Value of glint Bicoherence: ',MAX(DOUBLE(glintBicoherence^2.D))
END
