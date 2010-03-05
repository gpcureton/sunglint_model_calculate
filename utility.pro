;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Create tau abcissa with multile scales   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO createTauAbcissa, N,N_tau_samples,delta_x,tauLo,tauMid,tauHi,tau,tau_index,tauAbcissa
	
	PRINT,""
	PRINT,"<<< Input >>>"
	PRINT,"N = ",N
	PRINT,"N_tau_samples = ",N_tau_samples
	PRINT,"delta_x = ",delta_x
	PRINT,"tauLo = ",tauLo
	PRINT,"tauMid = ",tauMid
	PRINT,"tauHi = ",tauHi
	PRINT,""

	;;; Initialise the first half of the abcissa with an exponentially varying lag increment
	tauMid = MIN([DOUBLE(N/2L)*delta_x,tauMid])
	PRINT,"New tauMid = ",tauMid

	tau_index[1L:N_tau_samples/2L] = FIX((tauLo^(1L-FINDGEN(N_tau_samples/2L)/(N_tau_samples/2L-1L)) $
		*tauMid^(FINDGEN(N_tau_samples/2L)/(N_tau_samples/2L-1L)))/delta_x)

	;;; Initialise the second half of the abcissa with an constant lag increment
	tauHi = MIN([DOUBLE(N/2L)*delta_x,tauHi])
	PRINT,"New tauHi = ",tauHi
	PRINT,""

	tau_index[N_tau_samples/2L+1L:N_tau_samples] = tau_index[N_tau_samples/2L] + $
		FIX(((tauHi-tauMid)/(N_tau_samples/2L-1L))*FINDGEN(N_tau_samples/2L)/delta_x)

	;;; Remove duplicate indexed elements, reducing the number of samples accordingly
	tau_index_temp=tau_index[1L]
	FOR i=2L,N_tau_samples DO BEGIN
		IF (tau_index[i] EQ tau_index[i-1L]) THEN BEGIN
			FOR j=i,N_tau_samples-1L DO BEGIN
				tau_index[j] = tau_index[j+1L]
			ENDFOR
			N_tau_samples -= 1L
			i -= 1L
		ENDIF
		IF (N_tau_samples EQ i) THEN break
	ENDFOR

	tau_index = tau_index[0L:N_tau_samples]
	tauAbcissa = DBLARR(N_tau_samples+1L)
	tauAbcissa = tau[tau_index]

	RETURN

END



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Routine to return the glint cumulants from the glint moment(s), for 
;;; a number of geometries
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO glintCumulantsFromMoments,glintMoments,glintCumulants

	PRINT,'N_ELEMENTS: ',N_ELEMENTS(glintCumulants[*,0])
	FOR geometry=0L,N_ELEMENTS(glintCumulants[*,0])-1L DO BEGIN
		glintCumulants[geometry,0] = glintMoments[geometry]
		glintCumulants[geometry,1] = glintMoments[geometry]-glintMoments[geometry]^2.D
		glintCumulants[geometry,2] = glintMoments[geometry]-3.D*glintMoments[geometry]*glintMoments[geometry]+2.D*glintMoments[geometry]^3.D
	ENDFOR
;	glintCumulants[*,0] = glintMoments[*]
;	glintCumulants[*,1] = glintMoments[*]-glintMoments[*]^2.D
;	glintCumulants[*,2] = glintMoments[*]-3.D*glintMoments[*]*glintMoments[*]+2.D*glintMoments[*]^3.D
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Routine to return the slope cumulants from the slope moments
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO slopeCumulantsFromMoments,slopeMoments,slopeCumulants

	slopeCumulants[0] = slopeMoments[0]
	slopeCumulants[1] = slopeMoments[1]-slopeMoments[0]^2.D
	slopeCumulants[2] = slopeMoments[2]-3.D*slopeMoments[0]*slopeMoments[1]+2.D*slopeMoments[0]^3.D

END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Routine to use the symmetry properties of the bispectrum of a real 
;;; sequence to populate an entire NxN array from the primary octant. Takes
;;; as input an NxN complex array, and the array size N, and returns the 
;;; fully populated array in the input array
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO bispectrumSymmetry,bispectrum,N
	FOR j=0L,N/4L DO BEGIN
		FOR i=j,(N/2L-j) DO BEGIN
			bispectrum[(j LT 0L) ? N+(j) : j , (i LT 0L) ? N+(i) : i] = bispectrum[i,j]
			bispectrum[(j LT 0L) ? N+(j) : j , (-i-j LT 0L) ? N+(-i-j) : -i-j] = bispectrum[i,j]
			bispectrum[(-i-j LT 0L) ? N+(-i-j) : -i-j , (j LT 0L) ? N+(j) : j] = bispectrum[i,j]
			bispectrum[(-i-j LT 0L) ? N+(-i-j) : -i-j , (i LT 0L) ? N+(i) : i] = bispectrum[i,j]
			bispectrum[(i LT 0L) ? N+(i) : i , (-i-j LT 0L) ? N+(-i-j) : -i-j] = bispectrum[i,j]

			bispectrum[(-i LT 0L) ? N+(-i) : -i , (-j LT 0L) ? N+(-j) : -j   ] = CONJ(bispectrum[i,j])
			bispectrum[(-j LT 0L) ? N+(-j) : -j , (-i LT 0L) ? N+(-i) : -i   ] = CONJ(bispectrum[i,j])
			bispectrum[(-j LT 0L) ? N+(-j) : -j , (i+j LT 0L) ? N+(i+j) : i+j] = CONJ(bispectrum[i,j])
			bispectrum[(i+j LT 0L) ? N+(i+j) : i+j , (-j LT 0L) ? N+(-j) : -j] = CONJ(bispectrum[i,j])
			bispectrum[(i+j LT 0L) ? N+(i+j) : i+j , (-i LT 0L) ? N+(-i) : -i] = CONJ(bispectrum[i,j])
			bispectrum[(-i LT 0L) ? N+(-i) : -i , (i+j LT 0L) ? N+(i+j) : i+j] = CONJ(bispectrum[i,j])
		ENDFOR
	ENDFOR
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Routine to use the symmetry properties of the bicovariance of a real 
;;; sequence to populate an entire NxN array from the primary sextant. Takes
;;; as input an NxN complex array, and the array size N, and returns the 
;;; fully populated array in the input array
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO biCovarianceSymmetry,biCovariance,N
	FOR i=0L,N/2L DO BEGIN
		FOR j=0L,i DO BEGIN
			biCovariance[(j LT 0L) ? N+(j) : j , (i LT 0L) ? N+(i) : i]          = biCovariance[i,j]
			biCovariance[(-j LT 0L) ? N+(-j) : -j , (i-j LT 0L) ? N+(i-j) : i-j] = biCovariance[i,j]
			biCovariance[(i-j LT 0L) ? N+(i-j) : i-j , (-j LT 0L) ? N+(-j) : -j] = biCovariance[i,j]
			biCovariance[(j-i LT 0L) ? N+(j-i) : j-i , (-i LT 0L) ? N+(-i) : -i] = biCovariance[i,j]
			biCovariance[(-i LT 0L) ? N+(-i) : -i , (j-i LT 0L) ? N+(j-i) : j-i] = biCovariance[i,j]
		ENDFOR
	ENDFOR
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Open glint moment function data file for READ if exists, WRITE if not ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO getSlopeAndGlintRealisations,overWrite,Nruns,nlSwitch,SCALE,GEOM,POWER,NLCOUPLING,SLOPE,GLINT

	COMMON DEBUGGING, debugVar,windowIndex

	NN = SCALE.NN

	;;; Slope data arrays
	slopeData1= DBLARR(3L) 	;;; Slope Moments
	slopeData2 = DBLARR(SCALE.N/2L+1L) 	;;; Slope autocorrelation
	slopeData3= DBLARR(3L) 	;;; Slope Moments
	slopeData4 = DBLARR(NN/2L+1L,NN/2L+1L) 	;;; Slope bicovariances

	;;; Glint data arrays
	glintData1= DBLARR(GEOM.N_angles) 	;;; Glint First Moments
	glintData2 = DBLARR(GEOM.N_angles,SCALE.N/2L+1L) 	;;; Glint Autocorrelations
	glintData3= DBLARR(GEOM.N_angles) 	;;; Glint First Moments
	glintData4 = DBLARR(NN/2L+1L,NN/2L+1L) 	;;; Glint bicovariances

	;;; Spectrum Type
	spectrumType = POWER.spectrumType

	;;; Construct the data file names

	slopeDataFileName_M2 = 'M2slope_'+STRING(SCALE.N)+'_'+STRING(SCALE.NN)+'_'+$
	STRING(FIX(100*SCALE.delta_x))+'_'+STRING(Nruns)+'_'+spectrumType

	slopeDataFileName_M3 = 'M3slope_'+STRING(SCALE.N)+'_'+STRING(SCALE.NN)+'_'+$
	STRING(FIX(100*SCALE.delta_x))+'_'+STRING(Nruns)+'_'+spectrumType

	glintDataFileName_M2 = 'M2glint_'+STRING(SCALE.N)+'_'+STRING(SCALE.NN)+'_'+$
	STRING(FIX(100*SCALE.delta_x))+'_'+STRING(Nruns)+'_'+spectrumType

	glintDataFileName_M3 = 'M3glint_'+STRING(SCALE.N)+'_'+STRING(SCALE.NN)+'_'+$
	STRING(FIX(100*SCALE.delta_x))+'_'+STRING(Nruns)+'_'+spectrumType

	IF (nlSwitch EQ 1L) THEN BEGIN
		slopeDataFileName_M2 += '_nl.auto'
		slopeDataFileName_M3 += '_nl.auto'
		glintDataFileName_M2 += '_nl.auto'
		glintDataFileName_M3 += '_nl.auto'
	ENDIF ELSE BEGIN
		slopeDataFileName_M2 += '.auto'
		slopeDataFileName_M3 += '.auto'
		glintDataFileName_M2 += '.auto'
		glintDataFileName_M3 += '.auto'
	ENDELSE

	slopeDataFileName_M2 = STRCOMPRESS(slopeDataFileName_M2,/REMOVE_ALL)
	slopeDataFileName_M3 = STRCOMPRESS(slopeDataFileName_M3,/REMOVE_ALL)
	glintDataFileName_M2 = STRCOMPRESS(glintDataFileName_M2,/REMOVE_ALL)
	glintDataFileName_M3 = STRCOMPRESS(glintDataFileName_M3,/REMOVE_ALL)

	;;; Get the logical unit numbers for the data files
	PRINT,''
	GET_LUN,slopeDataFileLun_M2
	PRINT,'Got slope LUN:',slopeDataFileLun_M2
	GET_LUN,slopeDataFileLun_M3
	PRINT,'Got slope LUN:',slopeDataFileLun_M3
	GET_LUN,glintDataFileLun_M2
	PRINT,'Got glint LUN:',glintDataFileLun_M2
	GET_LUN,glintDataFileLun_M3
	PRINT,'Got glint LUN:',glintDataFileLun_M3

	slopeDataFileErr_M2 = 0
	slopeDataFileErr_M3 = 0

	glintDataFileErr_M2 = 0
	glintDataFileErr_M3 = 0

	;;; Open data files for reading
	OPENR, slopeDataFileLun_M2,slopeDataFileName_M2, ERROR = slopeDataFileErr_M2
	OPENR, slopeDataFileLun_M3,slopeDataFileName_M3, ERROR = slopeDataFileErr_M3
	OPENR, glintDataFileLun_M2,glintDataFileName_M2, ERROR = glintDataFileErr_M2
	OPENR, glintDataFileLun_M3,glintDataFileName_M3, ERROR = glintDataFileErr_M3

	IF (glintDataFileErr_M2 NE 0 OR slopeDataFileErr_M2 NE 0 OR $
		glintDataFileErr_M3 NE 0 OR slopeDataFileErr_M3 NE 0 OR $
		overWrite EQ 1) THEN BEGIN

		;;; Data file(s) do not exist or we want to create new ones

		IF (glintDataFileErr_M2 NE 0) THEN $
			PRINT,'There was a problem opening file '+glintDataFileName_M2
		IF (slopeDataFileErr_M2 NE 0) THEN $
			PRINT,'There was a problem opening file '+slopeDataFileName_M2
		PRINT,'Opening '+slopeDataFileName_M2+' for writing'
		PRINT,'Opening '+glintDataFileName_M2+' for writing'

		IF (glintDataFileErr_M3 NE 0) THEN $
			PRINT,'There was a problem opening file '+glintDataFileName_M3
		IF (slopeDataFileErr_M3 NE 0) THEN $
			PRINT,'There was a problem opening file '+slopeDataFileName_M3
		PRINT,'Opening '+slopeDataFileName_M3+' for writing'
		PRINT,'Opening '+glintDataFileName_M3+' for writing'

		CLOSE,slopeDataFileLun_M2
		CLOSE,glintDataFileLun_M2
		OPENW,slopeDataFileLun_M2,slopeDataFileName_M2
		OPENW,glintDataFileLun_M2,glintDataFileName_M2

		CLOSE,slopeDataFileLun_M3
		CLOSE,glintDataFileLun_M3
		OPENW,slopeDataFileLun_M3,slopeDataFileName_M3
		OPENW,glintDataFileLun_M3,glintDataFileName_M3

		;;; Simulate the first glint moment and second and third glint moment functions, 
		;;; and first slope moment, second slope moment and third slope moment
		;;; function,for one or more geometries definied in GEOM

;		Nruns = 50L 		;;; TEMP CHANGE

		FOR geometry=0L,GEOM.N_angles-1L DO BEGIN
			GLINT.geomIndex=geometry
			glintMomentFunctionSimulate,SCALE,GEOM,SLOPE,GLINT,POWER,NLCOUPLING,Nruns,nlSwitch
		ENDFOR

		;;; Write the slope moments, and the second slope moment function, to file
		PRINTF,slopeDataFileLun_M2,SLOPE.moments
		FOR i=0L,SCALE.N/2L DO BEGIN
			PRINTF,slopeDataFileLun_M2,SLOPE.secondMomentFunction[i]
		ENDFOR
		PRINT,'Closing file '+slopeDataFileName_M2
		CLOSE,slopeDataFileLun_M2
		PRINT,'Freeing LUN ',slopeDataFileLun_M2
		FREE_LUN,slopeDataFileLun_M2

		;;; Write the first glint moment, and the second glint moment function, to file
		;;; for each geometry.
		PRINTF,glintDataFileLun_M2,GLINT.firstMoment
		PRINTF,glintDataFileLun_M2,GLINT.secondMomentFunction
		PRINT,'Closing file '+glintDataFileName_M2
		CLOSE,glintDataFileLun_M2
		PRINT,'Freeing LUN ',glintDataFileLun_M2
		FREE_LUN,glintDataFileLun_M2

		;;; Write the third slope moment, and the third slope moment function, to file
		PRINTF,slopeDataFileLun_M3,SLOPE.moments
		formatStrng='('+STRING(SCALE.NN/2L+1L)+'E15.6)'
		formatStrng = STRCOMPRESS(formatStrng,/REMOVE_ALL)
		PRINT,'Format = ',formatStrng
		PRINTF,slopeDataFileLun_M3,DOUBLE(SLOPE.thirdMomentFunction),FORMAT=formatStrng
		PRINT,'Closing file '+slopeDataFileName_M3
		CLOSE,slopeDataFileLun_M3
		PRINT,'Freeing LUN ',slopeDataFileLun_M3
		FREE_LUN,slopeDataFileLun_M3

		;;; Write the third glint moment, and the third glint moment function, to file
		;;; for each geometry.
		PRINT,'NN = ',SCALE.NN
		formatStrng='('+STRING(SCALE.NN/2L+1L)+'E15.6)'
		formatStrng = STRCOMPRESS(formatStrng,/REMOVE_ALL)
		PRINT,'Format = ',formatStrng
		PRINTF,glintDataFileLun_M3,GLINT.firstMoment

;		FOR geometry=2L,2L DO BEGIN
		FOR geometry=0L,GEOM.N_angles-1L DO BEGIN
			GLINT.geomIndex=geometry
			PRINTF,glintDataFileLun_M3, $
				DOUBLE(GLINT.thirdMomentFunctionGeom[GLINT.geomIndex].thirdMomentFunction),FORMAT=formatStrng
		ENDFOR
		PRINT,'Closing file '+glintDataFileName_M3
		CLOSE,glintDataFileLun_M3
		PRINT,'Freeing LUN ',glintDataFileLun_M3
		FREE_LUN,glintDataFileLun_M3

	ENDIF ELSE BEGIN

		;;; Data file(s) exist and we do not want to create new ones
		PRINT,'Success opening file '+slopeDataFileName_M2
		PRINT,'Success opening file '+slopeDataFileName_M3
		PRINT,'Success opening file '+glintDataFileName_M2
		PRINT,'Success opening file '+glintDataFileName_M3

		;;; Read in the slope moments and the slope second moment function
		READF,slopeDataFileLun_M2,slopeData1,slopeData2
		SLOPE.moments = slopeData1
		SLOPE.secondMomentFunction = slopeData2
		PRINT,'Closing file '+slopeDataFileName_M2
		CLOSE,slopeDataFileLun_M2
		PRINT,'Freeing LUN ',slopeDataFileLun_M2
		FREE_LUN,slopeDataFileLun_M2

		;;; Read in the third slope moment and third slope moment function
		READF,slopeDataFileLun_M3,slopeData3,slopeData4
		SLOPE.moments = slopeData3
		SLOPE.thirdMomentFunction = slopeData4
		PRINT,'Closing file '+slopeDataFileName_M3
		CLOSE,slopeDataFileLun_M3
		PRINT,'Freeing LUN ',slopeDataFileLun_M3
		FREE_LUN,slopeDataFileLun_M3

		;;; Read in first glint moments and second moment functions for each geometry 
		READF,glintDataFileLun_M2,glintData1,glintData2
		GLINT.firstMoment = glintData1
		GLINT.secondMomentFunction = glintData2
		PRINT,'Closing file '+glintDataFileName_M2
		CLOSE,glintDataFileLun_M2
		PRINT,'Freeing LUN ',glintDataFileLun_M2
		FREE_LUN,glintDataFileLun_M2

		;;; Read in the third glint moments and the third glint moment functions for each geometry
		READF,glintDataFileLun_M3,glintData3
;		GLINT.firstMoment = glintData3

		FOR geometry=0L,GEOM.N_angles-1L DO BEGIN
			GLINT.geomIndex=geometry
			READF,glintDataFileLun_M3,glintData4
			GLINT.thirdMomentFunctionGeom[GLINT.geomIndex].thirdMomentFunction = glintData4
		ENDFOR
		PRINT,'Closing file '+glintDataFileName_M3
		CLOSE,glintDataFileLun_M3
		PRINT,'Freeing LUN ',glintDataFileLun_M3
		FREE_LUN,glintDataFileLun_M3

	ENDELSE

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;; Do a couple of plots                     ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;	tau_x = DBLARR(SCALE.NN)
;	tau_x = DINDGEN(SCALE.NN)*SCALE.delta_x
;	tau_y = DBLARR(SCALE.NN)
;	tau_y = DINDGEN(SCALE.NN)*SCALE.delta_x
;
;	newThirdMomentFunction = DBLARR(NN,NN)
;
;	DEVICE, DECOMPOSED=0, RETAIN=2
;	LOADCT,39
;	CONTOUR_LEVELS=20L
;
;	rotate_z=-20.
;	rotate_x=10.
;	z_axis=2
;	ysize=450
;	xsize=ROUND(1.2*ysize)
;; 	!P.MULTI=0
;	!P.MULTI = [0,2,1]
;;
;	newThirdMomentFunction[0L:NN/2L,0L:NN/2L] = SLOPE.thirdMomentFunction[0L:NN/2L,0L:NN/2L]
;	PRINT,"Value of Third Slope Moment Function at Zero Lag: ",newThirdMomentFunction[0L,0L]
;;
;	biCovarianceSymmetry,newThirdMomentFunction,NN
;	newSlopeBispectrum = COMPLEXARR(NN,NN)
;	newSlopeBispectrum = FFT(newThirdMomentFunction,/DOUBLE)
;
;	newThirdMomentFunction = SHIFT(newThirdMomentFunction,NN/2-1L,NN/2-1L)
;	windowIndex++
;	window,windowIndex,xsize = 2*xsize,ysize = ysize,title='Total New Uncentered Real Third Slope Moment Function',RETAIN=2
;	bispecMin = MIN([DOUBLE(newThirdMomentFunction),IMAGINARY(newThirdMomentFunction)])
;	bispecMax = MAX([DOUBLE(newThirdMomentFunction),IMAGINARY(newThirdMomentFunction)])
;	SURFACE,DOUBLE(newThirdMomentFunction),AZ=rotate_z,$
;;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax],$
;		zrange=[5.*bispecMin,5.*bispecMax], $
;		/xstyle,/ystyle,/zstyle,$
;		charsize=2.5,xtitle='k_x/k_N',ytitle='k_y/k_N',ZAXIS=z_axis
;	contourMin=MIN(DOUBLE(newThirdMomentFunction))
;	contourMax=MAX(DOUBLE(newThirdMomentFunction))
;	contourRange=contourMax-contourMin
;	CONTOUR,DOUBLE(newThirdMomentFunction), $
;;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax], $
;		xtitle='x',ytitle='y',/xstyle,/ystyle,/zstyle, $
;		/CELL_FILL,/ISOTROPIC,LEVELS=contourMin+contourRange*FINDGEN(20L)/DOUBLE(19L)
;	PRINT,'Maximum Value of Third Slope Moment Function: ',MAX(DOUBLE(newThirdMomentFunction))
;	PRINT,'Minimum Value of Third Slope Moment Function: ',MIN(DOUBLE(newThirdMomentFunction))
;
;
;	ysize=450
;	xsize=ROUND(1.2*ysize)
;; 	!P.MULTI=0
;	!P.MULTI = [0,2,1]
;
;	windowIndex++
;	window,windowIndex,xsize = 2*xsize,ysize = ysize,title='Total New Slope real bispectrum',RETAIN=2
;	bispecMin = MIN([DOUBLE(newSlopeBispectrum),IMAGINARY(newSlopeBispectrum)])
;	bispecMax = MAX([DOUBLE(newSlopeBispectrum),IMAGINARY(newSlopeBispectrum)])
;	SURFACE,DOUBLE(SHIFT(newSlopeBispectrum,NN/2L-1L,NN/2L-1L)),AZ=rotate_z,$
;;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax],$
;		zrange=[bispecMin,bispecMax], $
;		/xstyle,/ystyle,/zstyle,$
;		charsize=2.5,xtitle='k_x/k_N',ytitle='k_y/k_N',ZAXIS=z_axis
;	contourMin=MIN(DOUBLE(newSlopeBispectrum))
;	contourMax=MAX(DOUBLE(newSlopeBispectrum))
;	contourRange=contourMax-contourMin
;	CONTOUR,DOUBLE(SHIFT(newSlopeBispectrum,NN/2L-1L,NN/2L-1L)), $
;;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax], $
;		xtitle='k_x/k_N',ytitle='k_y/k_N',/xstyle,/ystyle,/zstyle, $
;		/CELL_FILL,/ISOTROPIC,LEVELS=contourMin+contourRange*FINDGEN(CONTOUR_LEVELS)/DOUBLE(CONTOUR_LEVELS-1L)
;	PRINT,'Maximum Value of Total Slope Real Bispectrum: ',MAX([DOUBLE(newSlopeBispectrum)])
;	PRINT,'Minimum Value of Total Slope Real Bispectrum: ',MIN([DOUBLE(newSlopeBispectrum)])
;
;	ysize=450
;	xsize=ROUND(1.2*ysize)
;; 	!P.MULTI=0
;	!P.MULTI = [0,2,1]
;
;	windowIndex++
;	window,windowIndex,xsize = 2*xsize,ysize = ysize,title='Total New Slope imaginary bispectrum',RETAIN=2
;	SURFACE,IMAGINARY(SHIFT(newSlopeBispectrum,NN/2L-1L,NN/2L-1L)),AZ=rotate_z,$
;;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax],$
;		zrange=[bispecMin,bispecMax], $
;		/xstyle,/ystyle,/zstyle,$
;		charsize=2.5,xtitle='k_x/k_N',ytitle='k_y/k_N',ZAXIS=z_axis
;	contourMin=MIN(IMAGINARY(newSlopeBispectrum))
;	contourMax=MAX(IMAGINARY(newSlopeBispectrum))
;	contourRange=contourMax-contourMin
;	CONTOUR,IMAGINARY(SHIFT(newSlopeBispectrum,NN/2L-1L,NN/2L-1L)), $
;;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax], $
;		xtitle='k_x/k_N',ytitle='k_y/k_N',/xstyle,/ystyle,/zstyle, $
;		/CELL_FILL,/ISOTROPIC,LEVELS=contourMin+contourRange*FINDGEN(CONTOUR_LEVELS)/DOUBLE(CONTOUR_LEVELS-1L)
;	PRINT,'Maximum Value of Total Slope Imaginary Bispectrum: ',MAX([IMAGINARY(newSlopeBispectrum)])
;	PRINT,'Minimum Value of Total Slope Imaginary Bispectrum: ',MIN([IMAGINARY(newSlopeBispectrum)])
;
;	ysize=450
;	xsize=ROUND(1.2*ysize)
;; 	!P.MULTI=0
;	!P.MULTI = [0,2,1]
;
;	newThirdMomentFunction[*,*] = 0.D
;
;;	FOR geometry=2L,2L DO BEGIN
;	FOR geometry=0L,GEOM.N_angles-1L DO BEGIN
;		GLINT.geomIndex=geometry
;		
;		newThirdMomentFunction[0L:NN/2L,0L:NN/2L] = $
;			GLINT.thirdMomentFunctionGeom[GLINT.geomIndex].thirdMomentFunction[0L:NN/2L,0L:NN/2L]
;		PRINT,"Value of Third Glint Moment Function at Zero Lag: ",newThirdMomentFunction[0L,0L]
;		biCovarianceSymmetry,newThirdMomentFunction,NN
;
;		newThirdMomentFunction = SHIFT(newThirdMomentFunction,NN/2L-1L,NN/2L-1L)
;		windowIndex++
;		window,windowIndex,xsize = 2*xsize,ysize = ysize, $
;			title='Total Uncentered Real Third Glint Moment Function: Geometry '+STRCOMPRESS(STRING(geometry),/REMOVE_ALL),RETAIN=2
;		bispecMin = MIN([DOUBLE(newThirdMomentFunction), IMAGINARY(newThirdMomentFunction)])
;		bispecMax = MAX([DOUBLE(newThirdMomentFunction), IMAGINARY(newThirdMomentFunction)])
;		SURFACE,DOUBLE(newThirdMomentFunction),AZ=rotate_z,$
;	;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax],$
;			zrange=[bispecMin,bispecMax], $
;			/xstyle,/ystyle,/zstyle,$
;			charsize=2.5,xtitle='k_x/k_N',ytitle='k_y/k_N',ZAXIS=z_axis
;		contourMin=MIN(DOUBLE(newThirdMomentFunction))
;		contourMax=MAX(DOUBLE(newThirdMomentFunction))
;		contourRange=contourMax-contourMin
;		CONTOUR,DOUBLE(newThirdMomentFunction), $
;	;		xrange=[k_PlotMin,k_PlotMax],yrange=[k_PlotMin,k_PlotMax], $
;			xtitle='x',ytitle='y',/xstyle,/ystyle,/zstyle, $
;			/CELL_FILL,/ISOTROPIC,LEVELS=contourMin+contourRange*FINDGEN(20L)/DOUBLE(19L)
;		PRINT,'Maximum Value of Third Glint Moment Function: ',MAX(DOUBLE(newThirdMomentFunction))
;		PRINT,'Minimum Value of Third Glint Moment Function: ',MIN(DOUBLE(newThirdMomentFunction))
;	ENDFOR
END
