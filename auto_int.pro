;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; The purpose of this program is to convert the retrieved slope autocorrelation
; into a power spectrum, and do the same with the simulated slope
; autocorrelation.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO auto_int,kappa_eta_2,c2_xi,c2_eta,steps,dtau

	COMMON DIFFBLOCK, deltaTau, dydx_tab, tau

	deltaTau = dtau

	;;; Assign the derivative
	dydx_tab = DBLARR(steps+1L)

	tau = DINDGEN(steps+1L)*dtau
	p = DBLARR(steps+1L)

	;FOR j=0L, 20L DO BEGIN
		dydx_tab = -c2_xi -0.0001001

		p[0L] = -0.00000999
		;p[0L] = -0.00001 + DOUBLE(j)*0.000002

		; Integrate the slope covariance function to obtain p
		FOR i=1L, steps DO BEGIN
			; Integrate over the interval (tau[i-1], tau[i]):
			p[i] = RK4(p[i-1], dydx_tab[i-1], tau[i-1], dtau, 'difftab',/DOUBLE)
		ENDFOR

		;;; Assign the derivative
		dydx_tab = p

		;;; The initial conditions for c2_eta
		c2_eta[0L] = 0.D

		; Integrate the slope covariance function to obtain p
		FOR i=1L, steps DO BEGIN
			; Integrate over the interval (tau[i-1], tau[i]):
			c2_eta[i] = RK4(c2_eta[i-1], dydx_tab[i-1], tau[i-1], dtau, 'difftab',/DOUBLE)
		ENDFOR

		c2_eta = c2_eta - c2_eta[steps]

		PRINT,'p[0]=',p[0],'   p[inf]=',p[steps],'  Elevation Variance: ',c2_eta[0]
	;END

	kappa_eta_2 = c2_eta[0]
	c2_eta = c2_eta/c2_eta[0L]

	xwinsize=800
	ywinsize=450

	; Plot the original and solution vectors
	WINDOW, 10, xsize=xwinsize,ysize=ywinsize,title="C2_xi_ret",RETAIN=2
	PLOT,tau,c2_xi/c2_xi[0],yrange=[-0.2,0.2],/xstyle,/ystyle

	WINDOW, 11, xsize=xwinsize,ysize=ywinsize,title="p",RETAIN=2
	PLOT,tau,p,/xstyle

	WINDOW, 12, xsize=xwinsize,ysize=ywinsize,title="C2_eta_ret",RETAIN=2
	PLOT, tau,c2_eta,/xstyle
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; This function takes as input the scalar value x, and 
;;; determines the indicies of xabciss which bound that value.
;;; The values of dydx_tab at the same indicies are then
;;; used to interpolate a value of dydx interior of these
;;; indices.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
FUNCTION difftab, X, Y
	COMMON DIFFBLOCK, H, dydx_tab, xabciss

	; Determine the indicies of the abcissa which bound the supplied value X 

	range = X-xabciss[0]
	fIncrements = range/H
	iIncrements = FLOOR(fIncrements)
	dInc = fIncrements-DOUBLE(iIncrements)

	lowDiff = 10.D^(-10.D)
	highDiff = 1.D - 1.D/(10.D^10.D)

	IF( (dInc LT lowDiff) OR (dInc GT highDiff) )THEN BEGIN
		IF (dInc LT lowDiff) THEN BEGIN
			dydx = dydx_tab(iIncrements)
			RETURN, dydx
		ENDIF
		IF (dInc GT highDiff) THEN BEGIN
			dydx = dydx_tab(iIncrements+1L)
			RETURN, dydx
		ENDIF
	ENDIF ELSE BEGIN
		; Determine the values of the derivative at the bounds
		dydx_low = dydx_tab[iIncrements]
		dydx_hi = dydx_tab[iIncrements+1L]
		
		; Slope of the derivative within the interval
		gradient = (dydx_hi-dydx_low)/H
		
		; Calculate the derivative at X
		dydx = dydx_low + gradient*(dInc*H)
		
		; Return the value of the derivative
		RETURN, dydx
	ENDELSE
END
