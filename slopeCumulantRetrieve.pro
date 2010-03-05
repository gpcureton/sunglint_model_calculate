;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; Subroutine to retrieve the slope cumulants, for glint data at the geometries
;;; specified in GEOM

	;;; GEOM: Structure containing angular dependent parameters
	;;;		for a range of geometries.
	;;; mu_L1: Vector containing the first glint moment for each 
	;;;		geometry in GEOM.
	;;; kappa_xi: Vector containing the retrieved slope cumulants.
	;;; sigma: The standard deviations of the retrieved slope cumulants.
	;;; chiSq: The reduced chi-square goodness-of-fit statistic.
	;;; mu_L1_fit: Vector containing the  fitted values of mu_L1 
	;;;		computed from the retrieved slope cumulants for each 
	;;;		geometry in GEOM.

PRO slopeCumulant_retrieve,GEOM,mu_L1,kappa_xi,fitted,sigma,chiSq,mu_L1_fit

	weights = REPLICATE(1.,GEOM.N_angles)

	mu_L1_fit = CURVEFIT(GEOM, mu_L1, weights,kappa_xi,$
	FITA=fitted,sigma,CHISQ=chiSq,/DOUBLE,FUNCTION_NAME='firstGlintMoment')

END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; Subroutine to retrieve slope cumulant functions, for glint data at the geometries
;;; specified in GEOM

	;;; GEOM: Structure containing angular dependent parameters
	;;;		for a range of geometries.
	;;; M2_L: Vector containing the second glint moment function for each 
	;;;		geometry in GEOM.
	;;; slopeCumulantFunctions: Vector containing the slope cumulant functions,
	;;;		for a single lag value
	;;; sigma: The standard deviations of the retrieved slope cumulant functions
	;;; chiSq: The reduced chi-square goodness-of-fit statistic.
	;;; M2_L_fit: Vector containing the  fitted values of M_L2 
	;;;		computed from the retrieved slope cumulant functions, for each 
	;;;		geometry in GEOM.

PRO slopeCumulantFunction_retrieve,GEOM,M2_L,slopeCumulantFunctions,$
	fitted,sigma,chiSq,M2_L_fit

	weights = REPLICATE(1.,GEOM.N_angles)

	M2_L_fit = CURVEFIT(GEOM,M2_L,weights,slopeCumulantFunctions,$
	FITA=fitted,sigma,CHISQ=chiSq,/DOUBLE,FUNCTION_NAME='secondGlintMomentFunction')

END

PRO slopeCumulantFunction_retrieve2,GEOM,M2_L,slopeCumulantFunctions,$
	fitted,sigma,chiSq,M2_L_fit

	weights = REPLICATE(1.,GEOM.N_angles)

	M2_L_fit = CURVEFIT(GEOM,M2_L,weights,slopeCumulantFunctions,$
	FITA=fitted,sigma,CHISQ=chiSq,/DOUBLE,/NODERIVATIVE,FUNCTION_NAME='secondGlintMomentFunction2')

END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Subroutine to retrieve third slope cumulant functions, for glint data at 
;;; the geometries specified in GEOM
;;;
;;;		GEOM: Structure containing angular dependent parameters
;;;			for a range of geometries.
;;;		M3_L: Vector containing the third glint moment function for each 
;;;			geometry in GEOM.
;;;		slopeCumulantFunctions: Vector containing the slope cumulant functions,
;;;			for a single lag value
;;;		sigma: The standard deviations of the retrieved slope cumulant functions
;;;		chiSq: The reduced chi-square goodness-of-fit statistic.
;;;		M3_L_fit: Vector containing the  fitted values of M_L3 
;;;			computed from the retrieved slope cumulant functions, for each 
;;;			geometry in GEOM.

PRO slopeCumulantFunction_retrieve21,GEOM,M21_L,slopeCumulantFunctions,$
	fitted,sigma,chiSq,M21_L_fit

	weights = REPLICATE(1.,GEOM.N_angles)
	M21_L_fit = CURVEFIT(GEOM,M21_L,weights,slopeCumulantFunctions,$
	FITA=fitted,sigma,CHISQ=chiSq,/DOUBLE,/NODERIVATIVE,FUNCTION_NAME='thirdGlintMomentFunction21')

END

;;;;;;;;;

PRO slopeCumulantFunction_retrieve12,GEOM,M12_L,slopeCumulantFunctions,$
	fitted,sigma,chiSq,M12_L_fit

	weights = REPLICATE(1.,GEOM.N_angles)

	M12_L_fit = CURVEFIT(GEOM,M12_L,weights,slopeCumulantFunctions,$
	FITA=fitted,sigma,CHISQ=chiSq,/DOUBLE,/NODERIVATIVE,FUNCTION_NAME='thirdGlintMomentFunction12')

END

;;;;;;;;

PRO slopeThirdCumulantFunction_retrieve,GEOM,M3_L,slopeCumulantFunctions,$
	fitted,sigma,chiSq,M3_L_fit

	weights = REPLICATE(1.,GEOM.N_angles)

	M3_L_fit = CURVEFIT(GEOM,M3_L,weights,slopeCumulantFunctions,$
	FITA=fitted,sigma,CHISQ=chiSq,/DOUBLE,/NODERIVATIVE,FUNCTION_NAME='thirdGlintMomentFunction2')
END
