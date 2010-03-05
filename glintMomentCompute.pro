
;;; Function to return the integration limits for xi2
FUNCTION PQ_Limits,x
	COMMON GEOMCOMMON, xiMin, xiMax
	RETURN,[xiMin,xiMax]
END

;;; Function to return the integration limits for xi3
FUNCTION UV_Limits,x,y
	COMMON GEOMCOMMON, xiMin, xiMax
	RETURN,[xiMin,xiMax]
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; Procedure to calculate the first glint moment and it's derivatives with
;;; respect to each of the slope cumulants in kappa_xi[]

PRO firstGlintMoment,GEOM,kappa_xi,mu_L1,pderiv

	;;; GEOM: Structure containing angular dependent parameters
	;;;		for a range of geometries.
	;;; kappa_xi: Vector containing the slope cumulants.
	;;; mu_L1: Vector containing the first glint moment for each 
	;;;		geometry in GEOM.
	;;; pderiv: Vector containing the derivatives of mu_L1 with 
	;;;		respect to each of the slope cumulants 
	;;;		in kappa_xi.
	;;; xi: The slope variable
	
	d2r = !DPI/180.D
	r2d = 180.D/!DPI	

	lambda_xi_3 = kappa_xi[2]/(kappa_xi[1]^(1.5D)) 

	nu_lo = (GEOM.xi_min-kappa_xi[0])/(SQRT(kappa_xi[1]))
	nu_hi = (GEOM.xi_max-kappa_xi[0])/(SQRT(kappa_xi[1]))
	
	;;; Calculate the glint moment mu_L1
	
	I_1 = 0.5D*(ERRORF(nu_hi/SQRT(2.D))-ERRORF(nu_lo/SQRT(2.D)))
	I_2=-(lambda_xi_3/(6.D*SQRT(2.D*!DPI)))*$
	(He(2,nu_hi)*EXP(-nu_hi*nu_hi/2.D)-He(2,nu_lo)*EXP(-nu_lo*nu_lo/2.D))
    
	mu_L1 = I_1 + I_2

	;;; Calculate the derivatives of mu_L1 with respect to kappa_xi[0],
	;;; kappa_xi[1] and kappa_xi[2]...
	
	mu_L1_D1=-(1.D/SQRT(!DPI*kappa_xi[1]))*(EXP(-nu_hi^2.D)-EXP(-nu_lo^2.D))- $
		(1.D/SQRT(2.D*!DPI*kappa_xi[1]))*(lambda_xi_3/6.D)*$
		(He(1,nu_hi)*(He(2,nu_hi)-2.D)*EXP(-nu_hi^2.D)-$
		He(1,nu_lo)*(He(2,nu_lo)-2.D)*EXP(-nu_lo^2.D))
		
	mu_L1_D2=-(1.D/(SQRT(4.D*!PI)*kappa_xi[1]))*$
		(He(1,nu_hi)*EXP(-nu_hi^2.D)-He(1,nu_lo)*EXP(-nu_lo^2.D))- $
		(lambda_xi_3/(6.D*SQRT(8.D*!DPI)*kappa_xi[1]))*$
		(He(1,nu_hi)*(He(2,nu_hi)-2.D)*EXP(-nu_hi^2.D)-$
		He(1,nu_lo)*(He(2,nu_lo)-2.D)*EXP(-nu_lo^2.D))
	
	mu_L1_D3 = -(1.D/(6.D*SQRT(2.D*!DPI)*(kappa_xi[1]^1.5D)))*$
		(He(2,nu_hi)*EXP(-nu_hi^2.D)-He(2,nu_lo)*EXP(-nu_lo^2.D))
	
	pderiv = [[mu_L1_D1],[mu_L1_D2],[mu_L1_D3]]

END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; Function to calculate the second glint moment function (nonlinear) for a 
;;; particular value of the lag

PRO secondGlintMomentFunction,GEOM,kappa_xi,M2_L,pderiv

	;;; GEOM: Structure containing angular dependent parameters
	;;;		for a range of geometries.
	;;; kappa_xi: Vector containing the slope cumulant functions,
	;;;		kappa_xi_11 and kappa_xi_21 for a single lag value.
	;;; M2_L: Vector containing the second glint moment function for each 
	;;;		geometry in GEOM, for a single lag value
	;;; pderiv: Vector containing the derivatives of M2_L with 
	;;;		respect to each of the slope cumulants and cumulant functions

	COMMON SLOPECUMULANTS,kappa_xi_1,kappa_xi_2,kappa_xi_3,kappa_xi_11,kappa_xi_21,kappa_xi_12, $
		kappa_xi_110,kappa_xi_101,kappa_xi_011, $
		kappa_xi_210,kappa_xi_120, $
		kappa_xi_201,kappa_xi_102, $
		kappa_xi_012,kappa_xi_021, $
		kappa_xi_111

	COMMON GLINTMOMENTS,mu_L_1,mu_L_2,mu_L_3,mu_L_11,mu_L_21,mu_L_12,mu_L_111
	COMMON LAMBDA,lambda2_xi_Matrix,lambda3_xi_Matrix,lambda_xi_Det,lambda2_xi_invMatrix,lambda3_xi_invMatrix
	COMMON GEOMCOMMON, xiMin, xiMax
	COMMON DERIV,derivIndex

	;;; Initialise the normalised slope cumulant functions from the input parameters
	;;; and the SLOPECUMULANTS common block. They only need to be computed once for each time
	;;; secondGlintMomentFunction is called

	kappa_xi_11 = kappa_xi[0]
	kappa_xi_21 = kappa_xi[1]

	lambda_xi_11 = kappa_xi_11/kappa_xi_2
	lambda_xi_21 = kappa_xi_21/(kappa_xi_2^(1.5D))

	lambda2_xi_Matrix = [[1.D,lambda_xi_11],[lambda_xi_11,1.D]]
	lambda_xi_Det = 1.D - lambda_xi_11^2.D
	lambda2_xi_invMatrix = (lambda_xi_Det GT 0.000001D) ? $
		[[1.D,-lambda_xi_11],[-lambda_xi_11,1.D]]/lambda_xi_Det : $
		[[0.D,-1.D],[-1.D,0.D]]/(2.D*lambda_xi_11)

	d2r = !DPI/180.D
	r2d = 180.D/!DPI	

	M2_L = DBLARR(GEOM.N_angles)

	;;; Calculate the glint moment function M2_L for each geometry

	M2_L_D11 = DBLARR(N_ELEMENTS(GEOM.xi_min))
	M2_L_D21 = DBLARR(N_ELEMENTS(GEOM.xi_min))

	FOR geometry = 0L,N_ELEMENTS(GEOM.xi_min)-1L DO BEGIN

		AB_Limits = [GEOM.xi_min[geometry],GEOM.xi_max[geometry]]
		xiMin = AB_Limits[0]
		xiMax = AB_Limits[1]

		M2_L[geometry] = (lambda_xi_Det GT 0.000001D) ? $
			(1.D/mu_L_2[geometry])*INT_2D('prob_xi_int',AB_Limits,'PQ_Limits',48,/DOUBLE): $
			1.D

		;;; Calculate the derivatives of M2_L with respect to kappa_xi[0],
		;;; kappa_xi[1] and kappa_xi[2], kappa_xi[3] and kappa_xi[4]...

		derivIndex=11
		M2_L_D11[geometry] = (1.D/mu_L_1[geometry]) $
			* INT_2D('prob_xi_deriv_int',AB_Limits,'PQ_Limits',48,/DOUBLE)

		derivIndex=21
		M2_L_D21[geometry] = (1.D/mu_L_1[geometry]) $
			* INT_2D('prob_xi_deriv_int',AB_Limits,'PQ_Limits',48,/DOUBLE)

		pderiv = [[M2_L_D11],[M2_L_D21]]
	ENDFOR
END


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; Function to calculate the second glint moment function (nonlinear) for a 
;;; particular value of the lag. This version sets the derivatives to unity,
;;; and the calling function computes the derivatives from differences

PRO secondGlintMomentFunction2,GEOM,kappa_xi,M2_L

	;;; GEOM: Structure containing angular dependent parameters
	;;;		for a range of geometries.
	;;; kappa_xi: Vector containing the slope cumulant functions,
	;;;		kappa_xi_11, kappa_xi_21 and kappa_xi_12 for a single lag value.
	;;; M2_L: Vector containing the second glint moment function for each 
	;;;		geometry in GEOM, for a single lag value

	COMMON SLOPECUMULANTS,kappa_xi_1,kappa_xi_2,kappa_xi_3,kappa_xi_11,kappa_xi_21,kappa_xi_12, $
		kappa_xi_110,kappa_xi_101,kappa_xi_011, $
		kappa_xi_210,kappa_xi_120, $
		kappa_xi_201,kappa_xi_102, $
		kappa_xi_012,kappa_xi_021, $
		kappa_xi_111

	COMMON GLINTMOMENTS,mu_L_1,mu_L_2,mu_L_3,mu_L_11,mu_L_21,mu_L_12,mu_L_111
	COMMON LAMBDA,lambda2_xi_Matrix,lambda3_xi_Matrix,lambda_xi_Det,lambda2_xi_invMatrix,lambda3_xi_invMatrix
	COMMON GEOMCOMMON, xiMin, xiMax
	COMMON DERIV,derivIndex

	;;; Initialise the normalised slope cumulant functions from the input parameters
	;;; and the SLOPECUMULANTS common block. They only need to be computed once for each time
	;;; secondGlintMomentFunction is called

	kappa_xi_11 = kappa_xi[0]
	kappa_xi_21 = kappa_xi[1]
	kappa_xi_12 = kappa_xi[2]

	lambda_xi_11 = kappa_xi_11/kappa_xi_2
	lambda_xi_21 = kappa_xi_21/(kappa_xi_2^(1.5D))
	lambda_xi_12 = kappa_xi_12/(kappa_xi_2^(1.5D))

	lambda_xi_Det = 1.D - lambda_xi_11^2.D

	lambda2_xi_invMatrix = 	[[1.D,-lambda_xi_11],[-lambda_xi_11,1.D]]/lambda_xi_Det
	
	d2r = !DPI/180.D
	r2d = 180.D/!DPI	

	M2_L = DBLARR(GEOM.N_angles)

	;;; Calculate the glint moment function M2_L for each geometry

	FOR geometry = 0L,N_ELEMENTS(GEOM.xi_min)-1L DO BEGIN

		AB_Limits = [GEOM.xi_min[geometry],GEOM.xi_max[geometry]]
		xiMin = AB_Limits[0]
		xiMax = AB_Limits[1]

;		M2_L[geometry] = (lambda_xi_Det GT 0.000001D) ? $
;			(1.D/mu_L_1[geometry])*INT_2D('prob2D_xi_int',AB_Limits,'PQ_Limits',6,/DOUBLE): $
;			1.D
		M2_L[geometry] = (1.D/mu_L_1[geometry])*INT_2D('prob2D_xi_int',AB_Limits,'PQ_Limits',6,/DOUBLE)

	ENDFOR
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO thirdGlintMomentFunction21,GEOM,kappa_xi,M2_L,M21_L

	;;; GEOM: Structure containing angular dependent parameters
	;;;		for a range of geometries.
	;;; kappa_xi: Vector containing the slope cumulant functions,
	;;;		kappa_xi_11, kappa_xi_21 and kappa_xi_12 for a single lag value.
	;;; M21_L: Vector containing the third glint moment function for each 
	;;;		geometry in GEOM, along the tau_1 axis, for a single lag value

	COMMON SLOPECUMULANTS,kappa_xi_1,kappa_xi_2,kappa_xi_3,kappa_xi_11,kappa_xi_21,kappa_xi_12, $
		kappa_xi_110,kappa_xi_101,kappa_xi_011, $
		kappa_xi_210,kappa_xi_120, $
		kappa_xi_201,kappa_xi_102, $
		kappa_xi_012,kappa_xi_021, $
		kappa_xi_111

	COMMON GLINTMOMENTS,mu_L_1,mu_L_2,mu_L_3,mu_L_11,mu_L_21,mu_L_12,mu_L_111
	COMMON LAMBDA,lambda2_xi_Matrix,lambda3_xi_Matrix,lambda_xi_Det,lambda2_xi_invMatrix,lambda3_xi_invMatrix
	COMMON GEOMCOMMON, xiMin, xiMax
	COMMON DERIV,derivIndex

	;;; Initialise the normalised slope cumulant functions from the input parameters
	;;; and the SLOPECUMULANTS common block. They only need to be computed once for each time
	;;; secondGlintMomentFunction is called

	kappa_xi_11 = kappa_xi[0]
	kappa_xi_21 = kappa_xi[1]
	kappa_xi_12 = kappa_xi[2]

	lambda_xi_11 = kappa_xi_11/kappa_xi_2
	lambda_xi_21 = kappa_xi_21/(kappa_xi_2^(1.5D))
	lambda_xi_12 = kappa_xi_12/(kappa_xi_2^(1.5D))

	lambda_xi_Det = 1.D - lambda_xi_11^2.D

	lambda2_xi_invMatrix = 	[[1.D,-lambda_xi_11],[-lambda_xi_11,1.D]]/lambda_xi_Det
	
	d2r = !DPI/180.D
	r2d = 180.D/!DPI	

	M21_L = DBLARR(GEOM.N_angles)

	;;; Calculate the glint moment function M21_L for each geometry

	FOR geometry = 0L,N_ELEMENTS(GEOM.xi_min)-1L DO BEGIN

		AB_Limits = [GEOM.xi_min[geometry],GEOM.xi_max[geometry]]
		xiMin = AB_Limits[0]
		xiMax = AB_Limits[1]

		linPart = M2_L[geometry]
		nlinPart = (1.D/mu_L_1[geometry])*INT_2D('prob2D_21_xi_int',AB_Limits,'PQ_Limits',6,/DOUBLE)
		M21_L[geometry] = linPart + nlinPart
		;M21_L[geometry] = linPart
		;M21_L[geometry] = nlinPart

	ENDFOR
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO thirdGlintMomentFunction12,GEOM,kappa_xi,M2_L,M12_L

	;;; GEOM: Structure containing angular dependent parameters
	;;;		for a range of geometries.
	;;; kappa_xi: Vector containing the slope cumulant functions,
	;;;		kappa_xi_11, kappa_xi_21 and kappa_xi_12 for a single lag value.
	;;; M21_L: Vector containing the third glint moment function for each 
	;;;		geometry in GEOM, along the tau_1 axis, for a single lag value

	COMMON SLOPECUMULANTS,kappa_xi_1,kappa_xi_2,kappa_xi_3,kappa_xi_11,kappa_xi_21,kappa_xi_12, $
		kappa_xi_110,kappa_xi_101,kappa_xi_011, $
		kappa_xi_210,kappa_xi_120, $
		kappa_xi_201,kappa_xi_102, $
		kappa_xi_012,kappa_xi_021, $
		kappa_xi_111

	COMMON GLINTMOMENTS,mu_L_1,mu_L_2,mu_L_3,mu_L_11,mu_L_21,mu_L_12,mu_L_111
	COMMON LAMBDA,lambda2_xi_Matrix,lambda3_xi_Matrix,lambda_xi_Det,lambda2_xi_invMatrix,lambda3_xi_invMatrix
	COMMON GEOMCOMMON, xiMin, xiMax
	COMMON DERIV,derivIndex

	;;; Initialise the normalised slope cumulant functions from the input parameters
	;;; and the SLOPECUMULANTS common block. They only need to be computed once for each time
	;;; secondGlintMomentFunction is called

	kappa_xi_11 = kappa_xi[0]
	kappa_xi_21 = kappa_xi[1]
	kappa_xi_12 = kappa_xi[2]

	lambda_xi_11 = kappa_xi_11/kappa_xi_2
	lambda_xi_21 = kappa_xi_21/(kappa_xi_2^(1.5D))
	lambda_xi_12 = kappa_xi_12/(kappa_xi_2^(1.5D))

	lambda_xi_Det = 1.D - lambda_xi_11^2.D

	lambda2_xi_invMatrix = 	[[1.D,-lambda_xi_11],[-lambda_xi_11,1.D]]/lambda_xi_Det
	
	d2r = !DPI/180.D
	r2d = 180.D/!DPI	

	M12_L = DBLARR(GEOM.N_angles)

	;;; Calculate the glint moment function M21_L for each geometry
	FOR geometry = 0L,N_ELEMENTS(GEOM.xi_min)-1L DO BEGIN

		AB_Limits = [GEOM.xi_min[geometry],GEOM.xi_max[geometry]]
		xiMin = AB_Limits[0]
		xiMax = AB_Limits[1]

		;M12_L[geometry] = (1.D/mu_L_1[geometry])*(M2_L[geometry] + INT_2D('prob2D_12_xi_int',AB_Limits,'PQ_Limits',6,/DOUBLE))

		linPart = M2_L[geometry]
		nlinPart = (1.D/mu_L_1[geometry])*INT_2D('prob2D_12_xi_int',AB_Limits,'PQ_Limits',6,/DOUBLE)
		M12_L[geometry] = linPart + nlinPart
		;M12_L[geometry] = linPart
		;M12_L[geometry] = nlinPart
	ENDFOR
END
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; Function to calculate the third glint moment function for a
;;; particular value of the lag. This version sets the derivatives to unity,
;;; and the calling function computes the derivatives from differences

PRO thirdGlintMomentFunction2,GEOM,kappa_xi,M3_L

	;;; GEOM: Structure containing angular dependent parameters
	;;;		for a range of geometries.
	;;; kappa_xi: Vector containing the slope cumulant functions,
	;;;		kappa_xi_210, kappa_xi_120 and kappa_xi_111 for a 
	;;;     for a single lag value pair
	;;; M3_L: vector containing the third glint moment function for each
	;;;		geometry in GEOM, for a single lag value pair

	COMMON SLOPECUMULANTS,kappa_xi_1,kappa_xi_2,kappa_xi_3,kappa_xi_11,kappa_xi_21,kappa_xi_12, $
		kappa_xi_110,kappa_xi_101,kappa_xi_011, $
		kappa_xi_210,kappa_xi_120, $
		kappa_xi_201,kappa_xi_102, $
		kappa_xi_012,kappa_xi_021, $
		kappa_xi_111

	COMMON GLINTMOMENTS,mu_L_1,mu_L_2,mu_L_3,mu_L_11,mu_L_21,mu_L_12,mu_L_111
	COMMON LAMBDA,lambda2_xi_Matrix,lambda3_xi_Matrix,lambda_xi_Det,lambda2_xi_invMatrix,lambda3_xi_invMatrix
	COMMON GEOMCOMMON, xiMin, xiMax
	COMMON DERIV,derivIndex

	;;; Initialise the normalised slope cumulant function from the input parameters
	;;; and the SLOPECUMULANTS common block. They only need to be computed once for each time
	;;; secondGlintMomentFunction is called

	;;; These values are from the Levenberg-Marquardt routine, pass to the COMMON block
	;kappa_xi_210 = kappa_xi[0]  ;-- These are already set
	;kappa_xi_120 = kappa_xi[1]  ;-- These are already set
	kappa_xi_111 = kappa_xi[2]

	;;; Set the second order normalised cumulants, so we can calc the determinant and inverse 
	;;; covariance matrix...

	lambda_xi_110 = kappa_xi_110 / kappa_xi_2
	lambda_xi_101 = kappa_xi_101 / kappa_xi_2
	lambda_xi_011 = kappa_xi_011 / kappa_xi_2

	lambda_xi_Det = 1.D - lambda_xi_110^2.D - lambda_xi_101^2.D - lambda_xi_011^2.D $
		+ 2.D*lambda_xi_110*lambda_xi_101*lambda_xi_011

	a = lambda_xi_110
	b = lambda_xi_101
	c = lambda_xi_011

	lambda3_xi_invMatrix = [                                                $
							[ (1.D - c^2.D) , (b*c-a)     , (a*c-b)     ],  $
							[ (b*c-a)     , (1.D - b^2.D) , (a*b-c)     ],  $
							[ (a*c-b)     , (a*b-c)     , (1.D - a^2.D) ]   $
						] / lambda_xi_Det

	d2r = !DPI/180.D
	r2d = 180.D/!DPI	

	M3_L = DBLARR(GEOM.N_angles)

	;;; Calculate the glint moment function M3_L for each geometry
	FOR geometry = 0L,N_ELEMENTS(GEOM.xi_min)-1L DO BEGIN

		AB_Limits = [GEOM.xi_min[geometry],GEOM.xi_max[geometry]]
		xiMin = AB_Limits[0]
		xiMax = AB_Limits[1]

		M3_L[geometry] = $
			(1.D/mu_L_1[geometry])*INT_3D('prob3D_xi_int',AB_Limits,'PQ_Limits','UV_Limits',6,/DOUBLE)

	ENDFOR

END

