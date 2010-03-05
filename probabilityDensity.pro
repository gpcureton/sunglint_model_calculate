;;; Function to compute 2D slope probability density function

FUNCTION prob2D_xi, xi1, xi2

	COMMON SLOPECUMULANTS,kappa_1,kappa_2,kappa_3,kappa_11,kappa_21,kappa_12, $
		kappa_110,kappa_101,kappa_011, $
		kappa_210,kappa_120, $
		kappa_201,kappa_102, $
		kappa_012,kappa_021, $
		kappa_111

	COMMON LAMBDA,lambda2_Matrix,lambda3_Matrix,lambda_Det,lambda2_invMatrix,lambda3_invMatrix

	nu1 = (xi1-kappa_1)/(SQRT(kappa_2))
	nu2 = (xi2-kappa_1)/(SQRT(kappa_2))

	nu = TRANSPOSE([nu1, nu2])

	;;; Compute the normalised cumulants
	lambda_3 = kappa_3/(kappa_2^1.5D)
	lambda_21 = kappa_21/(kappa_2^1.5D)
	lambda_12 = kappa_12/(kappa_2^1.5D)

	;;; Quadratic form nu^T . Inverse(Lambda) . nu

	quadForm = TRANSPOSE(nu) ## (lambda2_invMatrix ## nu)

	;;; Save Hermite polynomials
	He2_30 = He2([3,0],nu)
	He2_03 = He2([0,3],nu)
	He2_12 = He2([1,2],nu)
	He2_21 = He2([2,1],nu)

	;;; Compute the linear component of the pdf
	probDensLin = (1.D/(2.D*!DPI*SQRT(lambda_Det)*kappa_2))*EXP(-0.5D*quadForm)

	;;; Add the nonlinear correction
	probDens = probDensLin $
		*(1.D + (1.D/6.D)*(lambda_3*He2_30 + 3.D*lambda_21*He2_21 $
				 + 3.D*lambda_12*He2_12 + lambda_3*He2_03))

	RETURN,probDens

END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; Function to compute slope probability density function, for scalar xi1 
;;; and vector xi2. Since prob2D_xi() only accepts scalars, we must construct 
;;; a vector of values over xi2 for a single value of xi1.

FUNCTION prob2D_xi_int, xi1, xi2

	prob_xi_Return = DBLARR(N_ELEMENTS(xi2))
	FOR i=0,N_ELEMENTS(xi2)-1 DO prob_xi_Return[i]=prob2D_xi(xi1,xi2[i])
	RETURN,prob_xi_Return
END


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; Function to compute 3D slope probability density function

FUNCTION prob3D_xi, xi1, xi2, xi3

	COMMON SLOPECUMULANTS,kappa_1,kappa_2,kappa_3,kappa_11,kappa_21,kappa_12, $
		kappa_110,kappa_101,kappa_011, $
		kappa_210,kappa_120, $
		kappa_201,kappa_102, $
		kappa_012,kappa_021, $
		kappa_111

	COMMON LAMBDA,lambda2_Matrix,lambda3_Matrix,lambda_Det,lambda2_invMatrix,lambda3_invMatrix

	nu1 = (xi1-kappa_1)/(SQRT(kappa_2))
	nu2 = (xi2-kappa_1)/(SQRT(kappa_2))
	nu3 = (xi3-kappa_1)/(SQRT(kappa_2))

	nu = TRANSPOSE([nu1, nu2, nu3])

	;;; Quadratic form nu^T . Inverse(Lambda) . nu

	quadForm = TRANSPOSE(nu) ## (lambda3_invMatrix ## nu)

	;;; Compute the normalised cumulants
	lambda_300 = kappa_3 / (kappa_2^1.5D)
	lambda_030 = kappa_3 / (kappa_2^1.5D)
	lambda_003 = kappa_3 / (kappa_2^1.5D)

	lambda_210 = kappa_210  / (kappa_2^(1.5D))
	lambda_120 = kappa_120  / (kappa_2^(1.5D))

	lambda_201 = kappa_201  / (kappa_2^(1.5D))
	lambda_102 = kappa_102  / (kappa_2^(1.5D))

	lambda_021 = kappa_021  / (kappa_2^(1.5D))
	lambda_012 = kappa_012  / (kappa_2^(1.5D))

	lambda_111 = kappa_111 / (kappa_2^(1.5D))

	;;; Save Hermite polynomials
	He3_300 = He3([3,0,0],nu)
	He3_030 = He3([0,3,0],nu)
	He3_003 = He3([0,0,3],nu)

	He3_210 = He3([2,1,0],nu)
	He3_120 = He3([1,2,0],nu)
	He3_201 = He3([2,0,1],nu)
	He3_102 = He3([1,0,2],nu)
	He3_012 = He3([0,1,2],nu)
	He3_021 = He3([0,2,1],nu)

	He3_111 = He3([1,1,1],nu)

	;;; Compute the linear component of the pdf
	probDensLin = EXP(-0.5D*quadForm)/(((2.D*!DPI*kappa_2)^1.5D)*SQRT(lambda_Det))

	;;; Add the nonlinear correction
	probDens = probDensLin $
		*(1.D + (1.D/6.D)*(lambda_300*He3_300 + lambda_030*He3_030 + lambda_003*He3_003 + $ 
			3.D*(lambda_210 * He3_210  	     $
                + lambda_120 * He3_120 	     $
				+ lambda_201 * He3_201 	     $
				+ lambda_102 * He3_102 	     $
				+ lambda_012 * He3_012 	     $
				+ lambda_021 * He3_021)	     $
				+ 6.D*lambda_111 * He3_111   $
			))

	RETURN,probDens
	
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; Function to compute slope probability density function, for scalar xi1, xi2 
;;; and vector xi3. Since prob3D_xi() only accepts scalars, we must construct 
;;; a vector of values over xi3 for a single values of xi1 and xi2.

FUNCTION prob3D_xi_int, xi1, xi2, xi3
	prob_xi_Return = DBLARR(N_ELEMENTS(xi3))
	FOR i=0,N_ELEMENTS(xi3)-1 DO prob_xi_Return[i]=prob3D_xi(xi1,xi2,xi3[i])
	RETURN,prob_xi_Return
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Function to compute 2D slope probability density function for the case where
;;; lambda_101 goes to unity

FUNCTION prob2D_21_xi, xi1, xi2

	COMMON SLOPECUMULANTS,kappa_1,kappa_2,kappa_3,kappa_11,kappa_21,kappa_12, $
		kappa_110,kappa_101,kappa_011, $
		kappa_210,kappa_120, $
		kappa_201,kappa_102, $
		kappa_012,kappa_021, $
		kappa_111

	COMMON LAMBDA,lambda2_Matrix,lambda3_Matrix,lambda_Det,lambda2_invMatrix,lambda3_invMatrix

	nu1 = (xi1-kappa_1)/(SQRT(kappa_2))
	nu2 = (xi2-kappa_1)/(SQRT(kappa_2))

	nu = TRANSPOSE([nu1, nu2])

	;;; Compute the normalised cumulants
	lambda_3 = kappa_3/(kappa_2^1.5D)
	lambda_21 = kappa_21/(kappa_2^1.5D)
	lambda_12 = kappa_12/(kappa_2^1.5D)

	;;; Quadratic form nu^T . Inverse(Lambda) . nu

	quadForm = TRANSPOSE(nu) ## (lambda2_invMatrix ## nu)

	;;; Save Hermite polynomials
	He2_30 = He2([3,0],nu)
	He2_03 = He2([0,3],nu)
	He2_12 = He2([1,2],nu)
	He2_21 = He2([2,1],nu)

	;;; Compute the linear component of the pdf
	probDensLin = EXP(-0.5D*quadForm)/(2.D*!DPI*kappa_2*SQRT(lambda_Det))

	;;; Add the nonlinear correction
	probDens = probDensLin $
		*((1.D/6.D)*(  8.D*lambda_3*He2_30 $
		             + 1.D*lambda_3*He2_03 $
					 + 9.D*lambda_21*He2_21 $
					 + 3.D*lambda_12*He2_12 $
					 - 3.D*lambda_12*He2_21 $
					 - 3.D*lambda_21*He2_12 ))
					 ;+ 3.D*lambda_12*He2_21 $
					 ;+ 3.D*lambda_21*He2_12 ))

	RETURN,probDens

END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; Function to compute slope probability density function, for scalar xi1 
;;; and vector xi2. Since prob2D_xi() only accepts scalars, we must construct 
;;; a vector of values over xi2 for a single value of xi1.

FUNCTION prob2D_21_xi_int, xi1, xi2

	prob_xi_Return = DBLARR(N_ELEMENTS(xi2))
	FOR i=0,N_ELEMENTS(xi2)-1 DO prob_xi_Return[i]=prob2D_21_xi(xi1,xi2[i])
	RETURN,prob_xi_Return
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Function to compute 2D slope probability density function for the case where
;;; lambda_011 goes to unity

FUNCTION prob2D_12_xi, xi1, xi2

	COMMON SLOPECUMULANTS,kappa_1,kappa_2,kappa_3,kappa_11,kappa_21,kappa_12, $
		kappa_110,kappa_101,kappa_011, $
		kappa_210,kappa_120, $
		kappa_201,kappa_102, $
		kappa_012,kappa_021, $
		kappa_111

	COMMON LAMBDA,lambda2_Matrix,lambda3_Matrix,lambda_Det,lambda2_invMatrix,lambda3_invMatrix

	nu1 = (xi1-kappa_1)/(SQRT(kappa_2))
	nu2 = (xi2-kappa_1)/(SQRT(kappa_2))

	nu = TRANSPOSE([nu1, nu2])

	;;; Compute the normalised cumulants
	lambda_3 = kappa_3/(kappa_2^1.5D)
	lambda_21 = kappa_21/(kappa_2^1.5D)
	lambda_12 = kappa_12/(kappa_2^1.5D)

	;;; Quadratic form nu^T . Inverse(Lambda) . nu

	quadForm = TRANSPOSE(nu) ## (lambda2_invMatrix ## nu)

	;;; Save Hermite polynomials
	He2_30 = He2([3,0],nu)
	He2_03 = He2([0,3],nu)
	He2_12 = He2([1,2],nu)
	He2_21 = He2([2,1],nu)

	;;; Compute the linear component of the pdf
	probDensLin = EXP(-0.5D*quadForm)/(2.D*!DPI*kappa_2*SQRT(lambda_Det))

	;;; Add the nonlinear correction
	probDens = probDensLin $
		*((1.D/6.D)*(  1.D*lambda_3*He2_30 $
					 + 8.D*lambda_3*He2_03 $
					 + 6.D*lambda_21*He2_21 $
					 +12.D*lambda_12*He2_12 ))

	RETURN,probDens

END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; Function to compute slope probability density function, for scalar xi1 
;;; and vector xi2. Since prob2D_xi() only accepts scalars, we must construct 
;;; a vector of values over xi2 for a single value of xi1.

FUNCTION prob2D_12_xi_int, xi1, xi2

	prob_xi_Return = DBLARR(N_ELEMENTS(xi2))
	FOR i=0,N_ELEMENTS(xi2)-1 DO prob_xi_Return[i]=prob2D_12_xi(xi1,xi2[i])
	RETURN,prob_xi_Return
END
