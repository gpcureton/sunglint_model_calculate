;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Function to compute derivatives of slope probability density function, and
;;; the derivative specified by derivIndex

FUNCTION prob_xi_deriv, xi1, xi2

	COMMON SLOPECUMULANTS,kappa_1,kappa_2,kappa_3,kappa_11,kappa_21,kappa_12
	COMMON LAMBDA,lambda_Matrix,lambda_Det,lambda_invMatrix
	COMMON DERIV,derivIndex

	nu1 = (xi1-kappa_1)/(SQRT(kappa_2))
	nu2 = (xi2-kappa_1)/(SQRT(kappa_2))

	nu = [[nu1], [nu2]]
	
	lambda_3 = kappa_3/(kappa_2^1.5D)
	lambda_11 = kappa_11/kappa_2
	lambda_21 = kappa_21/(kappa_2^1.5D)
	lambda_12 = kappa_21/(kappa_2^1.5D)

	
	;;; Rows of inverse Lambda matrix
	b1 = lambda_invMatrix[*,0]
	b2 = lambda_invMatrix[*,1]
	;;; Column "unit" vector
	e = [[1],[1]]

	;;; Quadratic form nu^T . Inverse(Lambda) . nu

	quadForm = TRANSPOSE(nu) ## (lambda_InvMatrix ## nu)

	;;; Save Hermite polynomials

	He2_30 = He2([3,0],nu)
	He2_03 = He2([0,3],nu)
	He2_21 = He2([2,1],nu)
	He2_12 = He2([1,2],nu)
	He2_20 = He2([2,0],nu)
	He2_02 = He2([0,2],nu)
	He2_11 = He2([1,1],nu)
	He2_10 = He2([1,0],nu)
	He2_01 = He2([0,1],nu)

	probDensLin = (1.D/(2.D*!DPI*SQRT(lambda_Det)*kappa_2))*EXP(-0.5D*quadForm)

	probDens = probDensLin $
		*(1.D + (1.D/6.D)*(lambda_3*He2_30 + 3.D*lambda_21*He2_21 $
				 + 3.D*lambda_21*He2_12 + lambda_3*He2_03))

	CASE derivIndex OF

; 		1: BEGIN
; 			probDens_dkappa1 = (kappa_2^(-0.5D)) $
; 				*(TRANSPOSE(e) ## (lambda_InvMatrix ## nu))*probDens $
; 				+ probDensLin * (1.D/6.D) * (-3.D*lambda_3/(kappa_2^0.5D)) $
; 				* ((e # b1)*He2_20 $
; 				+ (e # b2)*He2_20 + 2.D*(e # b1)*He2_11 $
; 				+ (e # b1)*He2_02 + 2.D*(e # b2)*He2_11 $
; 				+ (e # b2)*He2_02)
; 
; 			RETURN,probDens_dkappa1
; 		END
; 
; 		2: BEGIN
; 			F_30 = He2_10*(He2_20-(2.D*lambda_invMatrix[0,1])^2.D)
; 
; 			F_21 = He2_01*(He2_02-(2.D*lambda_invMatrix[1,0])^2.D) $
; 				+ 2.D*He2_10*(He2_11 $
; 				- lambda_11*(2.D*lambda_invMatrix[0,1])^2.D - 2.D*lambda_invMatrix[0,1])
; 
; 			F_12 = He2_10*(He2_20-(2.D*lambda_invMatrix[0,1])^2.D) $
; 				+ 2.D*He2_01*(He2_11 $
; 				- lambda_11*(2.D*lambda_invMatrix[1,0])^2.D - 2.D*lambda_invMatrix[1,0])
; 
; 			F_03 = He2_01*(He2_02-(2.D*lambda_invMatrix[1,0])^2.D)
; 
; 			probDens_dkappa2 = -((1.D - lambda_11*lambda_invMatrix[0,1])*(1.D - quadForm) $
; 				+ 0.5D*(transpose(nu) ## nu)/lambda_Det)*probDens/kappa_2 $
; 				+ 3.D*probDens/(2.D*kappa_2) $
; 				+ 3.D*probDens/(2.D*kappa_2)*probDensLin $
; 				* (1.D - (1.D/6.D)*(lambda_3*F_30 + lambda_21*F_21 $
; 				+ lambda_12*F_12 + lambda_3*F_03))
; 
; 			RETURN,probDens_dkappa2
; 		END
; 
; 		3: BEGIN
; 			probDens_dkappa3 = probDensLin*(1.D/6.D)*(He2_30 - He2_03)/(kappa_2^1.5D)
; 			RETURN,probDens_dkappa3
; 		END

		11: BEGIN
			chi1 = He2_10
			chi2 = He2_01
			chi1_dkappa11 = (2.D*lambda_11*chi1 - nu2)/(kappa_2*lambda_Det)
			chi2_dkappa11 = (2.D*lambda_11*chi2 - nu1)/(kappa_2*lambda_Det)
			I_rot = [[0.D, 1.D],[1.D, 0.D]]
			lambda_invMatrix_dkappa11 = -(I_rot - 2.D*lambda_11*lambda_InvMatrix)/(kappa_2*lambda_Det)
			
			He2_30_dkappa11 = 3.D*He2_20*chi1_dkappa11 - 3.D*lambda_invMatrix_dkappa11[0,0]*He2_10

			He2_03_dkappa11 = 3.D*He2_02*chi2_dkappa11 - 3.D*lambda_invMatrix_dkappa11[1,1]*He2_01

			He2_21_dkappa11 = 2.D*chi1_dkappa11*He2_11 + chi2_dkappa11*He2_20 $
				- lambda_invMatrix_dkappa11[0,0]*He2_01 - 2.D*lambda_invMatrix_dkappa11[0,1]*He2_10

			He2_12_dkappa11 = 2.D*chi2_dkappa11*He2_11 + chi1_dkappa11*He2_02 $
				- lambda_invMatrix_dkappa11[1,1]*He2_10 - 2.D*lambda_invMatrix_dkappa11[1,0]*He2_01

			probDens_dkappa11 = (nu1*nu2 + lambda_11*(1.D - quadForm))*probDens $
				+ probDensLin * (lambda_3*He2_30_dkappa11 + 3.D * lambda_21*He2_21_dkappa11 $
				+ 3.D * lambda_12 * He2_12_dkappa11 + lambda_3 * He2_03_dkappa11)/6.D
			
			RETURN,probDens_dkappa11
		END

		21: BEGIN
			probDens_dkappa21 = 0.5D * probDensLin * (He2_21 + He2_12) / (kappa_2^1.5D)
			RETURN,probDens_dkappa21
		END
	ENDCASE

END
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;; Function to compute derivatives of slope probability density function, xi1 is a scalar,
;;; while xi2 is a vector. Must return a vector of the same length as xi2

FUNCTION prob_xi_deriv_int, xi1, xi2

	COMMON DERIV,derivIndex

	prob_xi_deriv_Return = DBLARR(N_ELEMENTS(xi2))
	FOR i=0,N_ELEMENTS(xi2)-1 DO prob_xi_deriv_Return[i]=prob_xi_deriv(xi1,xi2[i])
	RETURN,prob_xi_deriv_Return

END
