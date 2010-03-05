;;; Functions to calculate some Hermite Polynomials

FUNCTION He,index, x
    IF (index EQ 1) THEN RETURN, x
    IF (index EQ 2) THEN RETURN, x*x - 1.D    
    IF (index EQ 3) THEN RETURN, x*x*x - 3.D*x    
END

FUNCTION chi,order,row,x

	COMMON SLOPECUMULANTS,kappa_1,kappa_2,kappa_3,kappa_11,kappa_21,kappa_12, $
		kappa_110,kappa_101,kappa_011, $
		kappa_210,kappa_120, $
		kappa_201,kappa_102, $
		kappa_012,kappa_021, $
		kappa_111

	COMMON LAMBDA,lambda2_Matrix,lambda3_Matrix,lambda_Det,lambda2_invMatrix,lambda3_invMatrix

	;PRINT,x
	;PRINT,lambda2_invMatrix[row-1L,*]

	CASE order OF
		2: BEGIN
			RETURN, TRANSPOSE(lambda2_invMatrix[row-1L,*]) ## x
		END
		3: BEGIN
			RETURN, TRANSPOSE(lambda3_invMatrix[row-1L,*]) ## x
		END
	ENDCASE
END

FUNCTION He2,index, x
	COMMON LAMBDA,lambda2_Matrix,lambda3_Matrix,lambda_Det,lambda2_invMatrix,lambda3_invMatrix

	chi1_x = chi(2,1,x)
	chi2_x = chi(2,2,x)

	;;; He_10
	IF (index[0] EQ 1 AND index[1] EQ 0) THEN RETURN, chi1_x

	;;; He_01
	IF (index[0] EQ 0 AND index[1] EQ 1) THEN RETURN, chi2_x

	;;; He_20
	IF (index[0] EQ 2 AND index[1] EQ 0) THEN RETURN,$
		chi1_x^2.D - lambda2_invMatrix[0,0]
	
	;;; He_02
	IF (index[0] EQ 0 AND index[1] EQ 2) THEN RETURN,$
		chi2_x^2.D - lambda2_invMatrix[1,1]

	;;; He_30
	IF (index[0] EQ 3 AND index[1] EQ 0) THEN RETURN,$
		chi1_x^3.D - 3.D*lambda2_invMatrix[0,0]*chi1_x
	
	;;; He_03
	IF (index[0] EQ 0 AND index[1] EQ 3) THEN RETURN,$
		chi2_x^3.D - 3.D*lambda2_invMatrix[1,1]*chi2_x

	;;; He_21
	IF (index[0] EQ 2 AND index[1] EQ 1) THEN RETURN,$
		(chi1_x^2.D)*chi2_x - lambda2_invMatrix[0,0]*chi2_x $
		- 2.D*lambda2_invMatrix[0,1]*chi1_x
	
	;;; He_12
	IF (index[0] EQ 1 AND index[1] EQ 2) THEN RETURN,$
		(chi2_x^2.D)*chi1_x - lambda2_invMatrix[1,1]*chi1_x $
		- 2.D*lambda2_invMatrix[1,0]*chi2_x
END

FUNCTION He3,index, x
	COMMON LAMBDA,lambda2_Matrix,lambda3_Matrix,lambda_Det,lambda2_invMatrix,lambda3_invMatrix

	chi1_x = chi(3,1,x)
	chi2_x = chi(3,2,x)
	chi3_x = chi(3,3,x)

	;;; He_300
	IF (index[0] EQ 3 AND index[1] EQ 0 AND index[2] EQ 0) THEN RETURN, $
		chi1_x^3.D - 3.D*lambda3_invMatrix[0,0]*chi1_x

	;;; He_030
	IF (index[0] EQ 0 AND index[1] EQ 3 AND index[2] EQ 0) THEN RETURN, $
		chi2_x^3.D - 3.D*lambda3_invMatrix[1,1]*chi2_x

	;;; He_003
	IF (index[0] EQ 0 AND index[1] EQ 0 AND index[2] EQ 3) THEN RETURN, $
		chi3_x^3.D - 3.D*lambda3_invMatrix[2,2]*chi3_x

	;;; He_210
	IF (index[0] EQ 2 AND index[1] EQ 1 AND index[2] EQ 0) THEN RETURN, $
		(chi1_x^2.D)*chi2_x - lambda3_invMatrix[0,0]*chi2_x - 2.D*lambda3_invMatrix[0,1]*chi1_x

	;;; He_120
	IF (index[0] EQ 1 AND index[1] EQ 2 AND index[2] EQ 0) THEN RETURN, $
		chi1_x*(chi2_x^2.D) - lambda3_invMatrix[1,1]*chi1_x - 2.D*lambda3_invMatrix[1,0]*chi2_x

	;;; He_201
	IF (index[0] EQ 2 AND index[1] EQ 0 AND index[2] EQ 1) THEN RETURN, $
		(chi1_x^2.D)*chi3_x - lambda3_invMatrix[0,0]*chi3_x - 2.D*lambda3_invMatrix[0,2]*chi1_x

	;;; He_102
	IF (index[0] EQ 1 AND index[1] EQ 0 AND index[2] EQ 2) THEN RETURN, $
		chi1_x*(chi3_x^2.D) - lambda3_invMatrix[2,2]*chi1_x - 2.D*lambda3_invMatrix[2,0]*chi3_x

	;;; He_012
	IF (index[0] EQ 0 AND index[1] EQ 1 AND index[2] EQ 2) THEN RETURN, $
		chi2_x*(chi3_x^2.D) - lambda3_invMatrix[2,2]*chi2_x - 2.D*lambda3_invMatrix[2,1]*chi3_x

	;;; He_021
	IF (index[0] EQ 0 AND index[1] EQ 2 AND index[2] EQ 1) THEN RETURN, $
		(chi2_x^2.D)*chi3_x - lambda3_invMatrix[1,1]*chi3_x - 2.D*lambda3_invMatrix[1,2]*chi2_x

	;;; He_111
	IF (index[0] EQ 1 AND index[1] EQ 1 AND index[2] EQ 1) THEN RETURN, $
		chi1_x*chi2_x*chi3_x - lambda3_invMatrix[0,1]*chi3_x      $
				             - lambda3_invMatrix[0,2]*chi2_x      $
						     - lambda3_invMatrix[1,2]*chi1_x
END
