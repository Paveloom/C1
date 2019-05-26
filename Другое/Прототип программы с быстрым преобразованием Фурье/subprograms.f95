module subprograms
implicit none

     contains

     subroutine FFT(A,M,N)
     implicit none
     
     COMPLEX A(N), U, W, T
     
     N = 2**M
     NV2 = N/2
     NM1 = N - 1
     J ~ 1
     DO 7 I =	1, NM1
     IF(I .GE. J) GO TO 5
     T = A (J)
     A (J) « A (I) A(I) = T
     5	К = NV2
     6	IF(К .GE. J) GO TO 7
     J = J - К
     К = K/2
     GO TO 6
     7	J « J ♦ К
     PI = 3.141592653509793
     DO 20 L = 1,M
     LE = 2**L
     LE1 = LE/2
     U = (1.0,0.)
     W « CMPLX(COS(PI/LE1),SIN(PI/LE1) ) DO 20 J « 1,LE1
     DO 10 I = J,N,LE
     IP = I ♦ LE1
     T = A (IP) • U A(I₽» = A (I) - T 10 A (I) - A (I) ♦ T 20 U = U * W
     RETURN
     END

end
