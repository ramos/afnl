
Program GEVPprg

  USE NumTypes
  USE Linear
  USE LapackAPI

  Integer :: N = 2, Is
  Complex (kind=DPC), Allocatable :: A(:,:), B(:,:), v(:,:)
  Real (kind=DP), Allocatable :: l(:), R(:)

  ALLOCATE(A(N,N), B(N,N), v(N,N), l(N), R(N))

  A(1,1) = 1.0
  A(1,2) = 2.0
  A(2,1) = 2.0
  A(2,2) = 3.0

  B(1,1) = 2.0
  B(1,2) = 3.0
  B(2,1) = 3.0
  B(2,2) = 8.0

  Write(*,*)'Det(B)', Det(B)

  CALL GEVP(A, B, l, v, Is)
  Write(*,*)'Success: ', Is
  Write(*,*)l
  Write(*,*)
  Do I = 1, 2
     Do J = 1, 2
        Write(*,*)I,J, v(I,J)
     End Do
  End Do
  Write(*,*)

  Write(*,*)'Residuals: '
  Do J = 1, N
     Do I = 1, N
        R(I) = Abs(Sum( (A(I,:)-l(J)*B(I,:))*v(:,J)))
     End Do
     Write(*,*)R
  End Do


  Stop
End Program GEVPprg
