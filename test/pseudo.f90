Program SVDtest

  USE NumTypes
  USE LapackAPI

  Real (kind=DP), Allocatable :: A(:,:), B(:,:)
  Character (len=100) :: foo
  

  Write(*,*)'Enter matrix size: '
  Read(*,*)N

  ALLOCATE(A(N,N), B(N,N))
  CALL Random_Number(A)
  Write(*,*)'A = ['
  Do K1 = 1, N
     Write(foo,*)'(1X, ', N,'ES13.5,1A1)'
     Write(*,foo)(A(K1,K2), K2=1, N), ';'
  End Do
  Write(*,*)']'

  B = PseudoInverse(A)

  Write(*,*)'B = ['
  Do K1 = 1, N
     Write(foo,*)'(1X, ', N,'ES13.5,1A1)'
     Write(*,foo)(B(K1,K2), K2=1, N), ';'
  End Do
  Write(*,*)']'

  A = MatMul(A,B)
  Write(*,*)'A = ['
  Do K1 = 1, N
     Write(foo,*)'(1X, ', N,'ES13.5,1A1)'
     Write(*,foo)(A(K1,K2), K2=1, N), ';'
  End Do
  Write(*,*)']'

  Stop
End Program SVDtest
