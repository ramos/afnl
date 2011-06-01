
Program Test

  USE NumTypes
  USE Linear

  Integer, Parameter :: N = 47
  Real (kind=DP) :: M(N,N), L(N,N)

  CALL Random_Number(M)
  
  ! Try tomake it symmetric and positive definite
  M = M + Transpose(M)
  Do I = 1, N
     M(I,I) = M(I,I) + Real(N,kind=DP)
  End Do


  ! Compute the Cholesky decomposition and check it is ok
  L = Cholesky(M)
  Write(*,*)'Should be zero: ', Sum(Abs(M-MatMul(L,Transpose(L))))


  Stop
End Program Test
