
Program MultiN

  USE NumTypes
  USE Statistics

  Integer, Parameter :: N = 2, Ns = 20000000
  Real (kind=DP) :: X(N), Z(Ns, N), S(N,N), r, A(N,N)

  S(1,1) = 4.0_DP
  S(1,2) = 0.2_DP
  S(2,1) = 0.2_DP
  S(2,2) = 1.0_DP


  CALL MultiNormal(Z, (/1.0_DP, 2.0_DP/), S)

  r = 0.0_DP
  Do I = 1, Ns
     r = r + (Z(I,1)-1.0_DP)/2.0_DP * (Z(I,2)-2.0_DP)
  End Do
  r = r / Real(Ns-1,kind=DP)
  
  Write(*,*)"Pearson's correlation coefficient: ", r

  Stop
End Program MultiN
