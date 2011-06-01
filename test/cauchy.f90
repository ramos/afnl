
Program CacuchyTest

  USE NumTypes
  USE Statistics
  USE NonNumeric

  Integer, Parameter :: N = 5000000, Nd = 200
  Real (kind=DP) :: X(N), Y(N), V1, V2, h, Div(Nd), alph, b, s, m
  Integer :: Nt1(Nd), Nt2(Nd), Id
  

  Read(*,*)b, s, m
 
  CALL Cauchy(X, b, s, m)
  CALL Normal(Y, 0.0_DP, SR2_DP)

  h = 20.0_DP/Real(Nd,kind=DP)
  Do I = 1, Nd
     Div(I) = -10.0_DP + (I-1)*h
  End Do

  Nt1 = 0
  Nt2 = 0
  Do I = 1, N
     Id = Locate(Div, X(I))
     Nt1(Id) = Nt1(Id) +1
     Id = Locate(Div, Y(I))
     Nt2(Id) = Nt2(Id) +1
  End Do

  Do I = 2, Nd-1
     Write(*,*)-10.0_DP + (I-1)*h, Real(Nt1(I))/Real(N), Real(Nt2(I))/Real(N)
  End Do

  Write(0,*)Mean(X), Stddev(X)
  Write(0,*)Mean(Y), Stddev(Y)


  Stop
End Program CacuchyTest
