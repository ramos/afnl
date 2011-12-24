
Program CacuchyTest

  USE NumTypes
  USE Statistics
  USE NonNumeric

  Integer, Parameter :: N = 50000, Nd = 200
  Real (kind=DP) :: X(N), Y(N), Z(N), V1, V2, h, Div(Nd), alph, b, s, m
  Integer :: Nt1(Nd), Nt2(Nd), Nt3(Nd), Id
  

  Read(*,*)b, s, m
 
  CALL Cauchy(X, b, s, m)
  CALL Normal(Y, 0.0_DP, 1.0_DP)

  Z = 2*X/(s**2+X**2)

  h = 20.0_DP/Real(Nd,kind=DP)
  Do I = 1, Nd
     Div(I) = -10.0_DP + (I-1)*h
  End Do

  Nt1 = 0
  Nt2 = 0
  Nt3 = 0
  Do I = 1, N
     Id = Locate(Div, X(I))
     Nt1(Id) = Nt1(Id) +1
     Id = Locate(Div, Y(I))
     Nt2(Id) = Nt2(Id) +1
     Id = Locate(Div, Z(I))
     Nt3(Id) = Nt3(Id) +1
  End Do

  Do I = 2, Nd-1
     Write(*,'(100ES23.15)')-10.0_DP + I*h, &
          & Real(Nt1(I))/Real(N)/h, &
          & Real(Nt2(I))/Real(N)/h, &
          & Real(Nt3(I))/Real(N)/h
  End Do



  Write(0,*)Mean(X), Stddev(X)
  Write(0,*)Mean(Y), Stddev(Y)
  Write(0,*)Sum(Nt3(:)), N

  Stop
End Program CacuchyTest
