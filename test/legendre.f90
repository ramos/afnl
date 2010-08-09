
Program TestLeg


  USE NumTypes
  USE SpecialFunc

  Real (kind=DP) :: Y

  Read(*,*)L

  Do I = 0, L
     Y = Legendre(L, I, 0.50_DP)
     Write(*,*)I, Y, Legendre(L, I, 0.5_SP)
  End Do

  Stop
End Program TestLeg
