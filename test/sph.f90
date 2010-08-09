
Program TestLeg


  USE NumTypes
  USE SpecialFunc

  Complex (kind=DPC) :: Y

  Read(*,*)L

  Do I = 0, L
     Y = SphericalHarmonic(L, I, 0.50_DP, 0.35_DP)
     Write(*,'(1A,1I4,1A,1I4)')'Spherical Harmonic l=', L, ' m=+-', I
     Write(*,*)'  ', Y
     Write(*,*)'  ', SphericalHarmonic(L, -I, 0.5_DP, 0.35_DP)
  End Do

  Stop
End Program TestLeg
