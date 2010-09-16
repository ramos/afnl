
Program Test

  USE NumTypes
  USE SpecialFunc

  Real (kind=DP) :: X = 0.6783457843_DP, Y

  Write(*,'(ES33.25)')inverf(-0.998_DP)
  Write(*,'(ES33.25)')inverf(-0.5_DP)
  Write(*,'(ES33.25)')inverf(0.988_DP)

  Y = inverf(X)
  Write(*,*)'Same number two times!'
  Write(*,'(ES33.25)')X
  Write(*,'(ES33.25)')erf(Y)


  Stop
End Program Test
