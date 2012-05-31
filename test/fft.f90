
Program FFT


  USE NumTypes
  USE Statistics
  USE Fourier

  Integer, Parameter :: N = 4

  Type (Fourier_Serie) :: S1, S2, S3
  Complex (kind=DPC) :: CD(N)
  Real (kind=DP) :: dr(N), di(N)
  
  CALL Random_Number(dr)
  CALL Random_Number(di)

  CD(:) = Cmplx(dr, di)


  Write(*,10)CD(:)
  CALL Init_Serie(S1,N/2)
  CALL Init_Serie(S2,N/2)
  CALL Init_Serie(S3,N/2)
  

10 FORMAT(100ES13.5)

  Stop
End Program FFT
