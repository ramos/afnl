
Program FFTPR


  USE NumTypes
  USE Statistics
  USE Fourier
  USE TIME

  Integer, Parameter :: N = 4

  Type (Fourier_Serie) :: S1, S2, S3, S4
  Complex (kind=DPC) :: CD(N)
  Complex (kind=DPC), Allocatable :: CCC(:)
  Real (kind=DP) :: dr(N), di(N)
  Integer, Allocatable :: Ipt(:)

  CALL Random_Number(dr)
  CALL Random_Number(di)

  CD(:) = Cmplx(dr, di)


  CALL Init_Serie(S1,N/2)
  CALL Init_Serie(S2,N/2)
  CALL Init_Serie(S3,N/2)
  CALL Init_Serie(S4,N/2)
  
  Allocate(Ipt(-N/2:N/2))
  
  ForAll (I=1:N/2) Ipt(I) = I
  Ipt(-N/2) = N/2
  Ipt(0)    = N
  ForAll (I=1:N/2-1) Ipt(-N/2+I) = N/2+I
  
  Do I = -N/2, N/2
     S1%Coef(I) = Cd(Ipt(I))
  End Do

  Write(*,*)'GENERATION DONE'
  Write(*,*)asctime(gettime())
  S2 = FFT(CD)
  Write(*,*)'FFT'
  Write(*,*)asctime(gettime())
  S3 = FFT(S1)
  Write(*,*)'FASTFT'
  Write(*,*)asctime(gettime())
  S4 = DFT(CD)
  Write(*,*)'DFT'
  Write(*,*)asctime(gettime())
  Write(*,*)'Result: ', Sum(Abs(S2%Coef-S3%Coef)), MaxVal(Abs(S2%Coef-S3%Coef))
  Write(*,*)'Result: ', Sum(Abs(S2%Coef-S4%Coef)), MaxVal(Abs(S2%Coef-S4%Coef))

  Write(*,*)
  Write(*,*)"********* INVERSE TRANSFORMATION ************"
  Write(*,*)
  S3 = FFT(S2,1)
  Write(*,*)'Result: ', Sum(Abs(S1%Coef-S3%Coef)), MaxVal(Abs(S1%Coef-S3%Coef))  

  CALL Fourier2Data(S1, CCC)
  Write(*,*)'DONE', CCC(:), Size(CCC)
  
  CALL Data2Fourier(CCC, S3)
  Write(*,*)'Result transforms: ', Sum(Abs(S1%Coef-S3%Coef)), MaxVal(Abs(S1%Coef-S3%Coef))  


10 FORMAT(100ES13.5)

  Stop
End Program FFTPR
