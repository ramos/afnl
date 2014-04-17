
Program AA

  USE NumTypes
  USE NonNumeric
  USE ModBDIO
  USE Statistics

  Integer, Parameter :: Ns=5045507

  Type (BDIO) :: bb
  Character (len=200) :: fn, type
  Character (kind=1) :: ch
  Integer :: Isw, I(100)
  Integer (kind=8) :: di(3) = (/-1,-2,-3/)
  Real (kind=4) :: dd(10)
  Integer, Allocatable :: ibuf(:)
  Real (kind=DP), Allocatable :: buf(:)
  Integer :: Ia = 805306416, i4
  Real (kind=4) :: r4
  Real (kind=8) :: r8
  Complex (kind=DPC) :: z16

  Real (kind=DP) :: da(Ns), db(Ns)
  Complex (kind=DPC) :: zarr(Ns)
  Character (len=1) :: c(4)

  If (.not.GetOpt('-i', fn)) Stop


  CALL SetLuxLevel(5)
  CALL Putluxseed()
  I(1) = LUXrandomize()
  CALL Random_Number(dd)



  Forall (j=1:100) I(j) = j


  bb = BDIO_open(fn, 'w')!, 'Test file with complex numbers')
  CALL Normal(da)
  CALL Normal(db)
  zarr = Cmplx(da,db)

  CALL BDIO_start_record(bb,BDIO_BIN_F64LE, 2, .True.)
  Write(*,*)'Writting complex numbers'
  j = BDIO_write(bb, da, .true.)
  Write(*,*)'Written ', j, ' complex numbers'
  CALL BDIO_write_hash(bb)

  Stop
  CALL BDIO_start_record(bb,BDIO_BIN_F64LE, 2, .True.)
  j = BDIO_write(bb, zarr, .true.)
  CALL BDIO_write_hash(bb)



  Stop
End Program AA
