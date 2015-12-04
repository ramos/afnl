
Program AA

  USE NumTypes
  USE NonNumeric
  USE ModBDIO
!  USE Statistics

  USE Random

  Integer, Parameter :: Ns=5

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

  CALL rnlx_start()
  CALL rndm_seed((/1,2,3,4/))
  
  call BDIO_setuser('perico')
  call BDIO_sethost('casa.pepe')

  bb = BDIO_open(fn, 'w', 'Test file with RANLUX random numbers')
  Do i4 = 1, 20
     CALL BDIO_start_record(bb,BDIO_BIN_F64LE, 2, .True.)
     Do Isw = 1, 50
        CALL rndm(da)
        j = BDIO_write(bb, da,.true.)
     End Do
     CALL BDIO_write_hash(bb)
  End Do

  call BDIO_close(bb)

  bb = BDIO_open(fn, 'a')
  CALL BDIO_start_record(bb,BDIO_BIN_F64LE, 2, .True.)
  Do Isw = 1, 5000
     CALL rndm(da)
     j = BDIO_write(bb, da,.true.)
  End Do
  CALL BDIO_write_hash(bb)


  call BDIO_close(bb)
  
  bb = BDIO_open(fn, 'r')
  do
     if (BDIO_get_uinfo(bb)==2) then
        j = BDIO_read(bb, da)
        write(*,*)j, da
        stop
     end if
     if (.not.BDIO_seek(bb)) stop
  end do
  
  Stop
End Program AA
