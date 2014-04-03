
Program AA

  USE NumTypes
  USE NonNumeric
  USE ModBDIO
  USE Statistics

  Integer, Parameter :: Ns=58754

  Type (BDIO) :: bb
  Character (len=200) :: fn, type
  Character (kind=1) :: ch
  Integer :: Isw, I(100)
  Integer (kind=8) :: di(3) = (/-1,-2,-3/)
  Real (kind=4) :: dd(10)
  Integer, Allocatable :: ibuf(:)
  Real (kind=DP), Allocatable :: buf(:)
  Integer :: Ia = 805306416!805306368
  Real (kind=4) :: r4
  Real (kind=8) :: r8
  Complex (kind=DPC) :: z16

  Real (kind=DP) :: da(Ns), db(Ns)
  Complex (kind=DPC) :: zarr(Ns)


  If (.not.GetOpt('-i', fn)) Stop


  CALL SetLuxLevel(5)
  CALL Putluxseed()
  I(1) = LUXrandomize()
  CALL Random_Number(dd)



  Forall (j=1:100) I(j) = j


  bb = BDIO_open(fn, 'a')!, 'Test file with complex numbers')
!  CALL BDIO_show(bb)   alberto
  CALL BDIO_start_record(bb,BDIO_BIN_F64LE,2,.true.)
  CALL Normal(db)
  CALL Normal(da)
  zarr = Cmplx(da, db)
  Write(*,*)BDIO_write(bb,zarr), ' Complex written'

  CALL BDIO_start_record(bb,BDIO_BIN_INT32LE,2,.false.)
  Write(*,*)BDIO_write(bb,I)
  CALL BDIO_start_record(bb,BDIO_BIN_INT32LE,2,.true.)
  Write(*,*)BDIO_write(bb,I)

  Stop
  zarr = Cmplx(0.0_DP,0.0_DP)
  bb = BDIO_Open(fn,'r')
!  Write(*,*)bb%tcnt, bb%rcnt, bb%hcnt
  CALL BDIO_Seek(bb,bb%tcnt-1)
  Write(*,*)BDIO_read(bb,zarr(1:2)), ' Complex read!'
!  Write(*,*)zarr(1:4)


  Stop
!!$  CALL BDIO_start_record(bb,BDIO_BIN_INT32LE,7)
!!$  CALL BDIO_write(bb,I(:))
!!$  CALL BDIO_write(bb,I(:))
!!$  CALL BDIO_write(bb,I(:))
!!$  CALL BDIO_write(bb,I(:))
!!$  CALL BDIO_write(bb,I(2:3))
!!$
!!$  Forall (j=1:100) I(j) = j
!!$  CALL BDIO_start_record(bb,BDIO_BIN_INT32BE,7)
!!$  CALL BDIO_write(bb,I(:))
!!$  CALL BDIO_start_record(bb,BDIO_BIN_INT64LE,4)
!!$  CALL BDIO_write(bb,di(:))
!!$
!!$  CALL BDIO_start_record(bb,BDIO_BIN_F32LE,4)
!!$  CALL BDIO_write(bb,dd(:))
!!$
!!$
!!$


  CALL BDIO_show(bb)
  Allocate(ibuf(134), buf(343))
  If (GetOpt('-n', isw)) Then
     CALL BDIO_seek(bb, isw)
     Write(*,*)bb%current%rfmt, BDIO_BIN_F64LE
     if (bb%current%rfmt == BDIO_BIN_F64LE) Then
        Write(*,*)'INT'
        Write(*,*)'Leido: ', BDIO_read(bb, buf(1:15),.True.)
        Write(*,*)'Leido: ', BDIO_read(bb, buf(16:),.True.)
!        Write(*,'(100I19)')buf(1:)
        Write(*,'(100F19.8)')buf(1:)
        Write(*,*)hash(buf(:100)), bb%current%hash
     Else
!        CALL BDIO_read(bb, buf)
        Write(*,'(100F14.10)')buf(1:200)
     End if
  End If
  
      


  Stop
End Program AA
