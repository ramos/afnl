
Program AA

  USE NumTypes
  USE NonNumeric
  USE BDIO

  Type (BDIO_t) :: bb
  Character (len=200) :: fn, type
  Character (kind=1) :: ch
  Integer :: I(1), Isw
  Integer :: Ia = 805306416!805306368
  Real (kind=4) :: r4
  Real (kind=8) :: r8
  Complex (kind=DPC) :: z16


  If (.not.GetOpt('-i', fn)) Stop

  bb = BDIO_open(fn, 'r')

  Do 
     CALL BDIO_seek_record(bb)
     If (bb%Istate == BDIO_N_STATE) Exit
     Write(*,*)bb%ridx, bb%rfmt, bb%ruinfo, bb%rlen
  End Do


  
  Stop
End Program AA
