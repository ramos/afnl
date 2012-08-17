
Program Checksum

  USE NonNumeric
  use iso_c_binding, only: c_ptr

  Type data
     Real (kind=DP) :: C(200)
     Integer :: N
     Character (len=100) :: Comment
  End type data

  Type (data) :: ex(2)
  Integer, Allocatable :: buf(:)
  Integer :: I, J

  Allocate(buf(Size(Transfer(ex, buf))))
  buf = Transfer(ex, buf)

  Do I = 1, 10
     J = Hash(buf)
  End Do

  Write(*,'(1I8.8,2X,1Z9.9)')Size(buf), I

  CALL WriteBuffer(buf, 'data', 'Estos datos no valen ni pa cagarla!')

  Stop
End Program Checksum
