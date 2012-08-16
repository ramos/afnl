
Program Checksum

  USE NonNumeric
  use iso_c_binding, only: c_ptr

  Type data
     Real (kind=DP) :: C(100)
     Integer :: N
     Character (len=100) :: Comment
  End type data

  Type (data) :: ex, ex2
  Type (c_ptr) :: p
  Integer, Allocatable :: buf(:), dt(:)
  Integer :: I
  Real (kind=DP), Allocatable :: X(:)
  Logical :: Ierr

  Allocate(X(100))
  CALL Random_Number(X)
  Allocate(buf(Size(Transfer(X, buf))))
  buf = Transfer(X, buf)
  I = Hash(buf)
  Write(*,'(1I8.8,2X,1Z9.9)')Size(buf), I

  CALL WriteBuffer(buf, 'data', 'Estos datos no valen ni pa cagarla!')
  CALL ReadBuffer(dt,'data', Ierr)

  Write(*,*)Ierr

  Stop
End Program Checksum
