
Program OPT

  USE NumTypes
  USE NonNumeric

  Integer :: NN
  Character (len=10) :: aa
  Real (kind=SP) :: r
  Real (kind=DP) :: d


  If (GetOpt('-i', NN)) Then
     Write(*,*)'Valor del argumento -i: ', NN
  Else
     Write(*,*)'Argumento -i no esta presente'
  End If

  If (GetOpt('-c', aa)) Then
     Write(*,*)'Valor del argumento -c: ', aa
  Else
     Write(*,*)'Argumento -c no esta presente'
  End If

  If (GetOpt('-s', r)) Then
     Write(*,*)'Valor del argumento -s: ', r
  Else
     Write(*,*)'Argumento -s no esta presente'
  End If

  If (GetOpt('-d', d)) Then
     Write(*,*)'Valor del argumento -d: ', d
  Else
     Write(*,*)'Argumento -d no esta presente'
  End If

  If (GetOpt('-h')) Then
     Write(*,*)' Usage: opt.x [options]'
     Write(*,*)'   -i int'
     Write(*,*)'   -c character string'
     Write(*,*)'   -s real'
     Write(*,*)'   -d double'
     Write(*,*)'   -h'
  End If

  Stop
End Program OPT
