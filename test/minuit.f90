Program TestAPI

  USE MinuitAPI
  USE Statistics
  USE Constants

  Integer, Parameter :: N = 3
  Real (kind=8) :: X(N), Y(N), Ye(N), C(2), Ch

  Interface 
     Function Func(X)
       Real (kind=8), Intent (in) :: X(:)
       Real (kind=8) :: Func
     End Function Func
  End Interface
  

  X(:) = 20.0D0
  CALL Minimize(Func, X, Ch)
  Write(*,*)'= MINIMIZE ='
  Write(*,*)'Function value:  ', Ch
  Write(*,*)'Point of minima: ', X
  Write(*,*)

  X(:) = -20.0D0
  CALL Minimize(Func, X, Ch, (/1/))
  Write(*,*)'= MINIMIZE ONLY VAR 1 ='
  Write(*,*)'Function value:  ', Ch
  Write(*,*)'Point of minima: ', X
  Write(*,*)

  X(:) = -20.0D0
  CALL Minimize(Func, X, Ch, (/2/))
  Write(*,*)'= MINIMIZE ONLY VAR 2 ='
  Write(*,*)'Function value:  ', Ch
  Write(*,*)'Point of minima: ', X
  Write(*,*)

  X(:) = -20.0D0
  CALL Minimize(Func, X, Ch, (/1:3/))
  Write(*,*)'= MINIMIZE ONLY VAR 1,3 ='
  Write(*,*)'Function value:  ', Ch
  Write(*,*)'Point of minima: ', X
  Write(*,*)

  Stop
End Program TestAPI

Function Func(X)

  USE NumTypes

  Real (kind=8), Intent (in) :: X(:)
  Real (kind=8) :: Func

  Func = (X(1)-1.0_DP)**2 + (X(2) - 2.0_DP)**8 + (X(3) - 3.0_DP)**4

  Return
End Function Func
