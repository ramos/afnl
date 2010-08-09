Program TestAPI

  USE MinuitAPI
  USE Statistics
  USE Constants

  Integer, Parameter :: N = 2
  Real (kind=8) :: X(N), Y(N), Ye(N), C(2), Ch

  Interface 
     Function Func(X)
       Real (kind=8), Intent (in) :: X(:)
       Real (kind=8) :: Func
     End Function Func
  End Interface
  

  X(:) = -10.0D0
  CALL Minimize(Func, X, Ch)
  Write(*,*)'= MINIMIZE ='
  Write(*,*)'Function value:  ', Ch
  Write(*,*)'Point of minima: ', X
  Write(*,*)'Check:           ', Tan(X(1)), Cos(X(2)) 

  X(:) = -20.0D0
  CALL Migrad(Func, X, Ch)
  Write(*,*)'= MIGRAD ='
  Write(*,*)'Function value:  ', Ch
  Write(*,*)'Point of minima: ', X
  Write(*,*)'Check:           ', Tan(X(1)), Cos(X(2)) 

  X(:) = -20.0D0
  CALL Misimplex(Func, X, Ch)
  Write(*,*)'= SIMPLEX ='
  Write(*,*)'Function value:  ', Ch
  Write(*,*)'Point of minima: ', X
  Write(*,*)'Check:           ', Tan(X(1)), Cos(X(2)) 

  X(:) = -20.0D0
  CALL Miseek(Func, X, Ch)
  Write(*,*)'= SEEK ='
  Write(*,*)'Function value:  ', Ch
  Write(*,*)'Point of minima: ', X
  Write(*,*)'Check:           ', Tan(X(1)), Cos(X(2)) 

  X(:) = -20.0D0
  CALL Miseek(Func, X, Ch)
  Write(*,*)'= SCAN ='
  Write(*,*)'Function value:  ', Ch
  Write(*,*)'Point of minima: ', X
  Write(*,*)'Check:           ', Tan(X(1)), Cos(X(2)) 

  Stop
End Program TestAPI

Function Func(X)

  Real (kind=8), Intent (in) :: X(:)
  Real (kind=8) :: Func

  Func = (X(1)-tan(X(1)))**2 + (X(2) - Cos(X(2)))**2

  Return
End Function Func
