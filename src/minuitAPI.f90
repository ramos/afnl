
MODULE MinuitAPI

!  USE MODMinuit
  USE ISO_C_BINDING

  IMPLICIT NONE

  Interface Fit
     Module Procedure Fit1D, FitMD, Fit1DCorr, FitMDCorr
  End Interface

  Interface Minimize
     Module Procedure MinimizeMD, MinimizeMD_Bounds
  End Interface

  Real (kind=8), Allocatable :: ShX1d(:), ShY(:), ShXMd(:,:),&
       & ShYerr(:), ShInvC(:,:)

  Integer :: NcP

  Private Fit1D, FitMD, Fit1DCorr, FitMDCorr, MinimizeMD, Chisqr1d,&
       & ChisqrMd, Chisqr1dCorr, ChisqrMdCorr, Fm, MinimizeMD_Bounds

CONTAINS

! ***********************************************
! *
  Subroutine Fit1D(X, Y, Yerr, Func, Coef, Chisqr, logfile)
! *
! ***********************************************

    Real (kind=8), Intent (in) :: X(:), Y(:), Yerr(:)
    Real (kind=8), Intent (inout) :: Coef(:)
    Real (kind=8), Intent (out) :: Chisqr
    Character (len=*), Optional :: logfile
    
    Integer :: Idummy, I, Ierr, Ifoo
    Real (kind=8) :: foo
    Character (len=1) :: cfoo

    Interface 
       Function Func(X, C)
         Real (kind=8), Intent (in) :: X
         Real (kind=8), Intent (in) :: C(:)
         Real (kind=8) :: Func
       End Function Func
    End Interface

    If (Present(logfile)) Then
       Open (unit=69, File=Trim(logfile))
    Else
       Open (unit=69, File='minuit.log')
    End If
    CALL MnInit(5,69,69)
    CALL Mncomd(Chisqr1D,"set pri -1",Ierr, Func)
    CALL Mncomd(Chisqr1D,"set now",Ierr, Func)
    CALL Mncomd(Chisqr1D,"set str 2",Ierr, Func)
    
    Allocate(ShX1d(Size(X)), ShY(Size(Y)), ShYerr(Size(Y)))
    ShX1d(:) = X(:)
    ShY(:) = Y(:)
    ShYerr(:) = Yerr(:)

    Do I = 1, Size(Coef)
       CALL MnParm(I, 'X', Coef(I), 10.0D0, 0.0D0, 0.0D0, Ierr)
    End Do

    CALL Mncomd(Chisqr1D,"mini",Ierr,Func)
    CALL Mncomd(Chisqr1D,"seek",Ierr,Func)
    CALL Mncomd(Chisqr1D,"mini",Ierr,Func)
    
!    Output the Chisqr, and parameter value
    CALL MNstat(Chisqr, foo, foo, Ifoo, Ifoo, Ierr)
    Do I = 1, Size(Coef)
       CALL MnPout(I, cfoo, Coef(I), foo, foo, foo, Ierr)
    End Do

    DeAllocate(ShX1d, ShY, ShYerr)
    Close(69)

    Return
  End Subroutine Fit1D

  Subroutine Chisqr1d(Npar,Grad,Fval,Xval,Iflag,Func)
    
    Integer, Intent (in) :: Npar, Iflag
    Real (kind=8), Intent (out) :: Fval
    Real (kind=8), Intent (in) :: Xval(*), Grad(*)
        
    Integer :: I
    Real (kind=8) :: C(Npar)
    
    Interface 
       Function Func(X, C)
         Real (kind=8), Intent (in) :: X
         Real (kind=8), Intent (in) :: C(:)
         Real (kind=8) :: Func
       End Function Func
    End Interface

    Do I = 1, Npar
       C(I) = Xval(I)
    End Do
    
    Fval = 0.0D0
    Do I = 1, Size(ShX1d)
       FVal = Fval + (ShY(I) - Func(ShX1d(I),C))**2 / ShYerr(I)**2 
    End Do
    
    Return
  End Subroutine Chisqr1d


! ***********************************************
! *
  Subroutine FitMD(X, Y, Yerr, Func, Coef, Chisqr, logfile)
! *
! ***********************************************

    Real (kind=8), Intent (in) :: X(:,:), Y(:), Yerr(:)
    Real (kind=8), Intent (inout) :: Coef(:)
    Real (kind=8), Intent (out) :: Chisqr
    Character (len=*), Optional :: logfile
    
    Integer :: Idummy, I, Ierr, Ifoo
    Real (kind=8) :: foo
    Character (len=1) :: cfoo

    Interface 
       Function Func(X, C)
         Real (kind=8), Intent (in) :: X(:)
         Real (kind=8), Intent (in) :: C(:)
         Real (kind=8) :: Func
       End Function Func
    End Interface

    If (Present(logfile)) Then
       Open (unit=69, File=Trim(logfile))
    Else
       Open (unit=69, File='minuit.log')
    End If
    CALL MnInit(5,69,69)
    CALL Mncomd(ChisqrMD,"set pri -1",Ierr, Func)
    CALL Mncomd(ChisqrMD,"set now",Ierr, Func)
    CALL Mncomd(ChisqrMD,"set str 2",Ierr, Func)
    
    Allocate(ShXMd(Size(X,1),Size(X,2)), ShY(Size(Y)), ShYerr(Size(Y)))
    ShXMd(:,:) = X(:,:)
    ShY(:) = Y(:)
    ShYerr(:) = Yerr(:)

    Do I = 1, Size(Coef)
       CALL MnParm(I, 'X', Coef(I), 10.0D0, 0.0D0, 0.0D0, Ierr)
    End Do

    CALL Mncomd(ChisqrMD,"mini",Ierr,Func)
    CALL Mncomd(ChisqrMD,"seek",Ierr,Func)
    CALL Mncomd(ChisqrMD,"mini",Ierr,Func)
    
!    Output the Chisqr, and parameter value
    CALL MNstat(Chisqr, foo, foo, Ifoo, Ifoo, Ierr)
    Do I = 1, Size(Coef)
       CALL MnPout(I, cfoo, Coef(I), foo, foo, foo, Ierr)
    End Do

    DeAllocate(ShXMd, ShY, ShYerr)
    Close(69)

    Return
  End Subroutine FitMD

  Subroutine ChisqrMd(Npar,Grad,Fval,Xval,Iflag,Func)
    
    Integer, Intent (in) :: Npar, Iflag
    Real (kind=8), Intent (out) :: Fval
    Real (kind=8), Intent (in) :: Xval(*), Grad(*)
        
    Integer :: I
    Real (kind=8) :: C(Npar)
    
    Interface 
       Function Func(X, C)
         Real (kind=8), Intent (in) :: X(:)
         Real (kind=8), Intent (in) :: C(:)
         Real (kind=8) :: Func
       End Function Func
    End Interface

    Do I = 1, Npar
       C(I) = Xval(I)
    End Do
    
    Fval = 0.0D0
    Do I = 1, Size(ShXMD,1)
       FVal = Fval + (ShY(I) - Func(ShXMD(I,:),C))**2 / ShYerr(I)**2 
    End Do

    Return
  End Subroutine ChisqrMd

! ***********************************************
! *
  Subroutine Fit1DCorr(X, Y, InvCorr, Func, Coef, Chisqr, logfile)
! *
! ***********************************************

    Real (kind=8), Intent (in) :: X(:), Y(:), InvCorr(:,:)
    Real (kind=8), Intent (inout) :: Coef(:)
    Real (kind=8), Intent (out) :: Chisqr
    Character (len=*), Optional :: logfile
    
    Integer :: Idummy, I, Ierr, Ifoo
    Real (kind=8) :: foo
    Character (len=1) :: cfoo

    Interface 
       Function Func(X, C)
         Real (kind=8), Intent (in) :: X
         Real (kind=8), Intent (in) :: C(:)
         Real (kind=8) :: Func
       End Function Func
    End Interface

    If (Present(logfile)) Then
       Open (unit=69, File=Trim(logfile))
    Else
       Open (unit=69, File='minuit.log')
    End If
    CALL MnInit(5,69,69)
    CALL Mncomd(Chisqr1DCorr,"set pri -1",Ierr, Func)
    CALL Mncomd(Chisqr1DCorr,"set now",Ierr, Func)
    CALL Mncomd(Chisqr1DCorr,"set str 2",Ierr, Func)
    
    Allocate(ShX1d(Size(X)), ShY(Size(Y)), &
         & ShInvC(Size(InvCorr,1), Size(InvCorr,2)))
    ShX1d(:) = X(:)
    ShY(:) = Y(:)
    ShInvC(:,:) = InvCorr(:,:)

    Do I = 1, Size(Coef)
       CALL MnParm(I, 'X', Coef(I), 10.0D0, 0.0D0, 0.0D0, Ierr)
    End Do

    CALL Mncomd(Chisqr1DCorr,"mini",Ierr,Func)
    CALL Mncomd(Chisqr1DCorr,"seek",Ierr,Func)
    CALL Mncomd(Chisqr1DCorr,"mini",Ierr,Func)
    
!    Output the Chisqr, and parameter value
    CALL MNstat(Chisqr, foo, foo, Ifoo, Ifoo, Ierr)
    Do I = 1, Size(Coef)
       CALL MnPout(I, cfoo, Coef(I), foo, foo, foo, Ierr)
    End Do

    DeAllocate(ShX1d, ShY, ShInvC)
    Close(69)

    Return
  End Subroutine Fit1DCorr

  Subroutine Chisqr1dCorr(Npar,Grad,Fval,Xval,Iflag,Func)
    
    Integer, Intent (in) :: Npar, Iflag
    Real (kind=8), Intent (out) :: Fval
    Real (kind=8), Intent (in) :: Xval(*), Grad(*)
        
    Integer :: I, J
    Real (kind=8) :: C(Npar)
    
    Interface 
       Function Func(X, C)
         Real (kind=8), Intent (in) :: X
         Real (kind=8), Intent (in) :: C(:)
         Real (kind=8) :: Func
       End Function Func
    End Interface

    Do I = 1, Npar
       C(I) = Xval(I)
    End Do
    
    Fval = 0.0D0
    Do J = 1, Size(ShX1D)
       Do I = 1, Size(ShX1D)
          FVal = Fval + (ShY(I) - Func(ShX1D(I),C)) * ShInvC(I,J) * &
               & (ShY(J) - Func(ShX1D(J),C))
       End Do
    End Do

    Return
  End Subroutine Chisqr1dCorr

! ***********************************************
! *
  Subroutine FitMDCorr(X, Y, InvCorr, Func, Coef, Chisqr, logfile)
! *
! ***********************************************

    Real (kind=8), Intent (in) :: X(:,:), Y(:), InvCorr(:,:)
    Real (kind=8), Intent (inout) :: Coef(:)
    Real (kind=8), Intent (out) :: Chisqr
    Character (len=*), Optional :: logfile
    
    Integer :: Idummy, I, Ierr, Ifoo
    Real (kind=8) :: foo
    Character (len=1) :: cfoo

    Interface 
       Function Func(X, C)
         Real (kind=8), Intent (in) :: X(:)
         Real (kind=8), Intent (in) :: C(:)
         Real (kind=8) :: Func
       End Function Func
    End Interface

    If (Present(logfile)) Then
       Open (unit=69, File=Trim(logfile))
    Else
       Open (unit=69, File='minuit.log')
    End If
    CALL MnInit(5,69,69)
    CALL Mncomd(ChisqrMDCorr,"set pri -1",Ierr, Func)
    CALL Mncomd(ChisqrMDCorr,"set now",Ierr, Func)
    CALL Mncomd(ChisqrMDCorr,"set str 2",Ierr, Func)
    
    Allocate(ShXMd(Size(X,1),Size(X,2)), ShY(Size(Y)), &
         & ShInvC(Size(InvCorr,1), Size(InvCorr,2)))
    ShXMd(:,:) = X(:,:)
    ShY(:) = Y(:)
    ShInvC(:,:) = InvCorr(:,:)

    Do I = 1, Size(Coef)
       CALL MnParm(I, 'X', Coef(I), 10.0D0, 0.0D0, 0.0D0, Ierr)
    End Do

    CALL Mncomd(ChisqrMDCorr,"mini",Ierr,Func)
    CALL Mncomd(ChisqrMDCorr,"seek",Ierr,Func)
    CALL Mncomd(ChisqrMDCorr,"mini",Ierr,Func)
    
!    Output the Chisqr, and parameter value
    CALL MNstat(Chisqr, foo, foo, Ifoo, Ifoo, Ierr)
    Do I = 1, Size(Coef)
       CALL MnPout(I, cfoo, Coef(I), foo, foo, foo, Ierr)
    End Do

    DeAllocate(ShXMd, ShY, ShInvC)
    Close(69)

    Return
  End Subroutine FitMDCorr

  Subroutine ChisqrMdCorr(Npar,Grad,Fval,Xval,Iflag,Func)
    
    Integer, Intent (in) :: Npar, Iflag
    Real (kind=8), Intent (out) :: Fval
    Real (kind=8), Intent (in) :: Xval(*), Grad(*)
        
    Integer :: I, J
    Real (kind=8) :: C(Npar)
    
    Interface 
       Function Func(X, C)
         Real (kind=8), Intent (in) :: X(:)
         Real (kind=8), Intent (in) :: C(:)
         Real (kind=8) :: Func
       End Function Func
    End Interface

    Do I = 1, Npar
       C(I) = Xval(I)
    End Do
    
    Fval = 0.0D0
    Do J = 1, Size(ShX1D,1)
       Do I = 1, Size(ShX1D,1)
          FVal = Fval + (ShY(I) - Func(ShXMD(I,:),C)) * ShInvC(I,J) * &
               & (ShY(J) - Func(ShXMD(J,:),C))
       End Do
    End Do

    Return
  End Subroutine ChisqrMdCorr

! ***********************************************
! *
  Subroutine MinimizeMD(Func, X, Fval, Fparms, logfile)
! *
! ***********************************************

    Real (kind=8), Intent (inout) :: X(:)
    Real (kind=8), Intent (out) :: Fval
    Integer, Intent (in), Optional :: Fparms(0:,1:)
    Character (len=*), Intent (in), Optional :: logfile

    Integer :: Idummy, I, Ierr, Ifoo, J
    Real (kind=8) :: foo
    Character (len=1) :: cfoo
    Logical :: Empty

    Interface 
       Function Func(X)
         Real (kind=8), Intent (in) :: X(:)
         Real (kind=8) :: Func
       End Function Func
    End Interface

    If (Present(logfile)) Then
       Open (unit=69, File=Trim(logfile))
    Else
       Open (unit=69, File='minuit.log')
    End If
    CALL MnInit(5,69,69)
    CALL Mncomd(Fm,"set pri -1",Ierr, Func)
    CALL Mncomd(Fm,"set now",Ierr, Func)
    CALL Mncomd(Fm,"set str 2",Ierr, Func)
    
    Do I = 1, Size(X)
       If (Present(Fparms)) Then
!!$          Empty = .True.
!!$          Do J = 1, Fparms(0,1)
!!$             If (Fparms(1,J) == I) Then
!!$                empty = .False.
!!$             End If
!!$          End Do
!!$          If (Empty) Then
!!$             CALL MnParm(I, 'X', 0.0D0, 1.0D-4, 0.0D0, 0.0D0, Ierr)
!!$          Else
          CALL MnParm(I, 'X', X(I), 10.0D-3, 0.0D0, 0.0D0, Ierr)
!!$          End If
          CALL Mnfix(I)
       Else
          CALL MnParm(I, 'X', X(I), 10.0D-3, 0.0D0, 0.0D0, Ierr)
       End If
    End Do
    

    If (Present(Fparms)) Then
       Do I = 1, Size(Fparms,1)-1
          CALL MNstat(Fval, foo, foo, Ifoo, Ierr, Ierr)
          Do J = 1, Fparms(0,I)
             CALL MNfree(Fparms(I,J))
          End Do
          ! Determine current number of variable parameters
          CALL MNstat(Fval, foo, foo, NcP, Ierr, Ierr)
          ! Minimize
          If (I == 1) Then
             CALL Mncomd(Fm,"mini",Ierr,Func)
             CALL Mncomd(Fm,"seek",Ierr,Func)
          End If
          CALL Mncomd(Fm,"migrad",Ierr,Func)
          CALL Mncomd(Fm,"migrad",Ierr,Func)
       End Do
    Else
       ! Determine current number of variable parameters
       CALL MNstat(Fval, foo, foo, NcP, Ierr, Ierr)
       ! Minimize
       CALL Mncomd(Fm,"mini",Ierr,Func)
       CALL Mncomd(Fm,"seek",Ierr,Func)
       CALL Mncomd(Fm,"migrad",Ierr,Func)
       CALL Mncomd(Fm,"migrad",Ierr,Func)
    End If
    


!    Output function value, and position of the minima
    CALL MNstat(Fval, foo, foo, Ifoo, Ifoo, Ierr)
    Do I = 1, Size(X)
       CALL MnPout(I, cfoo, X(I), foo, foo, foo, Ierr)
    End Do

    Close(69)    

    Return
  End Subroutine MinimizeMD

! ***********************************************
! *
  Subroutine MinimizeMD_Bounds(Func, X, Fval, Bounds, Fparms, logfile)
! *
! ***********************************************

    Real (kind=8), Intent (inout) :: X(:)
    Real (kind=8), Intent (out) :: Fval
    Real (kind=8), Intent (in) :: Bounds(:)
    Integer, Intent (in), Optional :: Fparms(0:,1:)
    Character (len=*), Intent (in), Optional :: logfile

    Integer :: Idummy, I, Ierr, Ifoo, J, Ndim
    Real (kind=8) :: foo
    Character (len=1) :: cfoo
    Logical :: Empty

    Interface 
       Function Func(X)
         Real (kind=8), Intent (in) :: X(:)
         Real (kind=8) :: Func
       End Function Func
    End Interface

    If (Present(logfile)) Then
       Open (unit=69, File=Trim(logfile))
    Else
       Open (unit=69, File='minuit.log')
    End If
    CALL MnInit(5,69,69)
    CALL Mncomd(Fm,"set pri -1",Ierr, Func)
    CALL Mncomd(Fm,"set now",Ierr, Func)
    CALL Mncomd(Fm,"set str 2",Ierr, Func)
    

    Ndim = Size(X)
    Do I = 1, Ndim
       If (Present(Fparms)) Then
!!$          Empty = .True.
!!$          Do J = 1, Fparms(0,1)
!!$             If (Fparms(1,J) == I) Then
!!$                empty = .False.
!!$             End If
!!$          End Do
!!$          If (Empty) Then
!!$             CALL MnParm(I, 'X', 0.0D0, 1.0D-4, 0.0D0, 0.0D0, Ierr)
!!$          Else
          CALL MnParm(I, 'X', X(I), 10.0D-3, Bounds(2*I-1), Bounds(2*I), Ierr)
!!$          End If
          CALL Mnfix(I)
       Else
          CALL MnParm(I, 'X', X(I), 10.0D-3, Bounds(2*I-1), Bounds(2*I), Ierr)
       End If
    End Do
    

    If (Present(Fparms)) Then
       Do I = 1, Size(Fparms,1)-1
          CALL MNstat(Fval, foo, foo, Ifoo, Ierr, Ierr)
          Do J = 1, Fparms(0,I)
             CALL MNfree(Fparms(I,J))
          End Do
          ! Determine current number of variable parameters
          CALL MNstat(Fval, foo, foo, NcP, Ierr, Ierr)
          ! Minimize
          If (I == 1) Then
             CALL Mncomd(Fm,"mini",Ierr,Func)
             CALL Mncomd(Fm,"seek",Ierr,Func)
          End If
          CALL Mncomd(Fm,"migrad",Ierr,Func)
          CALL Mncomd(Fm,"migrad",Ierr,Func)
       End Do
    Else
       ! Determine current number of variable parameters
       CALL MNstat(Fval, foo, foo, NcP, Ierr, Ierr)
       ! Minimize
       CALL Mncomd(Fm,"mini",Ierr,Func)
       CALL Mncomd(Fm,"seek",Ierr,Func)
       CALL Mncomd(Fm,"migrad",Ierr,Func)
       CALL Mncomd(Fm,"migrad",Ierr,Func)
    End If
    


!    Output function value, and position of the minima
    CALL MNstat(Fval, foo, foo, Ifoo, Ifoo, Ierr)
    Do I = 1, Size(X)
       CALL MnPout(I, cfoo, X(I), foo, foo, foo, Ierr)
    End Do

    Close(69)    

    Return
  End Subroutine MinimizeMD_Bounds

  Subroutine Fm(Npar,Grad,Fval,Xval,Iflag,Func)
    
    Integer, Intent (in) :: Npar, Iflag
    Real (kind=8), Intent (out) :: Fval
    Real (kind=8), Intent (in) :: Xval(*), Grad(*)
    
    Integer :: I
    Real (kind=8), Allocatable :: C(:)
    
    Interface 
       Function Func(X)
         Real (kind=8), Intent (in) :: X(:)
         Real (kind=8) :: Func
       End Function Func
    End Interface

    Allocate(C(1:NcP))
    Do I = 1, NcP
       C(I) = Xval(I)
    End Do

    Fval = Func(C(:))
    Deallocate(C)

    Return
  End Subroutine Fm

  SUBROUTINE Mnfix(i)
! call minuit fix
! that fixes the real parameter
! in multifit
    INTEGER, INTENT(IN) :: i

    CHARACTER(len=20) :: mss
    INTEGER :: err

    WRITE (mss,'(I4)') i
    CALL Mncomd(Fm,"fix"//TRIM(mss),err,0)
  END SUBROUTINE Mnfix

  SUBROUTINE Mnfree(i)
! call minuit free
! that fixes the real parameter
! in multifit
    INTEGER, INTENT(IN) :: i

    CHARACTER(len=20) :: mss
    INTEGER :: err

    WRITE (mss,'(I4)') i
    CALL Mncomd(Fm,"rel"//TRIM(mss),err,0)
  END SUBROUTINE Mnfree

END MODULE MinuitAPI
