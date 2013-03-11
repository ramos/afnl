
Program Checksum

  USE NonNumeric

  IMPLICIT NONE

  Integer, Allocatable :: buf(:)
  Character (len=10000) :: fn, cmm
  Logical :: fex, chk
  Integer :: Ihsh, Isz, I, Narg


  Narg = Command_argument_count()
  If (Narg == 0) Then
     CALL Help()
     Stop
  End If

  I = 0
  Do 
     I = I+1
     CALL Get_command_argument(I, fn)

     Inquire(FILE=Trim(fn), EXIST=fex)
     If (fex) Then
        Write(*,'(1A)')'##'
        Write(*,'(1X,2A)')'File:    ', Trim(fn)
     Else
        CALL Abort('chkfile', 'File '//Trim(fn)//' does not exists')
     End If
     
     Write(*,'(1A)')'#'
     CALL ReadBuffer(buf, Trim(fn), chk, Ihsh)
     If (chk) Then
        Write(*,'(1X,1A,1Z9,1A)')'Hash:    ', Ihsh, ' checksum *OK*'
     Else
        Write(*,'(1X,1A,1Z9,1A)')'Hash:         ', Ihsh, ' checksum *FAILED*'
        Write(*,'(1X,1A,1Z9,1A)')'Hash on file: ', GetHash(fn)
     End If
     
     Write(*,'(1X,2A)')'Date:    ', GetDT(fn)
     
     Open (Unit=69, File=Trim(fn), Form='UNFORMATTED',&
          & Access='STREAM', Action='READ')
     Inquire(Unit=69, Size=ISZ)
     Close(69)
     If (Isz >= 1073741824) Then
        Write(*,'(1X,1A,1F9.2,1A)')'Size:    ', Real(Isz)/1073741824.&
             &, ' GBytes'
     Else If (Isz > 1048576) Then
        Write(*,'(1X,1A,1F9.2,1A)')'Size:    ', Real(Isz)/1048576., ' MBytes'
     Else If (Isz > 1024) Then
        Write(*,'(1X,1A,1F9.2,1A)')'Size:    ', Real(Isz)/1024., ' KBytes'
     Else 
        Write(*,'(1X,1A,1I12,1A)')'Size:    ', Isz, ' Bytes'
     End If
     
     CALL GetComment(fn, cmm)
     If (Trim(cmm) /= "") Then
        Write(*,'(1X,2A)')'Comment: ', Trim(cmm)
     End If
     
     Write(*,'(1A)')'##'
     Write(*,'(1A)')

     If (I==Narg) Stop
  End Do

  Stop

CONTAINS

  ! ***********************
  Subroutine help() 
  ! ***********************

    Write(*,*)
    Write(*,*)' Usage: chkfile -i <filename>'
    Write(*,*)'        Checks hash and output info about file.'
    Write(*,*)

    Return
  End Subroutine help

End Program Checksum
