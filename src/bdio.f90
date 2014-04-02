!
! A MODULE for non numerical routines (sorting and locate)
!
! Copyright (C) 2005  Alberto Ramos <alberto@martin.ft.uam.es>
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
! USA
! 

! ***************************************************
! *
MODULE ModBDio
! *
! ***************************************************
! *
! * NonNumeric routines for sorting and locating
! * data.
! *
! ***************************************************

  USE NonNumeric

  IMPLICIT NONE


  Integer, Parameter :: BDIO_R_MODE=0, BDIO_W_MODE=1, BDIO_A_MODE=2, &
       & MAXFNAME = 4096

  Integer :: BDIO_MAGIC, BDIO_VERSION, BDIO_BIN_GENERIC, &
       & BDIO_ASC_EXEC, BDIO_BIN_INT32BE, BDIO_BIN_INT32LE, &
       & BDIO_BIN_INT64BE, BDIO_BIN_INT64LE, BDIO_BIN_F32BE, &
       & BDIO_BIN_F32LE, BDIO_BIN_F64BE, BDIO_BIN_F64LE, &
       & BDIO_ASC_GENERIC, BDIO_ASC_XML, BDIO_BIN_INT32, &
       & BDIO_BIN_INT64, BDIO_BIN_F32, BDIO_BIN_F64, &
       & BDIO_H_STATE, BDIO_R_STATE, &
       & BDIO_N_STATE, BDIO_E_STATE, BDIO_LEND, BDIO_BEND, &
       & BDIO_MAX_RECORD_LENGTH, BDIO_MAX_LONG_RECORD_LENGTH, &
       & BDIO_BUF_SIZE, BDIO_MAX_HOST_LENGTH, BDIO_MAX_USER_LENGTH, &
       & BDIO_MAX_PINFO_LENGTH

  Logical :: DEFAULT_HASH_CHECK = .False.

  Data BDIO_MAGIC /Z'7ffbd07e'/
  Data BDIO_VERSION /1/
  
  ! Record data formats
  DATA BDIO_BIN_GENERIC /Z'00'/
  DATA BDIO_ASC_EXEC    /Z'01'/
  DATA BDIO_BIN_INT32BE /Z'02'/
  DATA BDIO_BIN_INT32LE /Z'03'/
  DATA BDIO_BIN_INT64BE /Z'04'/
  DATA BDIO_BIN_INT64LE /Z'05'/
  DATA BDIO_BIN_F32BE   /Z'06'/
  DATA BDIO_BIN_F32LE   /Z'07'/
  DATA BDIO_BIN_F64BE   /Z'08'/
  DATA BDIO_BIN_F64LE   /Z'09'/
  DATA BDIO_ASC_GENERIC /Z'0A'/
  DATA BDIO_ASC_XML     /Z'0B'/

  DATA BDIO_BIN_INT32   /Z'F0'/
  DATA BDIO_BIN_INT64   /Z'F1'/
  DATA BDIO_BIN_F32     /Z'F2'/
  DATA BDIO_BIN_F64     /Z'F3'/

  DATA  BDIO_H_STATE /1/
  DATA  BDIO_R_STATE /2/
  DATA  BDIO_N_STATE /3/
  DATA  BDIO_E_STATE /4/

  DATA  BDIO_LEND /0/
  DATA  BDIO_BEND /1/

  DATA  BDIO_MAX_RECORD_LENGTH /1048575/
  DATA  BDIO_MAX_LONG_RECORD_LENGTH /268435455 /

  DATA  BDIO_BUF_SIZE /7777/

  DATA  BDIO_MAX_HOST_LENGTH /256/

  DATA  BDIO_MAX_USER_LENGTH /33/

  DATA  BDIO_MAX_PINFO_LENGTH /3505/
  
  Type :: BDIO
     Integer :: ifn, rcnt=0, hcnt=0, tcnt=0
     Integer :: imode ! r/w/a
     Integer (kind=8) :: startfile=-1, endfile=-1, rwpos=-1
     Logical :: lendian, opened
     
     Type (BDIO_record), pointer :: first => null(), last => null(), &
          & current => null()
  End type BDIO

  Type :: BDIO_record
     logical :: ishdr, islong
     Integer (kind=8) :: rlen, rpos, rend
     Integer :: rfmt = -1, ruinfo = -1, rid = -1, iver = -1, &
          & hash = 314159265

     Type (BDIO_record), Pointer :: next => null(), prev => null()
  End type BDIO_RECORD


  Interface byteswap
     Module Procedure byteswap_int32, byteswap_int32V, byteswap_R, &
          & byteswap_RV, byteswap_DV, byteswap_D, byteswap_DZ, &
          & byteswap_DZV, byteswap_ZV, byteswap_Z, &
          & byteswap_int64, byteswap_int64V
  End Interface byteswap

  Interface BDIO_write
     Module Procedure BDIO_write_i32, BDIO_write_i64, &
          & BDIO_write_f32, BDIO_write_f64, BDIO_write_z32, &
          & BDIO_write_z64
  End Interface BDIO_write
  
  Interface BDIO_read
     Module Procedure BDIO_read_i32, BDIO_read_i64, &
          & BDIO_read_f32, BDIO_read_f64, BDIO_read_z32, &
          & BDIO_read_z64
  End Interface BDIO_read
  

  Private :: BDIO_MAGIC, BDIO_VERSION, BDIO_R_MODE, &
       & BDIO_W_MODE, BDIO_A_MODE, BDIO_LEND, BDIO_BEND, &
       & BDIO_MAX_RECORD_LENGTH, BDIO_MAX_LONG_RECORD_LENGTH, &
       & BDIO_BUF_SIZE, BDIO_MAX_HOST_LENGTH, BDIO_MAX_USER_LENGTH, &
       & BDIO_MAX_PINFO_LENGTH, byteswap_int32, byteswap_int32V, &
       & byteswap_R, byteswap_RV, byteswap_DV, byteswap_D, byteswap_DZ, &
       & byteswap_DZV, byteswap_ZV, byteswap_Z, byteswap_int64, &
       & byteswap_int64V, Rewrite_rlen, BDIO_read_f32, BDIO_read_f64, &
       & BDIO_read_i32, BDIO_read_i64, BDIO_read_z32, BDIO_read_z64, &
       & BDIO_write_i32, BDIO_write_i64, BDIO_write_f32, &
       & BDIO_write_f64, BDIO_write_z32, BDIO_write_z64

CONTAINS

! ********************************
! *
    Subroutine BDIO_close(fbd)
! *
! ********************************
      
      Type (BDIO), Intent (inout) :: fbd
      Type (BDIO_record), pointer :: p, n

      Integer :: kont

      p => fbd%first

      kont = 0
      Do 
         If (.not.Associated(p)) Exit
         n => p%next
         Deallocate(p)
         p => n
         kont = kont+1
      End Do
      Close(fbd%ifn)

      Return
    End Subroutine BDIO_close

! ********************************
! *
    Subroutine BDIO_error(fbd, routine, msg)
! *
! ********************************

      Type (BDIO), Intent (inout) :: fbd
      Character (len=*), Intent (in) :: routine, msg

      Write(0,*)'In '//Trim(routine)// ' :'
      Write(0,'(5X,1A)')Trim(msg)
      Write(0,*)
      Write(0,*)'Current Record: ', fbd%current%rid
      Stop
      
      Return
    End Subroutine BDIO_error

! ********************************
! *
    Function BDIO_write_i32(fbd,ibuf)
! *
! ********************************

      Type (BDIO), Intent (inout) :: fbd
      Integer (kind=4), Intent (inout) :: ibuf(:)
      Integer :: BDIO_write_i32, ios

      Integer (kind=8) :: iln
      Type (BDIO_record), pointer :: p

      BDIO_write_i32 = -1
      fbd%current => fbd%last
      fbd%rwpos = fbd%last%rend
      p => fbd%current
      If (  (p%rfmt /= BDIO_BIN_INT32LE).and.&
           &(p%rfmt /= BDIO_BIN_INT32BE).and.&
           &(p%rfmt /= BDIO_BIN_INT32  ) ) Then
         Call BDIO_error(fbd,'BDIO_Write_i32', &
              & 'Incorrect data type') 
      End If

      If (  ( (p%rfmt == BDIO_BIN_INT32BE).and.(fbd%lendian) ).or. &
           &( (p%rfmt == BDIO_BIN_INT32LE).and.(.not.fbd%lendian) ) ) &
           & Then
         CALL ByteSwap(ibuf)
         Write(fbd%ifn,Pos=fbd%rwpos,Iostat=ios)ibuf
         CALL ByteSwap(ibuf)
      Else
         Write(fbd%ifn,Pos=fbd%rwpos,Iostat=ios)ibuf
      End If
      
      iln = fbd%current%rlen + 4_8*Size(ibuf)
      CALL Rewrite_rlen(fbd, iln)

      If (ios == 0) BDIO_write_i32 = 4*Size(ibuf)

      Return
    End Function BDIO_write_i32

! ********************************
! *
    Function BDIO_write_f32(fbd,buf)
! *
! ********************************

      Type (BDIO), Intent (inout) :: fbd
      Real (kind=4), Intent (inout) :: buf(:)
      Integer :: BDIO_write_f32, ios

      Integer (kind=8) :: iln
      Type (BDIO_record), pointer :: p

      BDIO_write_f32 = -1
      fbd%current => fbd%last
      fbd%rwpos = fbd%last%rend
      p => fbd%current
      If (  (p%rfmt /= BDIO_BIN_F32LE).and.&
           &(p%rfmt /= BDIO_BIN_F32BE).and.&
           &(p%rfmt /= BDIO_BIN_F32  ) ) Then
         Call BDIO_error(fbd,'BDIO_Write_f32', &
              & 'Incorrect data type') 
      End If

      If (  ( (p%rfmt == BDIO_BIN_F32BE).and.(fbd%lendian) ).or. &
           &( (p%rfmt == BDIO_BIN_F32LE).and.(.not.fbd%lendian) ) ) &
           & Then
         CALL ByteSwap(buf)
         Write(fbd%ifn,Pos=fbd%rwpos,Iostat=ios)buf
         CALL ByteSwap(buf)
      Else
         Write(fbd%ifn,Pos=fbd%rwpos,Iostat=ios)buf
      End If

      iln = fbd%current%rlen + 4_8*Size(buf)
      CALL Rewrite_rlen(fbd, iln)

      If (ios == 0) BDIO_write_f32 = 4*Size(buf)

      Return
    End Function BDIO_write_f32

! ********************************
! *
    Function BDIO_write_i64(fbd,ibuf)
! *
! ********************************

      Type (BDIO), Intent (inout) :: fbd
      Integer (kind=8), Intent (inout) :: ibuf(:)
      Integer :: BDIO_write_i64, ios

      Integer (kind=8) :: iln
      Type (BDIO_record), pointer :: p

      BDIO_write_i64 = -1
      fbd%current => fbd%last
      fbd%rwpos = fbd%last%rend
      p => fbd%current
      If (  (p%rfmt /= BDIO_BIN_INT64LE).and.&
           &(p%rfmt /= BDIO_BIN_INT64BE).and.&
           &(p%rfmt /= BDIO_BIN_INT64  ) ) Then
         Call BDIO_error(fbd,'BDIO_Write_i64', &
              & 'Incorrect data type') 
      End If

      If (  ( (p%rfmt == BDIO_BIN_INT64BE).and.(fbd%lendian) ).or. &
           &( (p%rfmt == BDIO_BIN_INT64LE).and.(.not.fbd%lendian) ) ) &
           & Then
         CALL ByteSwap(ibuf)
         Write(fbd%ifn,Pos=fbd%rwpos,Iostat=ios)ibuf
         CALL ByteSwap(ibuf)
      Else
         Write(fbd%ifn,Pos=fbd%rwpos,Iostat=ios)ibuf
      End If
      
      iln = fbd%current%rlen + 8_8*Size(ibuf)
      CALL Rewrite_rlen(fbd, iln)

      If (ios == 0) BDIO_write_i64 = 8*Size(ibuf)

      Return
    End Function BDIO_write_i64

! ********************************
! *
    Function BDIO_write_f64(fbd,ibuf)
! *
! ********************************

      Type (BDIO), Intent (inout) :: fbd
      Real (kind=8), Intent (inout) :: ibuf(:)
      Integer :: BDIO_write_f64, ios

      Integer (kind=8) :: iln
      Type (BDIO_record), pointer :: p

      BDIO_write_f64 = -1
      fbd%current => fbd%last
      fbd%rwpos = fbd%last%rend
      p => fbd%current
      If (  (p%rfmt /= BDIO_BIN_F64LE).and.&
           &(p%rfmt /= BDIO_BIN_F64BE).and.&
           &(p%rfmt /= BDIO_BIN_F64  ) ) Then
         Call BDIO_error(fbd,'BDIO_Write_f64', &
              & 'Incorrect data type') 
      End If

      If (  ( (p%rfmt == BDIO_BIN_F64BE).and.(fbd%lendian) ).or. &
           &( (p%rfmt == BDIO_BIN_F64LE).and.(.not.fbd%lendian) ) ) &
           & Then
         CALL ByteSwap(ibuf)
         Write(fbd%ifn,Pos=fbd%rwpos,Iostat=ios)ibuf
         CALL ByteSwap(ibuf)
      Else
         Write(fbd%ifn,Pos=fbd%rwpos,Iostat=ios)ibuf
      End If
      
      iln = fbd%current%rlen + 8_8*Size(ibuf)
      CALL Rewrite_rlen(fbd, iln)

      If (ios == 0) BDIO_write_f64 = 8*Size(ibuf)

      Return
    End Function BDIO_write_f64
    
! ********************************
! *
    Function BDIO_write_z64(fbd,buf)
! *
! ********************************

      Type (BDIO), Intent (inout) :: fbd
      Complex (kind=8), Intent (inout) :: buf(:)
      Integer :: BDIO_write_z64, ios

      Integer (kind=8) :: iln
      Type (BDIO_record), pointer :: p

      BDIO_write_z64 = -1
      fbd%current => fbd%last
      fbd%rwpos = fbd%last%rend
      p => fbd%current
      If (  (p%rfmt /= BDIO_BIN_F64LE).and.&
           &(p%rfmt /= BDIO_BIN_F64BE).and.&
           &(p%rfmt /= BDIO_BIN_F64  ) ) Then
         Call BDIO_error(fbd,'BDIO_Write_z64', &
              & 'Incorrect data type') 
      End If

      If (  ( (p%rfmt == BDIO_BIN_F64BE).and.(fbd%lendian) ).or. &
           &( (p%rfmt == BDIO_BIN_F64LE).and.(.not.fbd%lendian) ) ) &
           & Then
         CALL ByteSwap(buf)
         Write(fbd%ifn,Pos=fbd%rwpos,Iostat=ios)buf
         CALL ByteSwap(buf)
      Else
         Write(fbd%ifn,Pos=fbd%rwpos,Iostat=ios)buf
      End If
      
      iln = fbd%current%rlen + 16_8*Size(buf)
      CALL Rewrite_rlen(fbd, iln)

      If (ios == 0) BDIO_write_z64 = 16*Size(buf)

      Return
    End Function BDIO_write_z64
    
! ********************************
! *
    Function BDIO_write_z32(fbd,buf)
! *
! ********************************

      Type (BDIO), Intent (inout) :: fbd
      Complex (kind=4), Intent (inout) :: buf(:)
      Integer :: BDIO_write_z32, ios

      Integer (kind=8) :: iln
      Type (BDIO_record), pointer :: p

      BDIO_write_z32 = -1
      fbd%current => fbd%last
      fbd%rwpos = fbd%last%rend
      p => fbd%current
      If (  (p%rfmt /= BDIO_BIN_F32LE).and.&
           &(p%rfmt /= BDIO_BIN_F32BE).and.&
           &(p%rfmt /= BDIO_BIN_F32  ) ) Then
         Call BDIO_error(fbd,'BDIO_Write_z32', &
              & 'Incorrect data type') 
      End If

      If (  ( (p%rfmt == BDIO_BIN_F32BE).and.(fbd%lendian) ).or. &
           &( (p%rfmt == BDIO_BIN_F32LE).and.(.not.fbd%lendian) ) ) &
           & Then
         CALL ByteSwap(buf)
         Write(fbd%ifn,Pos=fbd%rwpos,Iostat=ios)buf
         CALL ByteSwap(buf)
      Else
         Write(fbd%ifn,Pos=fbd%rwpos,Iostat=ios)buf
      End If

      iln = fbd%current%rlen + 8_8*Size(buf)
      CALL Rewrite_rlen(fbd, iln)

      If (ios == 0) BDIO_write_z32 = 8*Size(buf)

      Return
    End Function BDIO_write_z32

! ********************************
! *
    Subroutine BDIO_seek(fbd, nrec)
! *
! ********************************

      Type (BDIO), Intent (inout) :: fbd
      Integer, Intent (in), Optional :: nrec

      Type (BDIO_record), pointer :: p

      If (Present(nrec)) Then
         p => fbd%first
         Do 
            If (p%rid == nrec) Exit
            p => p%next
         End Do
      Else
         p => fbd%current%next
      End If
      
      fbd%current => p
      fbd%rwpos = fbd%current%rpos

      Return
    End Subroutine BDIO_seek

! ********************************
! *
    Subroutine BDIO_start_record(fbd,ifmt,iuinfo,long)
! *
! ********************************

      Type (BDIO), Intent (inout) :: fbd
      Integer, Intent (in) :: ifmt, iuinfo
      Logical, Intent (in), Optional :: long

      Type (BDIO_record), pointer :: newr, aux
      Integer :: irc, i4, ilong
      Integer (kind=8) :: ipos, iln
      logical :: lrec

      fbd%current => fbd%last
      Read(fbd%ifn,Pos=fbd%current%rend)
      ilong = 0
      If (Present(long)) Then
         lrec = long
         if (lrec) ilong = 1
      Else
         lrec = .False.
      End If

      irc = 1
      i4  = 0
      CALL MVbits(irc,0,1,i4,0)
      CALL MVBits(ilong,0,1,i4,3)
      CALL MVBits(ifmt,0,4,i4,4)
      CALL MVBits(iuinfo,0,4,i4,8)
      Write(fbd%ifn)i4
      i4 = 0
      If (lrec) Write(fbd%ifn)i4
      Inquire(fbd%ifn,Pos=ipos)
      fbd%rwpos = ipos

      iln = 0
      Allocate(newr)
      newr%ishdr = .False.
      newr%rlen  = iln
      newr%rpos  = ipos
      newr%rend  = ipos
      newr%rfmt  = ifmt
      newr%ruinfo = iuinfo
      newr%rid   = fbd%tcnt
      newr%islong = lrec
      If (Associated(fbd%first)) Then
         aux => fbd%last

         aux%next => newr
         newr%prev => aux
         fbd%last => newr
      Else
         Allocate(fbd%first)
         fbd%first => newr
         fbd%last  => newr
      End If

      fbd%rcnt = fbd%rcnt + 1
      fbd%tcnt = fbd%tcnt + 1

      Return
    End Subroutine BDIO_start_record

! ********************************
! *
    Function BDIO_read_i32(fbd, ibuf, do_chk)
! *
! ********************************

      Type (BDIO), Intent (inout) :: fbd
      Integer, Intent (inout) :: ibuf(:)
      Logical, Optional :: do_chk
      Integer :: BDIO_read_i32, ios

      Type (BDIO_record), pointer :: p
      Integer (kind=8) :: nmax
      logical :: chk

      BDIO_read_i32 = -1

      If (Present(do_chk)) Then
         chk = do_chk
      Else
         chk = DEFAULT_HASH_CHECK
      End If
      p => fbd%current
      If (  (p%rfmt /= BDIO_BIN_INT32LE).and.&
           &(p%rfmt /= BDIO_BIN_INT32BE).and.&
           &(p%rfmt /= BDIO_BIN_INT32  ) ) Then
         Call BDIO_error(fbd,'BDIO_Read_i32', &
              & 'Incorrect data type') 
      End If

      nmax = size(ibuf)
      If (4_8*size(ibuf) > (p%rend-fbd%rwpos)) nmax = (p%rend-fbd%rwpos)/4
      If (nmax .le. 0) Return
      Read(fbd%ifn,Pos=fbd%rwpos,iostat=ios)ibuf(:nmax)
      If (  ( (p%rfmt == BDIO_BIN_INT32BE).and.(fbd%lendian) ).or. &
           &( (p%rfmt == BDIO_BIN_INT32LE).and.(.not.fbd%lendian) ) ) &
           & CALL ByteSwap(ibuf)

      
      If (Present(do_chk)) p%hash = Hash(ibuf(:nmax), p%hash)
      fbd%rwpos = fbd%rwpos + 4*nmax

      If (ios == 0) BDIO_read_i32 = 4*Int(nmax,kind=4)
      
      Return
    End Function BDIO_read_i32

! ********************************
! *
    Function BDIO_read_i64(fbd, ibuf, do_chk)
! *
! ********************************

      Type (BDIO), Intent (inout) :: fbd
      Integer (kind=8), Intent (inout) :: ibuf(:)
      Logical, Optional :: do_chk
      Integer :: BDIO_read_i64, ios

      Type (BDIO_record), pointer :: p
      Integer (kind=8) :: nmax
      logical :: chk

      BDIO_read_i64 = -1

      If (Present(do_chk)) Then
         chk = do_chk
      Else
         chk = DEFAULT_HASH_CHECK
      End If
      p => fbd%current
      If (  (p%rfmt /= BDIO_BIN_INT64LE).and.&
           &(p%rfmt /= BDIO_BIN_INT64BE).and.&
           &(p%rfmt /= BDIO_BIN_INT64  ) ) Then
         Call BDIO_error(fbd,'BDIO_Read_i64', &
              & 'Incorrect data type') 
      End If

      nmax = size(ibuf)
      If (8_8*size(ibuf) > (p%rend-fbd%rwpos)) nmax = (p%rend-fbd%rwpos)/8
      If (nmax .le. 0) Return
      Read(fbd%ifn,Pos=fbd%rwpos,iostat=ios)ibuf(:nmax)
      If (  ( (p%rfmt == BDIO_BIN_INT64BE).and.(fbd%lendian) ).or. &
           &( (p%rfmt == BDIO_BIN_INT64LE).and.(.not.fbd%lendian) ) ) &
           & CALL ByteSwap(ibuf)
      
      If (Present(do_chk)) p%hash = Hash(ibuf(:nmax), p%hash)
      fbd%rwpos = fbd%rwpos + 8*nmax

      If (ios == 0) BDIO_read_i64 = 8*Int(nmax,kind=4)
      
      Return
    End Function BDIO_read_i64

! ********************************
! *
    Function BDIO_read_f32(fbd, buf, do_chk)
! *
! ********************************

      Type (BDIO), Intent (inout) :: fbd
      Real (kind=4), Intent (inout) :: buf(:)
      Logical, Optional :: do_chk
      Integer :: BDIO_read_f32, ios

      Type (BDIO_record), pointer :: p
      Integer (kind=8) :: nmax
      logical :: chk

      BDIO_read_f32 = -1

      If (Present(do_chk)) Then
         chk = do_chk
      Else
         chk = DEFAULT_HASH_CHECK
      End If
      p => fbd%current
      If (  (p%rfmt /= BDIO_BIN_F32LE).and.&
           &(p%rfmt /= BDIO_BIN_F32BE).and.&
           &(p%rfmt /= BDIO_BIN_F32  ) ) Then
         Call BDIO_error(fbd,'BDIO_Read_f32', &
              & 'Incorrect data type') 
      End If

      nmax = size(buf)
      If (4_8*size(buf) > (p%rend-fbd%rwpos)) nmax = (p%rend-fbd%rwpos)/4
      If (nmax .le. 0) Return
      Read(fbd%ifn,Pos=fbd%rwpos,iostat=ios)buf(:nmax)
      If (  ( (p%rfmt == BDIO_BIN_F32BE).and.(fbd%lendian) ).or. &
           &( (p%rfmt == BDIO_BIN_F32LE).and.(.not.fbd%lendian) ) ) &
           & CALL ByteSwap(buf)

      
      If (Present(do_chk)) p%hash = Hash(buf(:nmax), p%hash)
      fbd%rwpos = fbd%rwpos + 4*nmax

      If (ios == 0) BDIO_read_f32 = 4*Int(nmax,kind=4)
      
      Return
    End Function BDIO_read_f32

! ********************************
! *
    Function BDIO_read_f64(fbd, buf, do_chk)
! *
! ********************************

      Type (BDIO), Intent (inout) :: fbd
      Real (kind=8), Intent (inout) :: buf(:)
      Logical, Optional :: do_chk
      Integer :: BDIO_read_f64, ios

      Type (BDIO_record), pointer :: p
      Integer (kind=8) :: nmax
      logical :: chk

      BDIO_read_f64 = -1

      If (Present(do_chk)) Then
         chk = do_chk
      Else
         chk = DEFAULT_HASH_CHECK
      End If
      p => fbd%current
      If (  (p%rfmt /= BDIO_BIN_F64LE).and.&
           &(p%rfmt /= BDIO_BIN_F64BE).and.&
           &(p%rfmt /= BDIO_BIN_F64  ) ) Then
         Call BDIO_error(fbd,'BDIO_Read_f64', &
              & 'Incorrect data type') 
      End If

      nmax = size(buf)
      If (8_8*size(buf) > (p%rend-fbd%rwpos)) nmax = (p%rend-fbd%rwpos)/8
      If (nmax .le. 0) Return
      Read(fbd%ifn,Pos=fbd%rwpos,iostat=ios)buf(:nmax)
      If (  ( (p%rfmt == BDIO_BIN_F64BE).and.(fbd%lendian) ).or. &
           &( (p%rfmt == BDIO_BIN_F64LE).and.(.not.fbd%lendian) ) ) &
           & CALL ByteSwap(buf)
      
      If (Present(do_chk)) p%hash = Hash(buf(:nmax), p%hash)
      fbd%rwpos = fbd%rwpos + 8*nmax

      If (ios == 0) BDIO_read_f64 = 8*Int(nmax,kind=4)
      
      Return
    End Function BDIO_read_f64

! ********************************
! *
    Function BDIO_read_z32(fbd, buf, do_chk)
! *
! ********************************

      Type (BDIO), Intent (inout) :: fbd
      Complex (kind=4), Intent (inout) :: buf(:)
      Logical, Optional :: do_chk
      Integer :: BDIO_read_z32, ios

      Type (BDIO_record), pointer :: p
      Integer (kind=8) :: nmax
      logical :: chk

      BDIO_read_z32 = -1

      If (Present(do_chk)) Then
         chk = do_chk
      Else
         chk = DEFAULT_HASH_CHECK
      End If
      p => fbd%current
      If (  (p%rfmt /= BDIO_BIN_F32LE).and.&
           &(p%rfmt /= BDIO_BIN_F32BE).and.&
           &(p%rfmt /= BDIO_BIN_F32  ) ) Then
         Call BDIO_error(fbd,'BDIO_Read_f32', &
              & 'Incorrect data type') 
      End If

      nmax = size(buf)
      If (8_8*size(buf) > (p%rend-fbd%rwpos)) nmax = (p%rend-fbd%rwpos)/8
      If (nmax .le. 0) Return
      Read(fbd%ifn,Pos=fbd%rwpos,iostat=ios)buf(:nmax)
      If (  ( (p%rfmt == BDIO_BIN_F32BE).and.(fbd%lendian) ).or. &
           &( (p%rfmt == BDIO_BIN_F32LE).and.(.not.fbd%lendian) ) ) &
           & CALL ByteSwap(buf)

      
      If (Present(do_chk)) p%hash = Hash(buf(:nmax), p%hash)
      fbd%rwpos = fbd%rwpos + 8*nmax

      If (ios == 0) BDIO_read_z32 = 8*Int(nmax,kind=4)
      
      Return
    End Function BDIO_read_z32

! ********************************
! *
    Function BDIO_read_z64(fbd, buf, do_chk)
! *
! ********************************

      Type (BDIO), Intent (inout) :: fbd
      Complex (kind=8), Intent (inout) :: buf(:)
      Logical, Optional :: do_chk
      Integer :: BDIO_read_z64, ios

      Type (BDIO_record), pointer :: p
      Integer (kind=8) :: nmax
      logical :: chk

      BDIO_read_z64 = -1

      If (Present(do_chk)) Then
         chk = do_chk
      Else
         chk = DEFAULT_HASH_CHECK
      End If
      p => fbd%current
      If (  (p%rfmt /= BDIO_BIN_F64LE).and.&
           &(p%rfmt /= BDIO_BIN_F64BE).and.&
           &(p%rfmt /= BDIO_BIN_F64  ) ) Then
         Call BDIO_error(fbd,'BDIO_Read_f64', &
              & 'Incorrect data type') 
      End If

      nmax = size(buf)
      If (16_8*size(buf) > (p%rend-fbd%rwpos)) nmax = (p%rend-fbd%rwpos)/16
      If (nmax .le. 0) Return
      Read(fbd%ifn,Pos=fbd%rwpos,iostat=ios)buf(:nmax)
      If (  ( (p%rfmt == BDIO_BIN_F64BE).and.(fbd%lendian) ).or. &
           &( (p%rfmt == BDIO_BIN_F64LE).and.(.not.fbd%lendian) ) ) &
           & CALL ByteSwap(buf)
      
      If (Present(do_chk)) p%hash = Hash(buf(:nmax), p%hash)
      fbd%rwpos = fbd%rwpos + 16*nmax

      If (ios == 0) BDIO_read_z64 = 16*Int(nmax,kind=4)
      
      Return
    End Function BDIO_read_z64
    
! ********************************
! *
    Function BDIO_open(fname, mode, protocol_info) Result (fbd)
! *
! ********************************
      Character (len=*), Intent (in) :: fname
      Character (len=1), Intent (in) :: mode
      Character (len=*), Intent (in), Optional ::  protocol_info
      Integer (kind=8) :: ipos
      
      Type (BDIO) :: fbd

      Logical :: is_used

      Select Case (mode)
      Case ('r')
         fbd%imode = BDIO_R_MODE
      Case ('w')
         fbd%imode = BDIO_W_MODE
      Case ('a')
         fbd%imode = BDIO_A_MODE
      Case Default
         Call BDIO_error(fbd,'BDIO_Open', &
              & 'Incorrect mode') 
      End Select
      
      fbd%lendian = isLittleEndian()
      fbd%ifn = 69
      Do 
         Inquire(fbd%ifn, Opened=is_used)
         If (is_used) Then
            fbd%ifn = fbd%ifn + 1
         Else
            Exit
         End If
      End Do

      Select Case (fbd%imode)
      Case (BDIO_R_MODE) 
         Open (File=Trim(fname), unit=fbd%ifn, ACTION="READ", &
              & Form='UNFORMATTED', Access='STREAM')
         Inquire(fbd%ifn, Pos=ipos)

         fbd%opened = .True.
!         CALL BDIO_Read_header(fbd)
      CASE (BDIO_A_MODE) 
         Open (File=Trim(fname), unit=fbd%ifn, ACTION="READWRITE", &
              & Form='UNFORMATTED', Access='STREAM')
         Rewind(fbd%ifn)

         fbd%opened = .True.
      End Select

      CALL BDIO_parse(fbd)
      fbd%current => fbd%first

      Return
    End Function BDIO_open

! ********************************
! *
    Subroutine BDIO_parse(fbd)
! *
! ********************************

      Type (BDIO), Intent (inout) :: fbd

      Integer (kind=4) :: i4, irc
      Integer (kind=8) :: ipos
      
      
      Rewind(fbd%ifn)
      Do 
         Inquire(fbd%ifn, Pos=ipos)
         Read(fbd%ifn,END=20)i4
         irc = 0
         CALL MVbits(i4,0,1,irc,0)

         If (i4==BDIO_MAGIC) Then
            Read(fbd%ifn, Pos=ipos)
            CALL BDIO_Read_header(fbd)
         Else if (irc == 1) Then
            Read(fbd%ifn, Pos=ipos)
            CALL BDIO_Read_record(fbd)
         Else
            Write(0,*)'Corrupt file!'
         End If
      End Do
         
20    Continue

      Return
    End Subroutine BDIO_parse

! ********************************
! *
    Subroutine BDIO_read_record(fbd)
! *
! ********************************

      Type (BDIO), Intent (inout) :: fbd

      Integer (kind=4) :: i4, ilong, irc, j, ifmt, iuinfo
      Integer (kind=8) :: ipos, iend, iln
      Integer (kind=8) :: jlong
      Character (kind=1) :: ch

      Type (BDIO_Record), pointer :: newr, aux

      ! Record
      Read(fbd%ifn,END=20)i4
      irc = 0
      ilong = 0
      ifmt = 0
      iuinfo = 0
      CALL MVbits(i4,0,1,irc,0)
      CALL MVbits(i4,3,1,ilong,0)
      CALL MVbits(i4,4,4,ifmt,0)
      CALL MVbits(i4,8,4,iuinfo,0)
      
      j = 0
      CALL MVbits(i4,12,20,j,0)
      iln = Int(j,kind=8)
      If (ilong==1) Then
         Read(fbd%ifn)i4
         CALL MVBits(i4,0,32,j,0)
         jlong=Int(j,kind=8)
         CALL MVBits(jlong,20,32,iln,0)
      End If
      Inquire(fbd%ifn, Pos=ipos)
      Read(fbd%ifn)(ch, jlong=1, iln)
      Inquire(fbd%ifn, Pos=iend)
      
      Allocate(newr)
      newr%ishdr = .False.
      newr%rlen  = iln
      newr%rpos  = ipos
      newr%rend  = iend
      newr%rfmt  = ifmt
      newr%ruinfo = iuinfo
      newr%rid   = fbd%tcnt
      newr%islong = .False. 
      If (ilong==1) Then
         newr%islong = .True. 
      End If
      
      If (Associated(fbd%first)) Then
         aux => fbd%last

         aux%next => newr
         newr%prev => aux
         fbd%last => newr
      Else
         Allocate(fbd%first)
         fbd%first => newr
         fbd%last  => newr
      End If
      fbd%rcnt = fbd%rcnt + 1
      fbd%tcnt = fbd%tcnt + 1
         
20    Continue

      Return
    End Subroutine BDIO_read_record

! ********************************
! *
    Subroutine BDIO_read_header(ptf)
! *
! ********************************

      Type (BDIO), Intent (inout) :: ptf

      Integer (kind=4) :: i4, j, isp, iv, iln
      Character (kind=1) :: ch
      Integer (kind=8) :: ipos, iend

      Type (BDIO_Record), pointer :: newr, aux

      Read(ptf%ifn,END=20)i4
      if (.not.ptf%lendian) CALL ByteSwap(i4)
      If (i4/=BDIO_MAGIC) Then
         Write(0,*)'File not recognized as BDIO'
         Stop
      End If

      Read(ptf%ifn)i4
      if (.not.ptf%lendian) CALL ByteSwap(i4)

      iln = 0
      isp = 0
      iv  = 0
      CALL MVbits(i4,0,12,iln,0)
      CALL MVbits(i4,12,4,isp,0)
      CALL MVbits(i4,16,16,iv,0)
      Inquire(ptf%ifn, pos=ipos)
      Read(ptf%ifn)(ch, j=1, iln-1)
      Inquire(ptf%ifn, pos=iend)
      Read(ptf%ifn)ch
      
      Allocate(newr)
      newr%ishdr = .True.
      newr%rlen  = Int(iln,kind=8)
      newr%rpos  = ipos
      newr%rend  = iend
      newr%rid   = ptf%tcnt
      newr%iver  = iv
      If (Associated(ptf%first)) Then
         aux => ptf%last

         aux%next => newr
         newr%prev => aux
         ptf%last => newr
      Else
         Allocate(ptf%first)
         ptf%first => newr
         ptf%last  => newr
      End If
      ptf%hcnt = ptf%hcnt + 1
      ptf%tcnt = ptf%tcnt + 1

20    Continue

      Return
    End Subroutine BDIO_read_header

! ********************************
! *
    Subroutine BDIO_show_header(p)
! *
! ********************************
      
      Type (BDIO_record), pointer :: p
      
      Write(*,'(1I4,3X,1A,3X,1I4,1I10,1I4,1A15,1L4)')&
           & p%rid, 'header',p%iver, p%rlen, p%ruinfo, '', p%islong  

      Return
    End Subroutine BDIO_show_header

! ********************************
! *
    Subroutine BDIO_show_record(p)
! *
! ********************************
      
      Type (BDIO_record), pointer :: p
      
      Write(*,'(1I4,3X,1A,3X,1I4,1I10,1I6,1A15,1L4)')&
           & p%rid, 'record', p%rfmt, p%rlen, p%ruinfo, '', p%islong  

      Return
    End Subroutine BDIO_show_record

! ********************************
! *
    Subroutine BDIO_show(ptf)
! *
! ********************************

      Type (BDIO), Intent (in) :: ptf
      Type (BDIO_record), pointer :: p
      
      Write(*,'(1A,5X,1A,5X,1A,5X,1A,3X,1A,2X,1A)')'ID','record', 'type', 'size', 'uID', 'starts with', 'long' 

      p => ptf%first
      Do 
         If (.not.Associated(p)) Exit
         If (p%ishdr) Then
            CALL BDIO_show_header(p)
         Else
            CALL BDIO_show_record(p)
         End If
         p => p%next
      End Do

!!$      p => ptf%last
!!$      Do 
!!$         If (.not.Associated(p)) Exit
!!$         If (p%ishdr) Then
!!$            Write(*,*)'HDR', p%rid, p%rlen, p%rpos, p%rend
!!$         Else
!!$            Write(*,*)'RCD', p%rid, p%rlen, p%rpos, p%rend
!!$         End If
!!$         p => p%prev
!!$      End Do


      Return
    End Subroutine BDIO_show

! ********************************
! *
    Subroutine Rewrite_rlen(fbd, iln)
! *
! ********************************
      
      Type (BDIO), Intent (inout) :: fbd
      Integer (kind=8), Intent (in) :: iln
      
      Integer (kind=4) :: i4(2), j
      Integer (kind=8) :: jl
      Type (BDIO_record), pointer :: p

      p => fbd%current
      If (p%islong) Then
         CALL MVBits(iln, 0, 20, jl, 0)
         j = Int(jl,kind=4)
         CALL MVBits(j,0,20,i4(1),12)
         CALL MVBits(iln, 21, 32, jl, 0)
         j = Int(jl,kind=4)
         CALL MVBits(j,0,32,i4(2),0)
         Write(fbd%ifn,Pos=p%rpos-8)i4(1:2)
      Else
         Read(fbd%ifn,Pos=p%rpos-4)i4(1)
         CALL MVBits(iln, 0, 20, jl, 0)
         j = Int(jl,kind=4)
         CALL MVBits(j,0,20,i4(1),12)
         Write(fbd%ifn,Pos=fbd%current%rpos-4)i4(1)
      End If

      p%rlen = iln
      p%rend = fbd%current%rpos + iln

      Return
    End Subroutine Rewrite_rlen

! ********************************         
! *  
    Function isLittleEndian()
! * 
! ********************************    

      Logical :: isLittleEndian
      Character (kind=1) :: ch(4)
      Integer (kind=4) :: I(1), Isw

      I   = 0
      Isw = 0
      ch = Transfer(I, ch)
      ch(1) = achar(48)
      I   = Transfer(ch, I)

      Isw = 48
      CALL ByteSwap(Isw)

      isLittleEndian = .True.
      If (I(1) == 48) Then
         Return
      Else If (I(1) == Isw) Then
         isLittleEndian = .False.
         Return
      Else
         Write(0,*)'What a mess of Endianness!'
      End If

      Return
    End Function isLittleEndian

! ********************************
! *
    Subroutine byteswap_int32(Iii)
! *
! ********************************

      Integer (kind=4), Intent (inout) :: Iii
      
      Character (kind=1) :: isw(4)
      Integer (kind=4) :: J(1)

      isw  = Transfer(Iii,isw, 4)
      J    = Transfer(isw(4:1:-1), J)
      Iii  = J(1)

      Return
    End Subroutine byteswap_int32

! ********************************
! *
    Subroutine byteswap_int32V(Iii)
! *
! ********************************

      Integer (kind=4), Intent (inout) :: Iii(:)
      
      Character (kind=1) :: isw(4)
      Integer (kind=4) :: J(1), K

      Do K = 1, Size(Iii)
         isw  = Transfer(Iii(K),isw, 4)
         J    = Transfer(isw(4:1:-1), J)
         Iii(K)  = J(1)
      End Do

      Return
    End Subroutine byteswap_int32V

! ********************************
! *
    Subroutine byteswap_int64(Iii)
! *
! ********************************

      Integer (kind=8), Intent (inout) :: Iii
      
      Character (kind=1) :: isw(8)
      Integer (kind=8) :: J(1)

      isw  = Transfer(Iii,isw, 8)
      J    = Transfer(isw(8:1:-1), J)
      Iii  = J(1)

      Return
    End Subroutine byteswap_int64

! ********************************
! *
    Subroutine byteswap_int64V(Iii)
! *
! ********************************

      Integer (kind=8), Intent (inout) :: Iii(:)
      
      Character (kind=1) :: isw(8)
      Integer (kind=8) :: J(1)
      Integer :: K

      Do K = 1, Size(Iii)
         isw  = Transfer(Iii(K),isw, 8)
         J    = Transfer(isw(8:1:-1), J)
         Iii(K)  = J(1)
      End Do

      Return
    End Subroutine byteswap_int64V

! ********************************
! *
    Subroutine byteswap_R(d)
! *
! ********************************

      Real (kind=4), Intent (inout) :: d
      
      Character (kind=1) :: isw(4)
      Real (kind=4) :: dd(1)

      isw  = Transfer(d,isw, 4)
      dd   = Transfer(isw(4:1:-1), dd)
      d  = dd(1)

      Return
    End Subroutine byteswap_R

! ********************************
! *
    Subroutine byteswap_RV(d)
! *
! ********************************

      Real (kind=4), Intent (inout) :: d(:)
      
      Character (kind=1) :: isw(4)
      Real (kind=4) :: dd(1)
      Integer :: I

      Do I = 1, Size(d)
         isw  = Transfer(d(I),isw, 4)
         dd   = Transfer(isw(4:1:-1), dd)
         d(I) = dd(1)
      End Do

      Return
    End Subroutine byteswap_RV

! ********************************
! *
    Subroutine byteswap_D(d)
! *
! ********************************

      Real (kind=8), Intent (inout) :: d
      
      Character (kind=1) :: isw(8)
      Real (kind=8) :: dd(1)

      isw  = Transfer(d,isw, 8)
      dd   = Transfer(isw(8:1:-1), dd)
      d  = dd(1)

      Return
    End Subroutine byteswap_D

! ********************************
! *
    Subroutine byteswap_DV(d)
! *
! ********************************

      Real (kind=8), Intent (inout) :: d(:)
      
      Character (kind=1) :: isw(8)
      Real (kind=8) :: dd(1)
      Integer :: I

      Do I = 1, Size(d)
         isw  = Transfer(d(I),isw, 8)
         dd   = Transfer(isw(8:1:-1), dd)
         d(I) = dd(1)
      End Do

      Return
    End Subroutine byteswap_DV

! ********************************
! *
    Subroutine byteswap_Z(z)
! *
! ********************************

      Complex (kind=4), Intent (inout) :: z
      
      Real (kind=4) :: d(2)

      d = Transfer(z,d)
      CALL ByteSwap(d)
      z = Transfer(d, z)

      Return
    End Subroutine byteswap_Z

! ********************************
! *
    Subroutine byteswap_DZ(z)
! *
! ********************************

      Complex (kind=8), Intent (inout) :: z
      
      Real (kind=8) :: d(2)

      d = Transfer(z,d)
      CALL ByteSwap(d)
      z = Transfer(d, z)

      Return
    End Subroutine byteswap_DZ

! ********************************
! *
    Subroutine byteswap_DZV(z)
! *
! ********************************

      Complex (kind=8), Intent (inout) :: z(:)
      
      Real (kind=8) :: d(2)
      Integer :: I

      Do I = 1, Size(z)
         d = Transfer(z(I),d)
         CALL ByteSwap(d)
         z(I) = Transfer(d, z(I))
      End Do

      Return
    End Subroutine byteswap_DZV

! ********************************
! *
    Subroutine byteswap_ZV(z)
! *
! ********************************

      Complex (kind=4), Intent (inout) :: z(:)
      
      Real (kind=4) :: d(2)
      Integer :: I

      Do I = 1, Size(z)
         d = Transfer(z(I),d)
         CALL ByteSwap(d)
         z(I) = Transfer(d, z(I))
      End Do

      Return
    End Subroutine byteswap_ZV


End MODULE ModBDio
