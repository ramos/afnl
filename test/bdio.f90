
Program AA

  USE NumTypes
  USE BDIO

  Type (BDIO_t) :: bb
  Character (len=200) :: fn
  Character (kind=1) :: ch
  Integer :: I(1), Isw
  Integer :: Ia = 805306416!805306368
  Real (kind=4) :: r4
  Real (kind=8) :: r8
  Complex (kind=DPC) :: z16


!  If (.not.GetOpt('-i', fn)) Stop
!  bb = BDIO_open(fn, 'r', 'Perico de los palotes')
  
!  CALL BDIO_Read_header(bb)

  If (isLittleEndian()) Then
     Write(*,*)'Machine is Little Endian'
  Else
     Write(*,*)'Machine is Big Endian'
  End If
  
  Call Get_command_argument(1, fn)
  If (fn == '-w') Then
     If (isLittleEndian()) Then
        Open(File='int_one.le', unit=66, ACTION="Write", &
             & Form='UNFORMATTED', Access='STREAM')
        Isw = 1
        Write(66)Isw
        Close(66)
        
        Open(File='r4_one.le', unit=66, ACTION="Write", &
             & Form='UNFORMATTED', Access='STREAM')
        r4 = 1.0
        Write(66)r4
        Close(66)
        
        Open(File='r8_one.le', unit=66, ACTION="Write", &
             & Form='UNFORMATTED', Access='STREAM')
        r8 = 1.0D0
        Write(66)r8
        Close(66)
        
        Open(File='int_one.be', unit=66, ACTION="Write", &
             & Form='UNFORMATTED', Access='STREAM')
        Isw = 1
        CALL ByteSwap(Isw)
        Write(66)Isw
        Close(66)
        
        Open(File='r4_one.be', unit=66, ACTION="Write", &
             & Form='UNFORMATTED', Access='STREAM')
        r4 = 1.0
        CALL ByteSwap(r4)
        Write(66)r4
        Close(66)
        
        Open(File='r8_one.be', unit=66, ACTION="Write", &
             & Form='UNFORMATTED', Access='STREAM')
        r8 = 1.0D0
        CALL ByteSwap(r8)
        Write(66)r8
        Close(66)
     Else
        Open(File='int_one.le', unit=66, ACTION="Write", &
             & Form='UNFORMATTED', Access='STREAM')
        Isw = 1
        CALL ByteSwap(Isw)
        Write(66)Isw
        Close(66)
        
        Open(File='r4_one.le', unit=66, ACTION="Write", &
             & Form='UNFORMATTED', Access='STREAM')
        r4 = 1.0
        CALL ByteSwap(r4)
        Write(66)r4
        Close(66)
        
        Open(File='r8_one.le', unit=66, ACTION="Write", &
             & Form='UNFORMATTED', Access='STREAM')
        r8 = 1.0D0
        CALL ByteSwap(r8)
        Write(66)r8
        Close(66)
        
        Open(File='int_one.be', unit=66, ACTION="Write", &
             & Form='UNFORMATTED', Access='STREAM')
        Isw = 1
        Write(66)Isw
        Close(66)
        
        Open(File='r4_one.be', unit=66, ACTION="Write", &
             & Form='UNFORMATTED', Access='STREAM')
        r4 = 1.0
        Write(66)r4
        Close(66)
        
        Open(File='r8_one.be', unit=66, ACTION="Write", &
             & Form='UNFORMATTED', Access='STREAM')
        r8 = 1.0D0
        Write(66)r8
        Close(66)
     End If
  Else
     If (isLittleEndian()) Then
        Open(File='int_one.le', unit=66, ACTION="READ", &
             & Form='UNFORMATTED', Access='STREAM')
        Read(66)Isw
        Close(66)
        
        Open(File='r4_one.le', unit=66, ACTION="READ", &
             & Form='UNFORMATTED', Access='STREAM')
        Read(66)r4
        Close(66)
        
        Open(File='r8_one.le', unit=66, ACTION="READ", &
             & Form='UNFORMATTED', Access='STREAM')
        Read(66)r8
        Close(66)
        
        Write(*,*)'File int_one.le'
        Write(*,*)Isw
        Write(*,*)'File r4_one.le'
        Write(*,*)r4
        Write(*,*)'File r8_one.le'
        Write(*,*)r8
        
        Open(File='int_one.be', unit=66, ACTION="READ", &
             & Form='UNFORMATTED', Access='STREAM')
        Read(66)Isw
        CALL ByteSwap(Isw)
        Close(66)
        
        Open(File='r4_one.be', unit=66, ACTION="READ", &
             & Form='UNFORMATTED', Access='STREAM')
        Read(66)r4
        CALL ByteSwap(r4)
        Close(66)
        
        Open(File='r8_one.be', unit=66, ACTION="READ", &
             & Form='UNFORMATTED', Access='STREAM')
        Read(66)r8
        CALL ByteSwap(r8)
        Close(66)

        Write(*,*)'File int_one.be'
        Write(*,*)Isw
        Write(*,*)'File r4_one.be'
        Write(*,*)r4
        Write(*,*)'File r8_one.be'
        Write(*,*)r8
     Else
        Open(File='int_one.le', unit=66, ACTION="READ", &
             & Form='UNFORMATTED', Access='STREAM')
        Read(66)Isw
        CALL ByteSwap(Isw)
        Close(66)
        
        Open(File='r4_one.le', unit=66, ACTION="READ", &
             & Form='UNFORMATTED', Access='STREAM')
        Read(66)r4
        CALL ByteSwap(r4)
        Close(66)
        
        Open(File='r8_one.le', unit=66, ACTION="READ", &
             & Form='UNFORMATTED', Access='STREAM')
        Read(66)r8
        CALL ByteSwap(r8)
        Close(66)

        Write(*,*)'File int_one.le'
        Write(*,*)Isw
        Write(*,*)'File r4_one.le'
        Write(*,*)r4
        Write(*,*)'File r8_one.le'
        Write(*,*)r8
        
        Open(File='int_one.be', unit=66, ACTION="READ", &
             & Form='UNFORMATTED', Access='STREAM')
        Read(66)Isw
        Close(66)
        
        Open(File='r4_one.be', unit=66, ACTION="READ", &
             & Form='UNFORMATTED', Access='STREAM')
        Read(66)r4
        Close(66)
        
        Open(File='r8_one.be', unit=66, ACTION="READ", &
             & Form='UNFORMATTED', Access='STREAM')
        Read(66)r8
        Close(66)

        Write(*,*)'File int_one.be'
        Write(*,*)Isw
        Write(*,*)'File r4_one.be'
        Write(*,*)r4
        Write(*,*)'File r8_one.be'
        Write(*,*)r8
     End If
  End If

  z16 = Cmplx(1.0D0, 11.0D0)
  
  Open(File='z16.le', unit=66, ACTION="Write", &
       & Form='UNFORMATTED', Access='STREAM')
  Write(66)z16
  Close(66)
  
  Open(File='z16.be', unit=66, ACTION="Write", &
       & Form='UNFORMATTED', Access='STREAM')
  Write(*,*)z16
  CALL ByteSwap(z16)
  Write(*,*)z16
  Write(66)z16
  Close(66)
  
  Write(*,*)'sP', sP
  Write(*,*)'DP', DP
  Write(*,*)'sPC', sPC
  Write(*,*)'DPC', DPC



  Stop
End Program AA
