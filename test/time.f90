
Program TestTime

  USE NumTypes
  USE Time

  Type (tm) :: OneDay, moment

!  Write(*,*)Julian_Date(gettime())

  OneDay%hour = 12
  OneDay%min = 0
  OneDay%sec = 0
  OneDay%msec = 0
  OneDay%mday = 31
  OneDay%mon = 7
  OneDay%year = 2132
  OneDay%wday = 3

  OneDay = gettime()
  Write(*,*)'Year: ', OneDay%year
  Write(*,*)'Year: ', OneDay%mday
  Write(*,*)'Year: ', OneDay%mon
  Write(*,*)'Year: ', OneDay%wday
  Write(*,*)'Year: ', Modulo(OneDay%wday,7)
  Write(*,*)'Year: ', Mod(OneDay%wday,7)
  Write(*,*)asctime(OneDay)
  Stop
  Write(*,*)Len(asctime(gettime()))
  Write(*,*)asctime(gettime())
  Write(*,*)Julian_Date(OneDay)


  moment = gettime()
  Write(*,'(3I2.2,1I4.4)')moment%hour, moment%min, moment%sec,&
          & moment%msec

  Stop
End Program TestTime
