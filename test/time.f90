
Program TestTime

  USE NumTypes
  USE Time

  Type (tm) :: OneDay, moment

!  Write(*,*)Julian_Date(gettime())

  OneDay%hour = 12
  OneDay%min = 0
  OneDay%sec = 0
  OneDay%msec = 0
  OneDay%mday = 5
  OneDay%mon = 0
  OneDay%year = 2008
  OneDay%wday = 3

!  OneDay = gettime()
  Write(*,*)'Year: ', OneDay%year
  Write(*,*)'Year: ', OneDay%mday
  Write(*,*)'Year: ', OneDay%mon
  Write(*,*)'Year: ', OneDay%wday
  Write(*,*)'Year: ', Modulo(OneDay%wday,7)
  Write(*,*)'Year: ', Mod(OneDay%wday,7)
  Write(*,*)asctime(OneDay)
  Write(*,*)Len(asctime(gettime()))
  Write(*,*)asctime(gettime())
  Write(*,*)Julian_Date(OneDay)


  moment = gettime()
  Write(*,'(3I2.2,1I4.4)')moment%hour, moment%min, moment%sec,&
          & moment%msec
  Write(*,*)'Julian date: ', Julian_Date(moment)


  Stop
End Program TestTime
