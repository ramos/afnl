
Program TestTime

  USE NumTypes
  USE Time

  Type (tm) :: OneDay, moment

  Write(*,*)Julian_Date(gettime())

  OneDay%hour = 12
  OneDay%min = 0
  OneDay%sec = 0
  OneDay%msec = 0
  OneDay%mday = 31
  OneDay%mon = 7
  OneDay%year = 2132
  OneDay%wday = 3

  Write(*,*)asctime(OneDay)
  Write(*,*)Julian_Date(OneDay)


  moment = gettime()
  Write(*,'(3I2.2,1I4.4)')moment%hour, moment%min, moment%sec,&
          & moment%msec

  Stop
End Program TestTime
