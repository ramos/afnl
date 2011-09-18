
Program TestTime

  USE NumTypes
  USE Time

  Type (tm) :: OneDay

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


  Stop
End Program TestTime
