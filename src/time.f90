!
! MODULE to deal with date and time
!
! Copyright (C) 2006  Alberto Ramos <alberto@martin.ft.uam.es>
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

! $ v. 1.0; Released: 10/01/2007; $

MODULE Time

  USE NumTypes
  USE Error

  IMPLICIT NONE

  ! The year the Gregorian calendar was implemented 
  ! for more details take a look to
  ! http://en.wikipedia.org/wiki/Gregorian_calendar
  Integer, Parameter :: GREGORIAN_DATE = 1582

  Character (len=3), Dimension  (0:6), Parameter :: NAMES_OF_DAYS = &
       & (/'Sun', 'Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat'/)

  Character (len=3), Dimension (0:11), Parameter :: NAMES_OF_MONTHS = &
       & (/'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', &
       &   'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'  /)

  Type tm
     Integer :: hour  ! Hour of the day [0-23]
     Integer :: min   ! Minutes after the hour [0-59]
     Integer :: sec   ! Seconds after the minute [0-59]
     Integer :: msec  ! Miliseconds after the second [0-999]
     Integer :: mday  ! Day of the month [1-31]
     Integer :: mon   ! Month of the year [0-11]
     Integer :: year  ! Year 
     Integer :: wday  ! Day of the week since sunday [0-6]
  End Type tm
  

  Private NAMES_OF_DAYS, NAMES_OF_MONTHS, GREGORIAN_DATE

CONTAINS

! ******************************************************
! *
  Integer Function Day_of_Week(Day, Month, Year)
! *
! ******************************************************
! * Computes the day of the week using the Zeller's
! * congruence. More details at:
! * http://en.wikipedia.org/wiki/Zeller%27s_congruence
! ******************************************************
    Integer, Intent (in) :: Day, Month, Year
    Integer :: Mon, yy, cc, cyear
    
    If (Month == 0) Then
       Mon = 13
       cyear = year - 1
    Else If (Month == 1) Then
       Mon = 14
       cyear = year - 1
    Else 
       Mon = Month + 1
       cyear = year
    End If
    yy = Mod(cyear, 100)
    cc = Int((cyear - yy)/100)

    If (Year < GREGORIAN_DATE) Then
       Day_Of_Week = Day + Int(26*(mon+1)/10) + Int(5*yy/4) - &
            & cc + 5
    Else
       Day_Of_Week = Day + Int(26*(mon+1)/10) + Int(5*yy/4) + &
            & Int(cc/4) - 2*cc 
    End If
    Day_Of_Week = Mod(Day_Of_Week-1,7)


    Return
  End Function Day_of_Week

! ***************************************
! *
  Type (tm) Function gettime()
! *
! ***************************************
    
    Character (len=8)  :: Datestr
    Character (len=10) :: Timestr
    Character (len=1) :: ch

    CALL Date_and_Time(Datestr, Timestr)

    Read(Datestr,'(1I4,1I2,1I2)')Gettime%year, &
         & Gettime%mon, Gettime%mday
    Read(Timestr,'(1I2,1I2,1I2,1A1,1I3)')Gettime%hour, &
         & Gettime%min, Gettime%sec, ch, Gettime%msec

    Gettime%mon = Gettime%mon - 1
    
    Gettime%wday = Day_of_week(Gettime%mday, Gettime%mon, Gettime%year)

    Return
  End Function gettime

! ***************************************
! *
  Character (len=24) Function asctime(t)
! *
! ***************************************
    
    Type (tm), Intent (in) :: t
  
10  FORMAT((4A,1I2,1A,1I2,1A1,1I2.2,1A1,1I2.2,1A1,1I4))
    
    Write(asctime,10)NAMES_OF_DAYS(t%wday), ' ', &
         & NAMES_OF_MONTHS(t%mon), ' ', t%mday,' ', &
         & t%hour, ':', t%min, ':', t%sec, ' ', t%year
    
    Return
  End Function asctime

  


! ***************************************
! *
  Logical Function isleap(Nyear)
! *
! ***************************************
! * Returns .true. if the year is leap,
! * .false. in any other case
! ***************************************
    
    Integer, Intent (in) :: Nyear

    If (Nyear < GREGORIAN_DATE) Then
       If (Mod(Nyear,4) == 0) Then
          isleap = .True.
       Else
          isleap = .False.
       End If
    Else
       If ( (Mod(Nyear, 400) == 0) .or. &
            & ((Mod(Nyear,4) == 0) .and. (Mod(Nyear,100) /= 0))) Then
          isleap = .True.
       Else
          isleap = .False.
       End If
    End If
          
    Return
  End Function isleap


End MODULE Time
