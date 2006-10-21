!
! MODULE with useful constants
!
! Copyright (C) 2003  Alberto Ramos <alberto@martin.ft.uam.es>
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

! ***********************************************************
! *
MODULE Constants
! *
! *********************************************************** 

  USE NumTypes

  ! Mathematical constants related to pi

  Complex (Kind=DPC), Parameter :: UNITIMAG_DPC = (0.0_DP, 1.0_DP)
  Complex (Kind=DPC), Parameter :: PI_IMAG_DPC = (0.0_DP, 3.141592653589793238462_DP)
  Complex (Kind=DPC), Parameter :: TWOPI_IMAG_DPC = (0.0_DP, 6.283185307179586476925_DP)
  Complex (Kind=DPC), Parameter :: HALFPI_IMAG_DPC = (0.0_DP, 1.57079632679489661923_DP)

  Complex (Kind=SPC), Parameter :: UNITIMAG_SPC = (0.0_SP, 1.0_SP)
  Complex (Kind=SPC), Parameter :: PI_IMAG_SPC = (0.0_SP, 3.141592653589793238462_SP)
  Complex (Kind=SPC), Parameter :: TWOPI_IMAG_SPC = (0.0_SP, 6.283185307179586476925_SP)
  Complex (Kind=SPC), Parameter :: HALFPI_IMAG_SPC = (0.0_SP, 1.57079632679489661923_SP)

  Real (kind=DP), Parameter :: PI_DP = 3.141592653589793238462_DP
  Real (kind=DP), Parameter :: TWOPI_DP = 6.283185307179586476925_DP
  Real (kind=DP), Parameter :: HALFPI_DP = 1.57079632679489661923_DP

  Real (kind=SP), Parameter :: PI_SP = 3.141592653589793238462_SP
  Real (kind=SP), Parameter :: TWOPI_SP = 6.283185307179586476925_SP
  Real (kind=SP), Parameter :: HALFPI_SP = 1.57079632679489661923_SP

  
  ! Mathematical constants related with roots and log
  
  Real (kind=DP), Parameter :: SR2_DP  = 1.4142135623730950488_DP
  Real (kind=DP), Parameter :: SR3_DP  = 1.7320508075688772935_DP
  Real (kind=DP), Parameter :: SR5_DP  = 2.2360679774997896964_DP
  Real (kind=DP), Parameter :: SRe_DP  = 1.6487212707001281468_DP
  Real (kind=DP), Parameter :: SRpi_DP = 1.772453850905516027298167_DP

  Real (kind=DP), Parameter :: LG102_DP  = 0.3010299956639811952137389_DP
  Real (kind=DP), Parameter :: LG103_DP  = 0.4771212547196624372950279_DP
  Real (kind=DP), Parameter :: LG10e_DP  = 0.43429448190325182765_DP
  Real (kind=DP), Parameter :: LG10pi_DP = 0.4971498726941338543512683_DP

  Real (kind=DP), Parameter :: LGe2_DP  = 0.693147180559945309417232_DP
  Real (kind=DP), Parameter :: LGe3_DP  = 1.098612288668109691395245_DP
  Real (kind=DP), Parameter :: LGe10_DP = 2.302585092994045684017991_DP


  Real (kind=SP), Parameter :: SR2_SP  = 1.4142135623730950488_SP
  Real (kind=SP), Parameter :: SR3_SP  = 1.7320508075688772935_SP
  Real (kind=SP), Parameter :: SR5_SP  = 2.2360679774997896964_SP
  Real (kind=SP), Parameter :: SRe_SP  = 1.6487212707001281468_SP
  Real (kind=SP), Parameter :: SRpi_SP = 1.772453850905516027298167_SP

  Real (kind=SP), Parameter :: LG102_SP  = 0.3010299956639811952137389_SP
  Real (kind=SP), Parameter :: LG103_SP  = 0.4771212547196624372950279_SP
  Real (kind=SP), Parameter :: LG10e_SP  = 0.43429448190325182765_SP
  Real (kind=SP), Parameter :: LG10pi_SP = 0.4971498726941338543512683_SP

  Real (kind=SP), Parameter :: LGe2_SP  = 0.693147180559945309417232_SP
  Real (kind=SP), Parameter :: LGe3_SP  = 1.098612288668109691395245_SP
  Real (kind=SP), Parameter :: LGe10_SP = 2.302585092994045684017991_SP


  ! Other mathematical constants
  
  Real (kind=DP), Parameter :: GEULER_DP = 0.577215664901532860606512_DP

  Real (kind=SP), Parameter :: GEULER_SP = 0.577215664901532860606512_SP



End MODULE Constants
