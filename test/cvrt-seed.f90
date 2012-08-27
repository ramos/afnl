
Program Seed

  USE NumTypes
  USE Statistics
  USE Time

  Real (kind=DP), Allocatable :: X(:)
  Character (len=5000) :: arg
  Integer :: is

  CALL Get_command_argument(1, arg)
  
  CALL ReadLuxSeedOld(Trim(arg))
  CALL WriteLuxSeed(Trim(arg))

  Stop
End Program Seed
