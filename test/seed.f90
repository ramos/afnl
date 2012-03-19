
Program Seed

  USE NumTypes
  USE Statistics
  USE Time

  Real (kind=DP), Allocatable :: X(:)
  Character (len=50) :: CC(5)
  Integer :: is

  CC(1) = 'Follow Input'
  CC(2) = 'a=1'
  CC(3) = 'b=2'
  CC(4) = 'c=23'
  CC(5) = 'd=276278'

  CALL PutLuxSeed()
  is = LUXRandomize(1507370652)
  Write(*,*)'Reduced seed: ', is
  Allocate(X(1000))
  CALL Normal(X)
  CALL WriteLuxSeed('o.seed', CC)

  CALL Normal(X)
  Write(*,*)X(990:1000)

  Write(*,*)'Again: '
  CALL ReadLuxSeed('o.seed')
  CALL Normal(X)
  Write(*,*)X(990:1000)
  


  Stop
End Program Seed
