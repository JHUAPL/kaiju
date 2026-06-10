! Dragon King definitions/constants

module dkdefs
  use kdefs
  implicit none

  enum, bind(C)
     enumerator :: FEDDER=1,LINMRG,SUNNY
  end enum

  enum, bind(C)
     enumerator :: LOMPE=1
  end enum

  enum, bind(C)
     enumerator :: NONE=1,ADHOC,CMIT
  end enum

  enum, bind(C)
     enumerator :: KAE23=1,ROB87
  end enum

  enum, bind(C)
     enumerator :: AT_NoPre=0,AT_RMnoE,AT_RMfnE,AT_RMono
  end enum

!  real(rp), parameter :: mixeTINY = 1.D-8 ! Floor of average energy [keV]

end module dkdefs
