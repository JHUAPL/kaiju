! MIX definitions/constants

module mixdefs
  use kdefs
  implicit none

  integer, parameter :: nVars = 27 ! change together with the enumerator below
  enum, bind(C)
     enumerator :: POT=1,FAC,SIGMAP,SIGMAH,SOUND_SPEED,DENSITY,AVG_ENG,NUM_FLUX,NEUTRAL_WIND,EFIELD,IM_EAVG,IM_EFLUX,IM_IAVG,IM_IFLUX,Z_EAVG,Z_NFLUX,CRPOT,TPOT,IM_GTYPE,AUR_TYPE,IM_BETA,IM_EDEN,IM_EPRE,IM_ENFLX,IM_INFLX,DELTAE,IM_NPSP
  end enum

  ! enumerator for MHD->MIX variables
  enum, bind(C)
     enumerator :: MHDJ=1, MHDD, MHDC
  end enum

  ! enumerator for MIX->MHD variable(s)
  enum, bind(C)
     enumerator :: MHDPSI=1
  end enum

  enum, bind(C)
     enumerator :: LOMPE=1
  end enum

  enum, bind(C)
     enumerator :: FEDDER=1,LINMRG,SUNNY
  end enum

  enum, bind(C)
     enumerator :: NONE=1,ADHOC,CMIT
  end enum

  enum, bind(C)
     enumerator :: AT_NoPre=0,AT_RMnoE,AT_RMfnE,AT_RMono
  end enum

  ! enumerator for transform variables
  enum, bind(C)
    enumerator :: iSMtoGEO,iGEOtoSM,iSMtoGSM,iGSMtoSM,iGEOtoGSM,iGSMtoGEO,iAPEXtoSM,iSMtoAPEX
  end enum

  real(rp), parameter :: mixeTINY = 1.D-8 ! Floor of average energy [keV]

end module mixdefs
