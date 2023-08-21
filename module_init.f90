module soil_init

use module_mp_soil_init
use machine, only : kind_phys

contains

      subroutine noahmpsoilinit (lsm_cold_start, im, lsoil_lsm, lsoil,    & ! in
                                 zsi,zso, dzso,tskin_lnd,tg3, smc, slc, stc,    & ! in
                                 sh2o,tslb, smois, stype, vtype,          & ! in
                                 errmsg, errflg)

      use namelist_soilveg
      use noahmp_tables, only:smcref_table,smcdry_table

      implicit none

      logical,                                        intent(in   ) :: lsm_cold_start
      integer,                                        intent(in   ) :: im, lsoil
      integer,                                        intent(in   ) :: lsoil_lsm
      real (kind=kind_phys), dimension(im),           intent(in   ) :: tskin_lnd
      real (kind=kind_phys), dimension(im),           intent(in   ) :: tg3
      real (kind=kind_phys), dimension(1:lsoil),      intent(in   ) :: zsi
      real (kind=kind_phys), dimension(1:lsoil_lsm),  intent(in   ) :: zso
      real (kind=kind_phys), dimension(1:lsoil_lsm),  intent(in   ) :: dzso
      real (kind=kind_phys), dimension(im,lsoil),     intent(in   ) :: smc !  input
      real (kind=kind_phys), dimension(im,lsoil),     intent(in   ) :: stc !  input
      real (kind=kind_phys), dimension(im,lsoil),     intent(in   ) :: slc !  input

      integer,               dimension(im),    intent(in) :: stype, vtype

      real (kind=kind_phys), dimension(im,lsoil_lsm), intent(inout) :: smois! lsoil_lsm
      real (kind=kind_phys), dimension(im,lsoil_lsm), intent(inout) :: tslb ! lsoil_lsm
      real (kind=kind_phys), dimension(im,lsoil_lsm), intent(inout) :: sh2o ! lsoil_lsm

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!> local

      logical :: debug_print
      logical :: smadj    ! for soil mosture adjustment
      logical :: swi_init ! for initialization in terms of SWI (soil wetness index)

      real (kind=kind_phys),    dimension(1:lsoil_lsm)  :: factorsm
      real (kind=kind_phys),    dimension(im)           :: smcref2
      real (kind=kind_phys),    dimension(im)           :: smcwlt2

      integer,   parameter :: slcats = 30
      real                 :: refsmc1,wltsmc1
      real                 :: REFSMCnoah(slcats),WLTSMCnoah(slcats)


      integer , dimension( 1:im , 1:1 )      :: ivgtyp
      integer , dimension( 1:im , 1:1)       :: isltyp
      real (kind=kind_phys),    dimension( 1:im , 1:1 )       :: tsk
      real (kind=kind_phys),    dimension( 1:im , 1:1 )       :: tbot
      real (kind=kind_phys),    dimension( 1:im , 1:1 )       :: smtotn
      real (kind=kind_phys),    dimension( 1:im , 1:1 )       :: smtotr
      real (kind=kind_phys),    dimension( 1:im , 1:lsoil_lsm, 1:1 ) :: dumsm 
      real (kind=kind_phys),    dimension( 1:im , 1:lsoil_lsm, 1:1 ) :: dumt
!     real (kind=kind_phys),    dimension( 1:im , 1:lsoil_lsm, 1:1 ) :: smfr
      real (kind=kind_phys),    dimension( 1:im , 1:lsoil_lsm, 1:1 ) :: soilm
      real (kind=kind_phys),    dimension( 1:im , 1:lsoil_lsm, 1:1 ) :: soiltemp
      real (kind=kind_phys),    dimension( 1:im , 1:lsoil_lsm, 1:1 ) :: soilh2o

      real (kind=kind_phys) :: st_input(1:im,1:lsoil_lsm*3,1:1)
      real (kind=kind_phys) :: sm_input(1:im,1:lsoil_lsm*3,1:1)

!     integer,              dimension(1:lsoil)  :: st_levels_input ! 
!     integer,              dimension(1:lsoil)  :: sm_levels_input !

      integer :: i,j,k,l,ii,jj,num_soil_layers


      ! Initialize the CCPP error handling variables

          errmsg = ''
          errflg = 0

          num_soil_layers =  lsoil ! 4 - hard-wired for cold start from Noah lsm for now
          debug_print = .false.

!         st_levels_input = (/ 5, 25, 70, 150/)    ! Noah centers of soil layers
!         sm_levels_input = (/ 5, 25, 70, 150/)    ! Noah centers of soil layers

      if (lsm_cold_start) then     !need to change if warm start

          smadj = .true.
          swi_init = .true.

!       if (lsoil /= 4 ) then
!        write (0,*)'lsoil, lsoil =',lsoil
!         errflg = 1
!         return
!       endif


      if(debug_print) then
         write (0,*)'lsm_cold_start = ',lsm_cold_start
         write (0,*)'lsoil, lsoil_lsm =',lsoil, lsoil_lsm
      endif

      endif  !cold start


        do i=1,im ! i = horizontal loop; extra dimension index?

            tbot(i,1) = tg3(i)
            ivgtyp(i,1) = vtype(i)
            isltyp(i,1) = stype(i)
            tsk(i,1) = tskin_lnd(i)
        enddo

      ! Noah lsm input

         do i = 1,slcats

           if (SATDK(i) /= 0.0 .and. BB(i) > 0.0) then

           REFSMC1   = MAXSMC(I)*(5.79E-9/SATDK(I))  &
                        **(1.0/(2.0*BB(I)+3.0))
           REFSMCnoah(I) = REFSMC1 + (MAXSMC(I)-REFSMC1) / 6.
           WLTSMC1   = MAXSMC(I) *                       &
           (200.0/SATPSI(I))**(-1.0/BB(I))
           WLTSMCnoah(I) = WLTSMC1 - 0.5 * WLTSMC1

           endif
         end do


        do i=1,im ! i = horizontal loop

          st_input(i,1,1)=tsk(i,1)
          sm_input(i,1,1)=0.

          !--- initialize smcwlt2 and smcref2 using Noah values
            ii=stype(i)
            if(ii.le.1)ii=1
            smcref2 (i) = refsmcNoah(ii)
            smcwlt2 (i) = wltsmcNoah(ii)

          do k=1,lsoil
             st_input(i,k+1,1)=stc(i,k)
             ! convert volumetric soil moisture to SWI (soil wetness index)
             if( swi_init) then
               sm_input(i,k+1,1)=min(1.,max(0.,(smc(i,k) - smcwlt2(i))/  &
                                 (smcref2(i) - smcwlt2(i))))
             endif
          enddo

          do k=lsoil+2,lsoil_lsm * 3
             st_input(i,k,1)=0.
             sm_input(i,k,1)=0.
          enddo

        enddo ! i - horizontal loop

        CALL init_soil_3_mp ( im, tsk , tbot , dumsm , dumt ,           &
                                st_input , sm_input ,                   &
                                zsi, zso ,                              &
                                lsoil_lsm , num_soil_layers,            &
                                num_soil_layers,                        &
                                lsoil_lsm * 3 , lsoil_lsm * 3)

        do i=1,im
           do k=1,lsoil_lsm
           ! convert from SWI to Noah MP volumetric soil moisture
             if(swi_init) then
               soilm(i,k,1) = dumsm(i,k,1) *                            &
                 (smcref_table(isltyp(i,1))-smcdry_table(isltyp(i,1)))  &
                 + smcdry_table(isltyp(i,1))
             endif
             soiltemp(i,k,1) = dumt(i,k,1)
           enddo ! k
        enddo ! im loop

        ! smadj should be true when the Noah LSM is used

        if( smadj ) then

        ! With other LSMs as input, or when RUC soil moisture is cycled, it
        ! should be set to .false.

          do i=1,im

            ! initialize factor

            do k=1,lsoil_lsm
               factorsm(k)=1.
            enddo

            ! Noah MP soil moisture bucket

            smtotr(i,1)=0.

            do k=1,lsoil_lsm -1
              smtotr(i,1)=smtotr(i,1) + soilm(i,k,1) *dzso(k)
            enddo

            ! Noah soil moisture bucket 

!           smtotn(i,1)=smc(i,1)*0.1 + smc(i,2)*0.2 + smc(i,3)*0.7 + smc(i,4)*1.
            smtotn(i,1)=smc(i,1)*0.1 + smc(i,2)*0.3 + smc(i,3)*0.6 + smc(i,4)*1. ! the depths of 2 and 3 are corrected


            ! Noah MP soil moisture correction to match Noah soil moisture bucket

            do k=1,lsoil_lsm-1
              soilm(i,k,1) = max(0.02,soilm(i,k,1)*smtotn(i,1)/(0.9*smtotr(i,1)))
            enddo

            if( soilm(i,2,1) > soilm(i,1,1) .and. soilm(i,3,1) > soilm(i,2,1)) then
            ! typical for daytime, no recent precip
              factorsm(1) = 0.75
              factorsm(2) = 0.8
              factorsm(3) = 0.85
              factorsm(4) = 0.9
              factorsm(5) = 0.95
            endif

            do k=1,lsoil_lsm
               soilm(i,k,1) = factorsm(k) * soilm(i,k,1)
            enddo

               smtotr(i,1) = 0.

            do k=1,lsoil_lsm - 1
               smtotr(i,1)=smtotr(i,1) + soilm(i,k,1) *dzso(k)
            enddo

          enddo ! i loop

        endif ! smadj==.true.

      ! Initialize liquid soil moisture from total soil moisture
      ! and soil temperature

      call initslc(im,lsoil_lsm, isltyp, ivgtyp,               &
                 soilh2o, soiltemp, soilm)

      do i=1,im
        do k = 1, lsoil_lsm
          smois(i,k) = soilm(i,k,1)
          tslb(i,k)  = soiltemp(i,k,1)
          sh2o(i,k)  = soilh2o(i,k,1)
        enddo 
      enddo

      end subroutine noahmpsoilinit

end module soil_init
