module chem_module

      use fundconstants
      implicit none

      integer, parameter :: nspecies=3
      integer, parameter :: bgspecnum=1   !background gas
      integer, parameter :: especnum=2	  !electron species
      integer, parameter :: nreac=1
      real*8, parameter  :: cm_to_m = 0.01
      character (LEN=10), dimension(nspecies) :: specnames
      integer :: reactants(nreac,nspecies)
      integer :: products(nreac,nspecies)
      real*8  :: k_arrh(3,nreac)
      real*8  :: elecenergy(nreac)
      real*8  :: gasenergy(nreac)
      logical :: isratearrh(nreac)
      real*8  :: molmass(nspecies)
      real*8  :: spec_charge(nspecies)
      
      integer :: ionspecmin
      integer :: ionspecmax
      integer :: neutralspecmin
      integer :: neutralspecmax
      integer :: no_of_ions
      integer :: no_of_neutrals

      contains
!====================================================================
subroutine assignreactions()

	integer :: N,E,I
	integer :: rnum
	
	reactants  = 0
	products   = 0
	elecenergy = 0.d0
	gasenergy  = 0.d0
	k_arrh     = 0.d0
	isratearrh = .true.

	N   = 1
	E   = 2
	I   = 3

	!convention for energy is added is positive and
	!removed is negative.
	!eg: electron loses 19.8 eV when reaction 
	!1 happens

	rnum=1
	!***********************************
	!N + E --> I + 2E
	!***********************************
	reactants(rnum,N)   = 1
	reactants(rnum,E)   = 1
	
	products(rnum,I)    = 1	
	products(rnum,E)    = 2
	!-----------------------------------
	isratearrh(rnum) = .true.
	k_arrh(1,rnum)   = 2.5d-6*(cm_to_m**3)
	k_arrh(2,rnum)   = 0.d0
	k_arrh(3,rnum)   = 278508.0
	!-----------------------------------
	elecenergy(rnum) = -15.578
	gasenergy(rnum)  =  0.d0
	!+++++++++++++++++++++++++++++++++++

	if(rnum .ne. nreac) then
		print *,"rnum not equal to nreac"
		stop
	endif


end subroutine assignreactions
!====================================================================
subroutine setspecparams()

	      specnames(1) = 'N'
	      specnames(2) = 'E'
	      specnames(3) = 'I+'

	      molmass(1) = 4.d0*mass_prot
	      molmass(2) = mass_elec
	      molmass(3) = 4.d0*mass_prot

	      spec_charge(1) =  0.d0
	      spec_charge(2) = -1.d0
	      spec_charge(3) =  1.d0

	      ionspecmin = 3
	      ionspecmax = 3

	      no_of_ions = 1

	      neutralspecmin = 0
	      neutralspecmax = -1

	      no_of_neutrals = 0

end subroutine setspecparams
!====================================================================
subroutine initializechemistry()

	call setspecparams()
	call assignreactions()

end subroutine initializechemistry
!====================================================================
subroutine getspecproduction(specnum,Te,Tg,specden,specprod,efld)

	integer, intent(in)  :: specnum
	real*8, intent(in)   :: Te,Tg
	real*8, intent(in)   :: specden(nspecies)
	real*8,intent(inout) :: specprod
	real*8,optional :: efld
	real*8 :: efield
	real*8 :: k
	integer :: i,j
	integer :: nmoles
	integer :: nreacmoles
	real*8 :: specmult

	if(present(efld)) then
		efield=efld
	else
		efield=1.0
	endif
	
	specprod = 0.d0
	do i=1,nreac

		nmoles = products(i,specnum)-reactants(i,specnum)

		if(nmoles .ne. 0) then

			if(isratearrh(i) .eqv. .true.) then
				k=getarrhrate(k_arrh(:,i),Te)
			else
				k=getcustomrate(i,k_arrh(:,i),Te,efield)
			endif

			specmult = 1.d0

			do j=1,nspecies
				nreacmoles = reactants(i,j)
				specmult = specmult*(specden(j)**nreacmoles)
			enddo

			specprod = specprod + nmoles*k*specmult
			!print *,"reac:",i,"specprod:", nmoles*k*specmult
		endif
	enddo

end subroutine getspecproduction
!====================================================================
subroutine getelectroninelasticterm(Te,Tg,specden,inelterm,efld)
	
	real*8, intent(in)   :: Te,Tg
	real*8, intent(in)   :: specden(nspecies)
	real*8, intent(inout) :: inelterm
	real*8,optional :: efld
	integer :: i,j,nreacmoles
	real*8 :: k,specmult
	real*8 :: efield

	inelterm = 0.d0
	
	if(present(efld)) then
		efield=efld
	else
		efield=1.0
	endif
	
	do i=1,nreac

		if(elecenergy(i) .ne. 0) then
			
				if(isratearrh(i) .eqv. .true.) then
					k=getarrhrate(k_arrh(:,i),Te)
				else
					k=getcustomrate(i,k_arrh(:,i),Te,efield)
				endif

				specmult = 1.d0
				do j=1,nspecies
					nreacmoles = reactants(i,j)
					!print *,"nreacmoles at ",j,nreacmoles
					specmult = specmult*(specden(j)**nreacmoles)
				enddo

				inelterm = inelterm + k*specmult*elecenergy(i)
				!print *,"reac:",i,"inelterm:", k,specmult,elecenergy(i),inelterm

		endif
	enddo

end subroutine getelectroninelasticterm
!====================================================================
function getspecdcoeff(specnum,Te,Tg,Pg)  result(dcoeff)
	
	real*8, intent(in) :: Te,Tg,Pg
	integer,intent(in) :: specnum
	real*8 :: dcoeff
	real*8 :: Patm

	Patm = 101325.d0

	if(specnum .eq. 2) then
		dcoeff=100.d0
	else if(specnum .eq. 3) then
		dcoeff=0.01
	else
		write(*,*)"species does not exist"
	endif
end function getspecdcoeff
!====================================================================
function getspecmobility(specnum,Te,Tg,Pg)  result(mobility)
	
	real*8, intent(in) :: Te,Tg,Pg
	integer,intent(in) :: specnum
	real*8 :: mobility
	real*8 :: Patm

	Patm = 101325.d0

	if(specnum .eq. 2) then
		mobility = -20.d0
	else if(specnum .eq. 3) then
		mobility = 0.2
	else
		write(*,*)"species does not exist"
	endif

end function getspecmobility
!====================================================================
function getelectroncollisionfrequency(Te,Tg,Pg) result(collfreq)
	
	real*8, intent(in) :: Te,Tg,Pg
	real*8 :: collfreq
	real*8 :: mu
	real*8 :: m_e

	mu = getspecmobility(especnum,Te,Tg,Pg)
	m_e = molmass(especnum)

	collfreq = -echarge/m_e/mu

end function getelectroncollisionfrequency
!====================================================================
function getarrhrate(k,T) result(rateconst)

	real*8, intent(in) :: k(3)
	real*8, intent(in) :: T
	real*8 :: rateconst

	real*8 :: A,alpha,Ea

	A     = k(1)
	alpha = k(2)
	Ea    = k(3)

	rateconst = A*(T**alpha)*exp(-Ea/T)

end function getarrhrate
!====================================================================
function getcustomrate(reacnum,kparams,T,efield) result(rateconst)
	
	integer,intent(in) :: reacnum
	real*8, intent(in) :: kparams(3)
	real*8, intent(in) :: T,efield
	real*8 :: rateconst
	real*8 :: expterm

	real*8 :: T_in_eV

	T_in_eV = T/(eVtoK)
	expterm = 0.d0
	rateconst = 0.d0

	print *,"No custom rate for reaction",reacnum
	stop

end function getcustomrate
!====================================================================

end module chem_module
