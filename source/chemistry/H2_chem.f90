module chem_module

      use fundconstants
      implicit none

      integer, parameter :: nspecies=5
      integer, parameter :: bgspecnum=1   !background gas
      integer, parameter :: especnum=2	  !electron species
      integer, parameter :: nreac=8
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

	integer :: H2,E,Hp,H2p,H
	integer :: rnum
	
	reactants  = 0
	products   = 0
	elecenergy = 0.d0
	gasenergy  = 0.d0
	k_arrh     = 0.d0
	isratearrh = .true.

	H2  = 1
	E   = 2
	Hp  = 3
	H2p = 4
	H   = 5

	!convention for energy is added is positive and
	!removed is negative.
	!eg: electron loses 19.8 eV when reaction 
	!1 happens

	rnum=1
	!***********************************
	!H + E --> H^+ + 2E
	!***********************************
	reactants(rnum,H) = 1
	reactants(rnum,E) = 1
	
	products(rnum,Hp)  = 1	
	products(rnum,E)   = 2
	!-----------------------------------
	isratearrh(rnum) = .true.
	k_arrh(1,rnum)   = 6.5023d-9*(cm_to_m**3)/(eVtoK**0.48931)
	k_arrh(2,rnum)   = 0.48931
	k_arrh(3,rnum)   = 149624.36
	!-----------------------------------
	elecenergy(rnum) = -12.89365
	gasenergy(rnum)  =  0.d0
	!+++++++++++++++++++++++++++++++++++


	rnum=rnum+1
	!***********************************
	!H_2 + E --> H^+ + H + 2E
	!***********************************
	reactants(rnum,H2) = 1
	reactants(rnum,E ) = 1

	products(rnum,Hp)  = 1
	products(rnum,H)   = 1
	products(rnum,E)   = 2
	!-----------------------------------
	isratearrh(rnum) = .true.
	k_arrh(1,rnum)   = 2.9962d-8*(cm_to_m**3)/(eVtoK**0.44456)
	k_arrh(2,rnum)   = 0.44456
	k_arrh(3,rnum)   = 437818.75
	!-----------------------------------
	elecenergy(rnum) = -37.72836
        gasenergy(rnum)  = 0.0
	!+++++++++++++++++++++++++++++++++++


	rnum=rnum+1
	!***********************************
	!H_2^+ + E --> H^+ + H + E
	!***********************************
	reactants(rnum,H2p) = 1
	reactants(rnum,E)   = 1
	
	products(rnum,Hp)   = 1
	products(rnum,H)    = 1
	products(rnum,E)    = 1
	!-----------------------------------
	isratearrh(rnum) = .true.
	k_arrh(1,rnum)   = 1.0702d-7*(cm_to_m**3)/(eVtoK**0.04876)
	k_arrh(2,rnum)   = 0.04876
	k_arrh(3,rnum)   = 112450.85
	!-----------------------------------
	elecenergy(rnum) = -9.69028
	gasenergy(rnum)  =  0.d0
	!+++++++++++++++++++++++++++++++++++


	rnum=rnum+1
	!***********************************
	!H_2^+ + E --> H+ + H+ + 2E
	!***********************************
	reactants(rnum,H2p) = 1
	reactants(rnum,E)   = 1
	
	products(rnum,Hp)  = 2
	products(rnum,E)   = 2
	!-----------------------------------
	isratearrh(rnum) = .true.
	k_arrh(1,rnum)   = 2.1202d-9*(cm_to_m**3)/(eVtoK**0.31394)
	k_arrh(2,rnum)   = 0.31394
	k_arrh(3,rnum)   = 270371.5
	!-----------------------------------
	elecenergy(rnum) = -23.29885
	gasenergy(rnum)  =  0.d0
	!+++++++++++++++++++++++++++++++++++


	rnum=rnum+1
	!***********************************
	!H_2^+ + H --> H_2 + H^+
	!***********************************
	reactants(rnum,H2p) = 1
	reactants(rnum,H)   = 1
	
	products(rnum,H2)   = 1
	products(rnum,Hp)   = 1
	!-----------------------------------
	isratearrh(rnum) = .true.
	k_arrh(1,rnum)   = 9.0d-10*(cm_to_m**3)
	k_arrh(2,rnum)   = 0.0
	k_arrh(3,rnum)   = 0.0
	!-----------------------------------
	elecenergy(rnum) = 0.d0
	gasenergy(rnum)  = 0.d0
	!+++++++++++++++++++++++++++++++++++


	rnum=rnum+1
	!***********************************
	!H_2 + H+ --> H_2+ + H
	!***********************************
	reactants(rnum,H2) = 1
	reactants(rnum,Hp) = 1

	products(rnum,H2p)  = 1
	products(rnum,H)    = 1
	!-----------------------------------
	isratearrh(rnum) = .true.
	k_arrh(1,rnum)   = 1.19d-22*(cm_to_m**3)
	k_arrh(2,rnum)   = 0.d0
	k_arrh(3,rnum)   = 0.d0
	!-----------------------------------
	elecenergy(rnum) = 0.d0
	gasenergy(rnum)  = 0.d0
	!+++++++++++++++++++++++++++++++++++


	rnum=rnum+1
	!***********************************
	!H_2 + E --> H_2^+ + 2E
	!***********************************
	reactants(rnum,H2) = 1
	reactants(rnum,E)  = 1

	products(rnum,H2p)  = 1
	products(rnum,E)    = 2
	!-----------------------------------
	isratearrh(rnum) = .true.
	k_arrh(1,rnum)   = 3.1228d-8*(cm_to_m**3)/(eVtoK**0.17156)
	k_arrh(2,rnum)   = 0.17156
	k_arrh(3,rnum)   = 232987.49
	!-----------------------------------
	elecenergy(rnum) = -20.07734
	gasenergy(rnum)  = 0.d0
	!+++++++++++++++++++++++++++++++++++


	rnum=rnum+1
	!***********************************
	!H_2 + E --> 2H + E
	!***********************************
	reactants(rnum,H2) = 1
	reactants(rnum,E)  = 1

	products(rnum,H)   = 2
	products(rnum,E)   = 1
	!-----------------------------------
	isratearrh(rnum) = .true.
	k_arrh(1,rnum)   = 1.7527d-7*(cm_to_m**3)/(eVtoK**(-1.23668))
	k_arrh(2,rnum)   = -1.23668
	k_arrh(3,rnum)   = 146128.85
	!-----------------------------------
	elecenergy(rnum) = -12.59243
	gasenergy(rnum)  =   0.0
	!+++++++++++++++++++++++++++++++++++

	if(rnum .ne. nreac) then
		print *,"rnum not equal to nreac"
		stop
	endif


end subroutine assignreactions
!====================================================================
subroutine setspecparams()

	      specnames(1) = 'H2'
	      specnames(2) = 'E'
	      specnames(3) = 'H+'
	      specnames(4) = 'H2+'
	      specnames(5) = 'H'

	      molmass(1) = 2.d0*mass_prot
	      molmass(2) = mass_elec
	      molmass(3) = 1.d0*mass_prot
	      molmass(4) = 2.d0*mass_prot
	      molmass(5) = 1.d0*mass_prot

	      spec_charge(1) =  0.d0
	      spec_charge(2) = -1.d0
	      spec_charge(3) =  1.d0
	      spec_charge(4) =  1.d0
	      spec_charge(5) =  0.d0

	      ionspecmin = 3
	      ionspecmax = 4

	      no_of_ions = 2

	      neutralspecmin = 5
	      neutralspecmax = 5

	      no_of_neutrals = 1

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
		dcoeff=0.02
	else if(specnum .eq. 4) then
		dcoeff=0.01
	else if(specnum .eq. 5) then
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
		mobility = 0.4
	else if(specnum .eq. 4) then
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

	real*8 :: alpha,E0

	alpha     = kparams(1)
	E0        = kparams(3)

	rateconst = alpha*abs(efield)*exp(-E0/abs(efield))

end function getcustomrate
!====================================================================

end module chem_module
