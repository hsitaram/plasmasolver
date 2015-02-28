module chem_module

      use fundconstants
      implicit none

      integer, parameter :: nspecies=6
      integer, parameter :: bgspecnum=1   !background gas
      integer, parameter :: especnum=2	  !electron species
      integer, parameter :: nreac=7
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

	integer :: O2,E,Op,O2p,Om,O
	integer :: rnum
	
	reactants  = 0
	products   = 0
	elecenergy = 0.d0
	gasenergy  = 0.d0
	k_arrh     = 0.d0
	isratearrh = .true.

	O2  = 1
	E   = 2
	Op  = 3
	O2p = 4
	Om  = 5
	O   = 6

	!convention for energy is added is positive and
	!removed is negative.
	!eg: electron loses 19.8 eV when reaction 
	!1 happens

	rnum=1
	!***********************************
	!O_2 + E --> O_2^+ + 2E
	!***********************************
	reactants(rnum,O2)  = 1
	reactants(rnum,E)   = 1
	
	products(rnum,O2p)  = 1	
	products(rnum,E)    = 2
	!-----------------------------------
	isratearrh(rnum) = .false.
	k_arrh(1,rnum)   = 9.0d-10*(cm_to_m**3)/(eVtoK**2)
	k_arrh(2,rnum)   = 2.d0
	k_arrh(3,rnum)   = 146216.7
	!-----------------------------------
	elecenergy(rnum) = -12.6
	gasenergy(rnum)  =  0.d0
	!+++++++++++++++++++++++++++++++++++


	rnum=rnum+1
	!***********************************
	!O_2 + E --> 2 O + E
	!***********************************
	reactants(rnum,O2) = 1
	reactants(rnum,E ) = 1

	products(rnum,O)   = 2
	products(rnum,E)   = 1
	!-----------------------------------
	isratearrh(rnum) = .false.
	k_arrh(1,rnum)   = 4.23d-9*(cm_to_m**3)/(eVtoK**0.0)
	k_arrh(2,rnum)   = 0.d0
	k_arrh(3,rnum)   = 64521.02
	!-----------------------------------
	elecenergy(rnum) = -5.56
        gasenergy(rnum)  = 0.0
	!+++++++++++++++++++++++++++++++++++


	rnum=rnum+1
	!***********************************
	!O_2 + E --> O + O-
	!***********************************
	reactants(rnum,O2)  = 1
	reactants(rnum,E)   = 1
	
	products(rnum,O)    = 1
	products(rnum,Om)   = 1
	!-----------------------------------
	isratearrh(rnum) = .false.
	k_arrh(1,rnum)   = 4.6d-11*(cm_to_m**3)/(eVtoK**0.0)
	k_arrh(2,rnum)   = 0.0
	k_arrh(3,rnum)   = 33769.095
	!-----------------------------------
	elecenergy(rnum) = -3.6
	gasenergy(rnum)  =  0.d0
	!+++++++++++++++++++++++++++++++++++


	rnum=rnum+1
	!***********************************
	!O + E --> O+ + 2E
	!***********************************
	reactants(rnum,O)   = 1
	reactants(rnum,E)   = 1
	
	products(rnum,Op)   = 1
	products(rnum,E)    = 2
	!-----------------------------------
	isratearrh(rnum) = .true.
	k_arrh(1,rnum)   = 9.d-9*(cm_to_m**3)/(eVtoK**0.7)
	k_arrh(2,rnum)   = 0.7
	k_arrh(3,rnum)   = 157821.2
	!-----------------------------------
	elecenergy(rnum) = -13.6
	gasenergy(rnum)  =  0.d0
	!+++++++++++++++++++++++++++++++++++


	rnum=rnum+1
	!***********************************
	!O- + O_2^+ --> O + O_2
	!***********************************
	reactants(rnum,Om)    = 1
	reactants(rnum,O2p)   = 1
	
	products(rnum,O)      = 1
	products(rnum,O2)     = 1
	!-----------------------------------
	isratearrh(rnum) = .true.
	!this reaction rate is calculated at 300 K gas 
	!temperature from Doug's thesis
	k_arrh(1,rnum)   =  1.99d-7*(cm_to_m**3)
	k_arrh(2,rnum)   =  0.0
	k_arrh(3,rnum)   =  0.d0
	!-----------------------------------
	elecenergy(rnum) = 0.d0
	gasenergy(rnum)  = 10.69
	!+++++++++++++++++++++++++++++++++++


	rnum=rnum+1
	!***********************************
	!O- + O+ --> O + O
	!***********************************
	reactants(rnum,Om) = 1
	reactants(rnum,Op) = 1

	products(rnum,O)   = 2
	!-----------------------------------
	isratearrh(rnum) = .true.
	!this reaction rate is calculated at 300 K gas 
	!temperature from Doug's thesis
	k_arrh(1,rnum)   = 2.66d-7*(cm_to_m**3)
	k_arrh(2,rnum)   = 0.0
	k_arrh(3,rnum)   = 0.d0
	!-----------------------------------
	elecenergy(rnum) = 0.d0
	gasenergy(rnum)  = 12.69
	!+++++++++++++++++++++++++++++++++++


	rnum=rnum+1
	!***********************************
	!O- + E --> O + 2E
	!***********************************
	reactants(rnum,Om) = 1
	reactants(rnum,E)  = 1

	products(rnum,O)   = 1
	products(rnum,E)   = 2
	!-----------------------------------
	isratearrh(rnum) = .true.
	k_arrh(1,rnum)   = 2.10d-10*(cm_to_m**3)
	k_arrh(2,rnum)   = 0.5
	k_arrh(3,rnum)   = 39400.0
	!-----------------------------------
	elecenergy(rnum) = -0.92
	gasenergy(rnum)  = 0.d0
	!+++++++++++++++++++++++++++++++++++

	if(rnum .ne. nreac) then
		print *,"rnum not equal to nreac"
		stop
	endif


end subroutine assignreactions
!====================================================================
subroutine setspecparams()

	      specnames(1) = 'O2'
	      specnames(2) = 'E'
	      specnames(3) = 'O+'
	      specnames(4) = 'O2+'
	      specnames(5) = 'O-'
	      specnames(6) = 'O'

	      molmass(1) = 32.d0*mass_prot
	      molmass(2) = mass_elec
	      molmass(3) = 16.d0*mass_prot
	      molmass(4) = 32.d0*mass_prot
	      molmass(5) = 16.d0*mass_prot
	      molmass(6) = 16.d0*mass_prot

	      spec_charge(1) =  0.d0
	      spec_charge(2) = -1.d0
	      spec_charge(3) =  1.d0
	      spec_charge(4) =  1.d0
	      spec_charge(5) =  -1.d0
	      spec_charge(6) =  0.d0

	      ionspecmin = 3
	      ionspecmax = 5

	      no_of_ions = 3

	      neutralspecmin = 6
	      neutralspecmax = 6

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
		dcoeff=0.02
	else if(specnum .eq. 4) then
		dcoeff=0.01
	else if(specnum .eq. 5) then
		dcoeff=0.02
	else if(specnum .eq. 6) then
		dcoeff=0.02
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
	else if(specnum .eq. 5) then
		mobility = -0.4
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

	if(reacnum .eq. 1) then
		!Fit from shankar's thesis	
		rateconst = 1.043d-13
		expterm = -2.74d5/T - 2.0015d8/(T**2) + 5.33d13/(T**3) - 4.001d17/(T**4)
		rateconst=rateconst*exp(expterm)

	else if(reacnum .eq. 2) then
		
		!Fit from shankar's thesis	
		rateconst = 9.577d-16
		expterm = 1.246d5/T - 8.647d9/(T**2) + 1.381d14/(T**3) - 6.953d17/(T**4)
		rateconst=rateconst*exp(expterm)

	else if(reacnum .eq. 3) then
		!Fit from shankar's thesis	
		rateconst = 1.46d-14
		expterm = 6.21d4/T - 7.27d9/(T**2) + 1.25d14/(T**3) - 6.57d17/(T**4)
		rateconst=rateconst*exp(expterm)
	else
		print *,"No custom rate for reaction",reacnum
		stop
	endif

end function getcustomrate
!====================================================================

end module chem_module
