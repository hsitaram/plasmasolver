module plasma_solver

	use fundconstants
	use chem_module 
	use solvergmres
	use convdiff
	implicit none

      real*8 :: P_gas
      real*8 :: T_gas
      real*8 :: N_gas
      real*8 :: collfreq_e
      real*8 :: Tefixed
      real*8 :: Temin
      real*8 :: voltage_L
      real*8 :: voltage_R
      real*8 :: sec_elec_coeff
      
      real*8 :: nscale
      real*8 :: Eescale
      real*8 :: Tescale
      
      real*8 :: k_i
      real*8 :: E_i
      real*8 :: deltaE_i
      real*8,allocatable :: Te(:)
     
      real*8  :: length
      integer :: np
      real*8  :: dx
      real*8  :: dt
      real*8  :: finaltime
      real*8  :: maxtimesteps
      integer :: printit
      integer :: printfileit

      real*8, allocatable :: numden(:,:)
      real*8, allocatable :: phi(:)
      real*8, allocatable :: efield(:)
      real*8, allocatable :: Ee(:)

      real*8 :: n0
      real*8 :: nemin

      integer :: kspdim
      integer :: itmax_restart

      contains
!=============================================================================
subroutine init()

	integer :: i
	integer :: restart_option
	logical :: restart
	character(LEN=100) :: temp
	character(LEN=100) :: restartfname

	integer :: restartfptr,inpfptr,logfileptr
	real*8  :: specinit(nspecies)

	inpfptr  = 23
	restartfptr = 24
	logfileptr = 25

	call initializechemistry()

	open(unit=inpfptr,file="plasma_inputs")
	open(unit=logfileptr,file="plsolver.log")
	
	!discharge gap,time step size,scales
	read(inpfptr,*),temp
	read(inpfptr,*),temp,length
	read(inpfptr,*),temp,np
	read(inpfptr,*),temp,dt
	read(inpfptr,*),temp,maxtimesteps
	read(inpfptr,*),temp,printfileit
	read(inpfptr,*),temp,printit
	read(inpfptr,*),temp,finaltime
	read(inpfptr,*),temp,nscale

	dx = length/(np-1)

	!Restart parameters
	read(inpfptr,*),temp
	read(inpfptr,*),temp,restart_option
	
	if(restart_option .eq. 1) then
		restart = .true.
	else
		restart = .false.
	endif
	read(inpfptr,*),temp,restartfname

	!neutral gas parameters
	read(inpfptr,*),temp	
	read(inpfptr,*),temp,P_gas
	read(inpfptr,*),temp,T_gas
	N_gas = P_gas/k_B/T_gas       ! neutral density (#/m3)
	specinit(1) = N_gas

	!potential equation bcs
	read(inpfptr,*),temp	
	read(inpfptr,*),temp,voltage_L
	read(inpfptr,*),temp,voltage_R

	!Electron temperature,sec elec emission
	read(inpfptr,*),temp
	read(inpfptr,*),temp,Tefixed
	read(inpfptr,*),temp,Temin
	read(inpfptr,*),temp,sec_elec_coeff
	
	!number densities
	read(inpfptr,*),temp
	do i=1,nspecies
		read(inpfptr,*),temp,specinit(i)
	enddo

	nemin = 1.0d11                ! #/m3
	Tescale = echarge/k_B         ! scale Te to be eV
	Eescale = Tescale*nscale*k_B        
	
	!GMRES parameters
	kspdim = 5
	itmax_restart = 5
	
	allocate(numden(np,nspecies))
	allocate(phi(np))
	allocate(efield(np))
	allocate(Te(np))
	allocate(Ee(np))

	if(restart .eqv. .false.) then
		
		phi = 0.d0
		Te  = Tefixed/Tescale

		do i=1,nspecies
			numden(:,i)=specinit(i)/nscale
		enddo
	
		Ee(:)  = 1.5*numden(:,especnum)*Te(:)
	else
		print *,"Restarting not implemented yet"
		stop
		!open(unit=restartfptr,file=trim(restartfname))
		!do i=1,np
		!enddo
		!close(restartfptr)
	endif
	
	write(logfileptr,'(A,E20.10)') "Length:"         ,length
	write(logfileptr,'(A,I20)')    "npoints:"        ,np
	write(logfileptr,'(A,E20.10)') "dx:"             ,dx
	write(logfileptr,'(A,E20.10)') "dt:"             ,dt
	write(logfileptr,'(A,E20.10)') "maxtsteps:"      ,maxtimesteps
	write(logfileptr,'(A,E20.10)') "finaltime:"      ,finaltime
	write(logfileptr,'(A,E20.10)') "nscale:"         ,nscale
	write(logfileptr,*)
	write(logfileptr,'(A,I20)') "restart option:"    ,restart_option
	write(logfileptr,'(A,A)') "restart file:"        ,restartfname
	write(logfileptr,*)
	write(logfileptr,'(A,E20.10)') "Pressure:"       ,P_gas
	write(logfileptr,'(A,E20.10)') "Temperature:"    ,T_gas
	write(logfileptr,'(A,E20.10)') "Voltage_L:"      ,voltage_L
	write(logfileptr,'(A,E20.10)') "Voltage_R:"      ,voltage_R
	write(logfileptr,'(A,E20.10)') "Tefixed:"        ,Tefixed
	write(logfileptr,'(A,E20.10)') "Temin:"          ,Temin
	write(logfileptr,'(A,E20.10)') "Sec elec coeff:" ,sec_elec_coeff
	write(logfileptr,*)
	write(logfileptr,*)"Number densities"
	do i=1,nspecies
		write(logfileptr,'(A,E20.10)')specnames(i),&
			       specinit(i)
	enddo	       
	
	close(logfileptr)
	close(inpfptr)
	
	call printfile(0)


end subroutine init
!=============================================================================
subroutine potentialsolve(printflag,residual)
	
	real*8, intent(out) :: residual
	logical, intent(in) :: printflag
	integer :: i,j
	logical :: dirc_bc_flags(2),flux_bc_flags(2)
	real*8  :: dircvals(2),fluxvals(2)
	real*8  :: timederivfactor
	real*8 :: dcoeff(np),vel(np),reac(np),source(np)
	real*8 :: b(np),phiold(np)
	real*8 :: const_1,Zi_ni_sum,sterm
	
	dcoeff=0.d0
	vel=0.d0
	reac=0.d0
	source=0.d0

	if(printflag .eqv. .true.) then
		print *,"______solving potential equation______"
	endif
	
	timederivfactor=0.d0
	const_1 = echarge*nscale/epsilon_0
      	
	do i=1,np

	 	vel(i)    =  0.d0
        	dcoeff(i) =  1.d0
      	 	reac(i)   =  0.d0
		
		sterm = numden(i,especnum)
		Zi_ni_sum =  0.d0
		do j=ionspecmin,ionspecmax
			Zi_ni_sum = Zi_ni_sum &
			+ spec_charge(j)*numden(i,j)
		enddo
		sterm = sterm - Zi_ni_sum
		source(i) = -const_1*sterm
      	enddo

	dirc_bc_flags(1) = .true.
	dirc_bc_flags(2) = .true.
	flux_bc_flags(1) = .false.
	flux_bc_flags(2) = .false.

	dircvals(1) = voltage_L
	dircvals(2) = voltage_R
	fluxvals(1) = 0.d0
	fluxvals(2) = 0.d0

	phiold=phi

	call findrhs(b,phiold,timederivfactor,source,dirc_bc_flags, &
			 flux_bc_flags,dircvals,fluxvals,dx,dt,np)
 	
	call performgmres(b,phiold,phi,timederivfactor,&
			vel,dcoeff,reac,dirc_bc_flags,&
			flux_bc_flags,dircvals,fluxvals,dx,dt,&
			kspdim,np,itmax_restart,findAX,mgridprecond,printflag,residual)
	
	call compute_efield()

end subroutine potentialsolve
!=============================================================================
subroutine compute_efield()

	integer :: i
		
	do i=2,np-1
		efield(i)=-(phi(i+1)-phi(i-1))/(2.d0*dx)
	enddo
	efield(1)  = -(phi(2) -phi(1)   )/dx
	efield(np) = -(phi(np)-phi(np-1))/dx 

end subroutine compute_efield
!=============================================================================
subroutine ionsolve(ispecnum,printflag,residual)
	
	integer,intent(in) :: ispecnum
	logical,intent(in) :: printflag
	real*8, intent(out) :: residual
 	
	integer :: i,j

	logical :: dirc_bc_flags(2),flux_bc_flags(2)
	real*8  :: dircvals(2),fluxvals(2)
	real*8  :: timederivfactor
	real*8  ::  dcoeff(np),vel(np),reac(np),source(np)
	real*8  ::  b(np),niold(np),ninew(np)

	real*8  ::  mu_i
	real*8  :: specarray(nspecies)
	real*8  :: specprod
	
	dcoeff=0.d0
	vel=0.d0
	reac=0.d0
	source=0.d0

	if(printflag .eqv. .true.) then
		print *,"______solving ",trim(specnames(ispecnum))," ion density equation______"
	endif

	timederivfactor=1.d0
      	do i=1,np
		mu_i      = getspecmobility(ispecnum,Te(i)*Tescale,T_gas,P_gas)
	 	vel(i)    = mu_i*efield(i)
        	dcoeff(i) = getspecdcoeff(ispecnum,Te(i)*Tescale,T_gas,P_gas)
      	 	
		do j=1,nspecies
			specarray(j)=numden(i,j)*nscale
		enddo

		call getspecproduction(ispecnum,Te(i)*Tescale,T_gas,&
				specarray,specprod,efield(i))
		
		source(i) = specprod/nscale

      	enddo

	dirc_bc_flags(1) = .false.
	dirc_bc_flags(2) = .false.
	flux_bc_flags(1) = .true.
	flux_bc_flags(2) = .true.

	dircvals(1) = 0.d0
	dircvals(2) = 0.d0
	mu_i        = getspecmobility(ispecnum,Te(1)*Tescale,T_gas,P_gas)
	if(mu_i*efield(1) .lt. 0) then
		fluxvals(1) = mu_i*numden(1,ispecnum)*efield(1)
        else
		fluxvals(1) = 0.d0
	endif

	mu_i        = getspecmobility(ispecnum,Te(np)*Tescale,T_gas,P_gas)
	if(mu_i*efield(np) .gt. 0) then
		fluxvals(2) = mu_i*numden(np,ispecnum)*efield(np)
	else
		fluxvals(2) = 0.d0
	endif

	niold=numden(:,ispecnum)

	call findrhs(b,niold,timederivfactor,source,dirc_bc_flags, &
			 flux_bc_flags,dircvals,fluxvals,dx,dt,np)
 	
	call performgmres(b,niold,ninew,timederivfactor,&
			vel,dcoeff,reac,dirc_bc_flags,&
			flux_bc_flags,dircvals,fluxvals,dx,dt,&
			kspdim,np,itmax_restart,findAX,mgridprecond,printflag,residual)

	numden(:,ispecnum)=ninew

end subroutine ionsolve
!=============================================================================
subroutine neutralsolve(nspecnum,printflag,residual)
	
	integer,intent(in) :: nspecnum
	logical,intent(in) :: printflag
	real*8,intent(out)  :: residual
 	
	integer :: i,j

	logical :: dirc_bc_flags(2),flux_bc_flags(2)
	real*8  :: dircvals(2),fluxvals(2)
	real*8  :: timederivfactor
	real*8  ::  dcoeff(np),vel(np),reac(np),source(np)
	real*8  ::  b(np),nold(np),nnew(np)

	real*8  :: specarray(nspecies)
	real*8  :: specprod,cbar1,cbar2
	
	dcoeff=0.d0
	vel=0.d0
	reac=0.d0
	source=0.d0

	if(printflag .eqv. .true.) then
		print *,"______solving neutral ",trim(specnames(nspecnum))," density equation______"
	endif

	timederivfactor=1.d0
      	do i=1,np
        	
		dcoeff(i) = getspecdcoeff(nspecnum,Te(i)*Tescale,T_gas,P_gas)
		
		do j=1,nspecies
			specarray(j)=numden(i,j)*nscale
		enddo

		call getspecproduction(nspecnum,Te(i)*Tescale,T_gas,&
				specarray,specprod,efield(i))
		
		source(i) = specprod/nscale
      	enddo

	dirc_bc_flags(1) = .false.
	dirc_bc_flags(2) = .false.
	flux_bc_flags(1) = .true.
	flux_bc_flags(2) = .true.

	dircvals(1) = 0.d0
	dircvals(2) = 0.d0
	
	cbar1 = sqrt(8.d0*k_B*T_gas /pi/molmass(nspecnum))  ! in m/s
	cbar2 = sqrt(8.d0*k_B*T_gas /pi/molmass(nspecnum))  ! in m/s

	fluxvals(1) = -0.25 * numden(1,nspecnum)  * cbar1
	fluxvals(2) =  0.25 * numden(np,nspecnum) * cbar2
	
	nold=numden(:,nspecnum)

	call findrhs(b,nold,timederivfactor,source,dirc_bc_flags, &
			 flux_bc_flags,dircvals,fluxvals,dx,dt,np)
 	
	call performgmres(b,nold,nnew,timederivfactor,&
			vel,dcoeff,reac,dirc_bc_flags,&
			flux_bc_flags,dircvals,fluxvals,dx,dt,&
			kspdim,np,itmax_restart,findAX,mgridprecond,printflag,residual)

	numden(:,nspecnum)=nnew

end subroutine neutralsolve
!=============================================================================
subroutine electronsolve(printflag,residual)
	
	logical, intent(in) :: printflag
	real*8, intent(out) :: residual

	integer :: i,j
	logical :: dirc_bc_flags(2),flux_bc_flags(2)
	real*8  :: dircvals(2),fluxvals(2)
	real*8  :: timederivfactor
	real*8 :: dcoeff(np),vel(np),reac(np),source(np)
	real*8 :: b(np),neold(np)
	
	real*8 :: mu_e,mu_i
	real*8 :: cbar1,cbar2
	real*8 :: nenew(np)
	real*8 :: specprod
	real*8 :: specarray(nspecies)
	real*8 :: iflux
	real*8 :: m_e

	dcoeff=0.d0
	vel=0.d0
	reac=0.d0
	source=0.d0

	m_e = molmass(especnum)
	
	if(printflag .eqv. .true.) then
		print *,"______solving electron density equation______"
	endif

	timederivfactor=1.d0
      	do i=1,np
		mu_e      =   getspecmobility(especnum,Te(i)*Tescale,T_gas,P_gas)
	 	vel(i)    =   mu_e*efield(i)
        	dcoeff(i) =   getspecdcoeff(especnum,Te(i)*Tescale,T_gas,P_gas)
      	
		source(i) =   0.d0
		do j=1,nspecies
			specarray(j)=numden(i,j)*nscale
		enddo

		call getspecproduction(especnum,Te(i)*Tescale,T_gas,&
				specarray,specprod,efield(i))

		source(i)   =   source(i) + specprod/nscale  
      	enddo

	dirc_bc_flags(1) = .true.
	dirc_bc_flags(2) = .false.
	flux_bc_flags(1) = .false.
	flux_bc_flags(2) = .true.

	dircvals(1) = nemin/nscale
	dircvals(2) = nemin/nscale

	cbar1 = sqrt(8.d0*k_B*Te(1) *Tescale /pi/m_e)  ! in m/s
	cbar2 = sqrt(8.d0*k_B*Te(np)*Tescale /pi/m_e)  ! in m/s

	fluxvals(1) = -0.25 * numden(1,especnum)  * cbar1
	fluxvals(2) =  0.25 * numden(np,especnum) * cbar2

	!sec electron emission, assuming all ions have the
	!same sec electron emission coefficient
	iflux=0.d0
	do j=ionspecmin,ionspecmax

		mu_i = getspecmobility(j,Te(1)*Tescale,T_gas,P_gas)

		if(mu_i*efield(1) .lt. 0) then
			iflux = iflux + mu_i*numden(1,j)*efield(1)
		endif
	enddo
	fluxvals(1) = fluxvals(1) - sec_elec_coeff*iflux

	iflux=0.d0
	do j=ionspecmin,ionspecmax

		mu_i = getspecmobility(j,Te(np)*Tescale,T_gas,P_gas)

		if(mu_i*efield(np) .gt. 0) then
			iflux = iflux + mu_i*numden(np,j)*efield(np)
		endif
	enddo
	fluxvals(2) = fluxvals(2) - sec_elec_coeff*iflux

	neold=numden(:,especnum)

	call findrhs(b,neold,timederivfactor,source,dirc_bc_flags, &
			 flux_bc_flags,dircvals,fluxvals,dx,dt,np)
 	
	call performgmres(b,neold,nenew,timederivfactor,&
			vel,dcoeff,reac,dirc_bc_flags,&
			flux_bc_flags,dircvals,fluxvals,dx,dt,&
			kspdim,np,itmax_restart,findAX,mgridprecond,printflag,residual)

	!floor electron density
	do i=1,np
		if(nenew(i)*nscale .lt. nemin) then
			nenew(i)=nemin/nscale
		endif
	enddo

	numden(:,especnum) = nenew

end subroutine electronsolve
!=============================================================================
subroutine elecenergysolve(printflag,residual)
	
	logical, intent(in) :: printflag
	real*8, intent(out) :: residual

	integer :: i
	logical :: dirc_bc_flags(2),flux_bc_flags(2)
	real*8  :: dircvals(2),fluxvals(2)
	real*8  :: timederivfactor
	real*8 :: dcoeff(np),vel(np),reac(np),source(np)
	real*8 :: mu_e,jheating(np),elastic_col(np)
	real*8 :: b(np),Eeold(np)
	real*8 :: inelastic_col(np)
	real*8 :: cbar1,cbar2
	real*8 :: fivebythree
	real*8 :: m_e
	
	dcoeff=0.d0
	vel=0.d0
	reac=0.d0
	source=0.d0

	fivebythree = 5.d0/3.d0
	m_e = molmass(especnum)


	if(printflag .eqv. .true.) then
		print *,"______solving electron energy equation______"
	endif
	
	call jouleheatingterm(jheating)
	call elastic_colterm(elastic_col)
	call inelastic_colterm(inelastic_col)

	timederivfactor=1.d0
	do i=1,np
		mu_e      =  getspecmobility(especnum,Te(i)*Tescale,T_gas,P_gas)
	 	vel(i)    =  fivebythree*mu_e*efield(i)
        	dcoeff(i) =  fivebythree* &
				getspecdcoeff(especnum,Te(i)*Tescale,T_gas,P_gas)
      	
      	 	source(i) =  jheating(i)
      	 	source(i) =  source(i) -   elastic_col(i)
      	 	source(i) =  source(i) + inelastic_col(i)
      	enddo

	dirc_bc_flags(1) = .true.
	dirc_bc_flags(2) = .false.
	flux_bc_flags(1) = .false.
	flux_bc_flags(2) = .true.

	!fixed temperature of 0.5 eV
	dircvals(1) = 1.5*numden(1,especnum)*(Tefixed/Tescale)
	dircvals(2) = 1.5*numden(np,especnum)*(Tefixed/Tescale)

	cbar1 = sqrt(8.d0*k_B*Te(1) *Tescale /pi/m_e)
	cbar2 = sqrt(8.d0*k_B*Te(np)*Tescale /pi/m_e)

	fluxvals(1) =  -0.25 *numden(1,especnum)  * cbar1* (2.d0 * Te(1) )
	fluxvals(2) =   0.25 *numden(np,especnum) * cbar2* (2.d0 * Te(np))

	Eeold=Ee

	call findrhs(b,Eeold,timederivfactor,source,dirc_bc_flags, &
			 flux_bc_flags,dircvals,fluxvals,dx,dt,np)
 	
	call performgmres(b,Eeold,Ee,timederivfactor,&
			vel,dcoeff,reac,dirc_bc_flags,&
			flux_bc_flags,dircvals,fluxvals,dx,dt,&
			kspdim,np,itmax_restart,findAX,mgridprecond,printflag,residual)

	!update electron temperature
	do i=1,np
	   if(numden(i,especnum)*nscale .lt. nemin) then
		Te(i) = Temin/Tescale
	   else
		Te(i)=Ee(i)/1.5/numden(i,especnum)
		
		if(Te(i)*Tescale .lt. Temin) then
			Te(i)=Temin/Tescale
		endif
	   endif
	enddo

	do i=1,np
		Ee(i)=1.5*numden(i,especnum)*Te(i)
	enddo
		

end subroutine elecenergysolve
!=============================================================================
subroutine jouleheatingterm(jheating)

	real*8, intent(out) :: jheating(np)
	integer :: i
	real*8 :: mu_e,D_e,grad_ne
	real*8 :: scaling
	real*8 :: Gamma_left,Gamma_right
	real*8 :: elec_field_left,elec_field_right,e_den
	real*8 :: ne(np)

	scaling = nscale/Eescale
	ne = numden(:,especnum)
	
	do i=1,np

		Gamma_left  = 0.d0
		Gamma_right = 0.d0
		elec_field_left = 0.d0
		elec_field_right = 0.d0
		
		mu_e    = getspecmobility(especnum,Te(i)*Tescale,T_gas,P_gas)
		D_e     = getspecdcoeff  (especnum,Te(i)*Tescale,T_gas,P_gas)
		
		if(i .gt. 1) then

			grad_ne =  (ne(i)-ne(i-1))/dx
			elec_field_left  = (efield(i)+efield(i-1))*0.5
	
			if(mu_e*elec_field_left .gt. 0) then
				e_den = ne(i-1)
			else
				e_den = ne(i)
			endif
			Gamma_left = (mu_e*e_den*elec_field_left - D_e*grad_ne)
		endif
		
		if(i .lt. np) then
			
			grad_ne =  (ne(i+1)-ne(i))/dx
			elec_field_right  = (efield(i)+efield(i+1))*0.5
			
			if(mu_e*elec_field_right .gt. 0) then
				e_den = ne(i)
			else
				e_den = ne(i+1)
			endif
			Gamma_right = (mu_e*e_den*elec_field_right - D_e*grad_ne)

		endif
		
	
		if(i .eq. 1) then	
			jheating(i) = -echarge*elec_field_right*Gamma_right*scaling
		else if(i .eq. np) then
			jheating(i) = -echarge*elec_field_left*Gamma_left*scaling
		else
			jheating(i) = -echarge*0.5*(elec_field_left+elec_field_right)*&
					0.5*(Gamma_left+Gamma_right)*scaling
		endif

		if(jheating(i) .lt. 0.d0) then
			jheating(i)=0.d0
		endif
	enddo
		

end subroutine jouleheatingterm
!=============================================================================
subroutine elastic_colterm(elastic_col)

	real*8, intent(out) :: elastic_col(np)
	real*8 :: nu
	integer :: i

	real*8 :: ne(np)
	real*8 :: m_e
	real*8 :: m_bg

	m_e  = molmass(especnum)
	m_bg = molmass(bgspecnum)

	ne = numden(:,especnum)
	
	do i=1,np
		nu = getelectroncollisionfrequency(Te(i)*Tescale,T_gas,P_gas)
		elastic_col(i)=&
		3.d0*ne(i)*(Te(i)-(T_gas/Tescale))*(m_e/m_bg)*nu
	enddo
	
end subroutine elastic_colterm
!=============================================================================
subroutine inelastic_colterm(inelastic_col)
	
	real*8, intent(out) :: inelastic_col(np)

	integer :: i,j
	real*8 :: specarray(nspecies)
	real*8 :: inelterm

	do i=1,np
	
      	  	do j=1,nspecies
  			specarray(j)=numden(i,j)*nscale
	 	 enddo	
	  
	  	call getelectroninelasticterm(Te(i)*Tescale,T_gas,&
			  specarray,inelterm,efield(i))

	  	inelastic_col(i) = inelterm/nscale
	  
	  enddo 

end subroutine inelastic_colterm
!=============================================================================
subroutine timestepping()

	integer :: it,j,rescounter,nresiduals,fnum
	real*8 :: time
	logical :: printflag
	real*8,allocatable :: residual(:)
	
	fnum=13

	nresiduals=1+1+1+no_of_ions+no_of_neutrals
	allocate(residual(nresiduals))
	open(unit=fnum,file='residual.dat')
	
	it = 0
	time = 0.d0
	
	printflag=.false.
	do while((it .lt. maxtimesteps) .and. (time .lt. finaltime))
		
		rescounter=1
		if(mod(it,printfileit) .eq. 0) then
			call printfile(it)
		endif
		
		if(mod(it,printit) .eq. 0) then
			printflag=.true.
		endif
		
		call potentialsolve(printflag,residual(rescounter))
		rescounter=rescounter+1
		
		call elecenergysolve(printflag,residual(rescounter))
		rescounter=rescounter+1
		
		call electronsolve(printflag,residual(rescounter))
		rescounter=rescounter+1

		do j=ionspecmin,ionspecmax
			call ionsolve(j,printflag,residual(rescounter))
			rescounter=rescounter+1
		enddo

		call compute_efield()
		
		do j=neutralspecmin,neutralspecmax
			call neutralsolve(j,printflag,residual(rescounter))
			rescounter=rescounter+1
		enddo

		it=it+1
		time=time+dt
		if(printflag .eqv. .true.) then
			print *,"it=",it,"time=",time
			print *,"==============================&
			 ==============================&
			 =============================="
			 write(fnum,*) residual
		endif
		printflag = .false.

	enddo

	call printfile(it)

	close(fnum)

end subroutine timestepping
!=============================================================================
subroutine printfile(it)
	
	integer, intent(in) :: it
	character(LEN=50) :: fname
	character(LEN=7) :: itstr
	integer :: i,j
	real*8 :: jheating(np),elastic_col(np),inelastic_col(np)
	real*8 :: icurr(np),ecurr(np)

	write(itstr,'(I7.7)') it

	fname="soln_"//trim(itstr)//".dat"
        open(unit=5,file=trim(fname))
	
	call jouleheatingterm(jheating)
	call elastic_colterm(elastic_col)
	call inelastic_colterm(inelastic_col)
	call findconductioncurrent(icurr,ecurr)
      
	do i=1,np
      		write(5,'(E20.10,E20.10,E20.10)',advance='no') (i-1)*dx,&
			phi(i),efield(i)
		do j=1,nspecies
			write(5,'(E20.10)',advance='no') numden(i,j)*nscale
		enddo
		write(5,'(E20.10,E20.10,E20.10,E20.10,E20.10,E20.10,E20.10)',advance='yes') & 
		Ee(i),Te(i),&
		jheating(i),elastic_col(i),inelastic_col(i),ecurr(i),icurr(i)
        enddo

	close(5)

end subroutine printfile
!=============================================================================
function generategaussianvalue(A,b,x0,x) result(val)

		real*8,intent(in) :: A,b,x0,x
		real*8 :: val

		val=A*exp(-b*(x-x0)*(x-x0))

end function generategaussianvalue
!=============================================================================
subroutine findconductioncurrent(ioncurrent,electroncurrent)

	real*8, intent(inout) :: ioncurrent(np)
	real*8, intent(inout) :: electroncurrent(np)
	integer :: i,j
	real*8 :: mu_e,mu_i,D_e,D_i
	real*8 :: grad_ne,grad_ni
	real*8 :: flux

	ioncurrent=0.d0
	electroncurrent=0.d0

	do i=1,np

		mu_e      =   getspecmobility(especnum,Te(i)*Tescale,T_gas,P_gas)
		D_e       =   getspecdcoeff  (especnum,Te(i)*Tescale,T_gas,P_gas)
		
		if(i .eq. 1) then
			grad_ne = (numden(i+1,especnum) - numden(i,especnum))/dx
		else if(i .eq. np) then
			grad_ne = (numden(i,especnum) - numden(i-1,especnum))/dx
		else
			grad_ne = (numden(i+1,especnum) - numden(i-1,especnum))/(2.d0*dx)
		endif
		
		flux = (mu_e*numden(i,especnum)*efield(i) - D_e*grad_ne)*nscale
		electroncurrent(i) = -echarge*flux
		
		do j=ionspecmin,ionspecmax

			mu_i      =   getspecmobility(j,Te(i)*Tescale,T_gas,P_gas)
			D_i       =   getspecdcoeff  (j,Te(i)*Tescale,T_gas,P_gas)
		
			if(i .eq. 1) then
				grad_ni = (numden(i+1,j) - numden(i,j))/dx
			else if(i .eq. np) then
				grad_ni = (numden(i,j) - numden(i-1,j))/dx
			else
				grad_ni = (numden(i+1,j) - numden(i-1,j))/(2.d0*dx)
			endif
		
			flux = (mu_i*numden(i,j)*efield(i) - D_i*grad_ni)*nscale
			ioncurrent(i) = ioncurrent(i) + spec_charge(j)*echarge*flux
		enddo

	enddo

end subroutine findconductioncurrent
!=============================================================================

end module plasma_solver
