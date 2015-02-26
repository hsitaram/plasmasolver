module convdiff

implicit none

integer,parameter :: minsize=5
real*8 :: k
real*8 :: c
real*8 :: h
real*8 :: phi0
real*8 :: phiL

contains
!===============================================================
subroutine findupwindconvflux(cL,cR,uL,uR,flux)

	real*8,intent(in) :: cL,cR
	real*8,intent(in) :: uL,uR
	real*8,intent(out) :: flux

	real*8 :: c_half

	c_half=0.5*(cL+cR)
	
	if(c_half .ge. 0.d0) then
	    flux=c_half*uL
	else
            flux=c_half*uR
	endif
	     
end subroutine findupwindconvflux
!===============================================================
subroutine finddiffusiveflux(DL,DR,uL,uR,dx,flux)

	real*8,intent(in) :: DL,DR
	real*8,intent(in) :: uL,uR
	real*8,intent(in) :: dx
	real*8,intent(out) :: flux

	real*8 :: D_half

	D_half = 0.5*(DL+DR)
	flux = D_half*(uR-uL)/dx

end subroutine finddiffusiveflux
!===============================================================
subroutine findAX(AX,X,timederivfactor,vel,dcoeff,reac,dirc_bc_flags,&
			flux_bc_flags,dircvals,fluxvals,dx,dt,n)

      integer, intent(in) :: n
      real*8, intent(inout) :: AX(n)
      real*8, intent(in)  :: X(n)

      real*8,intent(in)  :: vel(n),dcoeff(n),reac(n)
      real*8,intent(in)  :: dx,dt
      logical,intent(in) :: dirc_bc_flags(2),flux_bc_flags(2)
      real*8, intent(in) :: dircvals(2),fluxvals(2)
      real*8, intent(in) :: timederivfactor
      
      integer :: i
      real*8 :: flux
      real*8 :: dx2

      dx2 = dx*dx

      AX(:) = timederivfactor*X(:)/dt
      
      !is the volume half of dx at the boundary?
      AX(1) = 0.5*timederivfactor*X(1)/dt
      AX(n) = 0.5*timederivfactor*X(n)/dt

      !convection and diffusion terms
      do i=1,n-1
      	call findupwindconvflux(vel(i),vel(i+1),X(i),X(i+1),flux)
	AX(i)   = AX(i)   + flux/dx
	AX(i+1) = AX(i+1) - flux/dx
	call finddiffusiveflux(dcoeff(i),dcoeff(i+1),X(i),X(i+1),dx,flux)
	AX(i)   = AX(i)   - flux/dx
	AX(i+1) = AX(i+1) + flux/dx
      enddo

      !print *,"AX inside findAX:",AX

      !reaction term
      do i=2,n-1
      	AX(i) = AX(i) - reac(i)*X(i)
      enddo
      AX(1) = AX(1) - 0.5*reac(i)*X(i)
      AX(n) = AX(n) - 0.5*reac(i)*X(i)

      !boundary conditions
      if(dirc_bc_flags(1) .eqv. .true.) then
	      AX(1)=X(1)
      endif
      if(dirc_bc_flags(2) .eqv. .true.) then
	     AX(n)=X(n)
      endif

      !note: for flux bc, the terms go into b and not AX
      
      !if the volume is half of dx at the boundaries the reac(i)*X(i)
      !should be 0.5*reac(i)*X(i) 


end subroutine findAX
!===============================================================
subroutine findrhs(b,xold,timederivfactor,&
		source,dirc_bc_flags,flux_bc_flags,dircvals,&
				fluxvals,dx,dt,n)

      integer, intent(in) :: n
      real*8, intent(inout) :: b(n)
      real*8, intent(in) :: xold(n),source(n)
      logical, intent(in) :: dirc_bc_flags(2),flux_bc_flags(2)
      real*8, intent(in) :: dircvals(2),fluxvals(2)
      real*8, intent(in) :: dx,dt
      real*8, intent(in) :: timederivfactor
      integer :: i

      do i=1,n
      	b(i) = timederivfactor*xold(i)/dt + source(i)
      enddo	
      !is the volume half of dx at the boundary?
      b(1) = 0.5*(timederivfactor*xold(1)/dt + source(1))
      b(n) = 0.5*(timederivfactor*xold(n)/dt + source(n))

      if(dirc_bc_flags(1) .eqv. .true.) then
	      b(1) = dircvals(1)
      endif
      if(dirc_bc_flags(2) .eqv. .true.) then
	      b(n) = dircvals(2)
      endif
      if(flux_bc_flags(1) .eqv. .true.) then
	      b(1) = b(1) + fluxvals(1)/(dx)
      endif
      if(flux_bc_flags(2) .eqv. .true.) then
	      b(n) = b(n) - fluxvals(2)/(dx)
      endif

      !if the volume is half of dx then the fluxvals should be
      !divided by 0.5*dx

end subroutine findrhs
!===============================================================
subroutine noprecond(MinvX,X,timederivfactor,vel,dcoeff,reac,dirc_bc_flags,&
				flux_bc_flags,dircvals,fluxvals,&
				dx,dt,n)
      integer, intent(in) :: n
      real*8, intent(inout) :: MinvX(n)
      real*8, intent(in)  :: X(n)
	
      logical, intent(in) :: dirc_bc_flags(2),flux_bc_flags(2)
      real*8, intent(in) :: dircvals(2),fluxvals(2)
      real*8,intent(in)  :: vel(n),dcoeff(n),reac(n)
      real*8,intent(in)  :: dx,dt
      real*8,intent(in)  :: timederivfactor

      MinvX = X

end subroutine noprecond
!===============================================================
subroutine gauss_seidel_smoothing(res,b,X,timederivfactor,vel,dcoeff,reac,&
		dirc_bc_flags,flux_bc_flags,dircvals,&
				fluxvals,dx,dt,n,maxiter)

	integer, intent(in)   :: n
	real*8, intent(in)    :: b(n)
	real*8, intent(inout) :: X(n)
	real*8, intent(inout)   :: res(n)
	integer, intent(in) :: maxiter

        logical, intent(in) :: dirc_bc_flags(2),flux_bc_flags(2)
        real*8, intent(in) :: dircvals(2),fluxvals(2)
        real*8,intent(in)  :: vel(n),dcoeff(n),reac(n)
        real*8,intent(in)  :: dx,dt
	real*8,intent(in)  :: timederivfactor

	integer :: i,it

	real*8 :: diag
	real*8 :: cL,cR,chalf
	real*8 :: dL,dR,dhalf
	real*8 :: offdiag
	real*8 :: AX(n)
	real*8 :: dx2

	dx2 = dx*dx

	!gauss seidel iterations
	do it=1,maxiter
	     
	     do i=1,n

	        diag = -reac(i) + timederivfactor*1.d0/dt

		!is the volume half of dx at the boundary?
		if((i .eq. 1) .or. (i .eq. n)) then
			diag = 0.5*diag
		else
			diag = 0.5*diag
		endif

		offdiag=0.d0

		!right face
		if(i .lt. n) then
			!convection term
			cL = vel(i)
			cR = vel(i+1)
			chalf = 0.5*(cL+cR)
			if(chalf .ge. 0) then
				diag = diag + chalf/dx
			else
				offdiag = offdiag + chalf*X(i+1)/dx
			endif
		
			!diffusion term
			dL = dcoeff(i)
			dR = dcoeff(i+1)
			dhalf = 0.5*(dL + dR)
		
			diag = diag + dhalf/dx2
			offdiag = offdiag - X(i+1)*dhalf/dx2
		endif
	     	
	     	!left face
		if(i .gt. 1) then

			!convection term
			cL = vel(i-1)
			cR = vel(i)
			chalf = 0.5*(cL+cR)
			if(chalf .ge. 0) then
				offdiag = offdiag - chalf*X(i-1)/dx
			else
				diag = diag - chalf/dx
			endif
		
			!diffusion term
			dL = dcoeff(i-1)
			dR = dcoeff(i)
			dhalf = 0.5*(dL+dR)
		
			diag = diag + dhalf/dx2
			offdiag = offdiag - X(i-1)*dhalf/dx2
		endif

		if(i .eq. 1) then
			if(dirc_bc_flags(1) .eqv. .true.) then
				diag    = 1.d0
				offdiag = 0.d0
			endif
		endif	

		if(i .eq. n) then
			if(dirc_bc_flags(2) .eqv. .true.) then
				diag    = 1.d0
				offdiag = 0.d0
			endif
		endif

		X(i) = (b(i)-offdiag)/diag

		enddo
	enddo

	call  findAX(AX,X,timederivfactor,vel,dcoeff,reac,dirc_bc_flags,&
			flux_bc_flags,dircvals,fluxvals,dx,dt,n)

	res = b - AX 

end subroutine gauss_seidel_smoothing
!===============================================================
subroutine restriction(Xh,X2h,n2h) 

	integer, intent(in) :: n2h
	real*8, intent(inout) :: Xh(2*n2h-1)
	real*8, intent(inout) :: X2h(n2h)

	integer :: i

	X2h(1)     = Xh(1)
	X2h(n2h) = Xh(2*n2h-1)

	do i=2,n2h-1
	   X2h(i)=0.25*(Xh(2*i-2) + 2.0*Xh(2*i-1) + Xh(2*i))
	enddo


end subroutine restriction
!===============================================================
subroutine prolong(Xh,X2h,n2h)

	integer, intent(in) :: n2h
	real*8, intent(inout) :: Xh(2*n2h-1)
	real*8, intent(inout) :: X2h(n2h)

	integer :: i

	do i=1,n2h-1
	   Xh(2*i-1) = X2h(i)
	   Xh(2*i)   = 0.5*(X2h(i) + X2h(i+1))
	enddo

	Xh(2*n2h-1) = X2h(n2h)
	
end subroutine prolong
!===============================================================
recursive subroutine dovcycle(X,b,timederivfactor,vel,dcoeff,reac,&
		dirc_bc_flags,flux_bc_flags,dircvals,&
				fluxvals,dx,dt,n)

	integer, intent(in)   :: n
	real*8, intent(inout) :: X(n)
	real*8, intent(inout) :: b(n)
        
	logical, intent(in)  :: dirc_bc_flags(2),flux_bc_flags(2)
        real*8, intent(in)   :: dircvals(2),fluxvals(2)
        real*8, intent(in)   :: vel(n),dcoeff(n),reac(n)
        real*8 ,intent(in)   :: dx,dt
	real*8 ,intent(in)   :: timederivfactor

	real*8  :: resh(n)
	real*8 :: eh(n)
	real*8 :: velh(n)
	real*8 :: dcoeffh(n)
	real*8 :: reach(n)
	
	real*8 :: res2h(n/2+1)
	real*8 :: e2h(n/2+1)
	real*8 :: vel2h(n/2+1)
	real*8 :: dcoeff2h(n/2+1)
	real*8 :: reac2h(n/2+1)

	velh    = vel
	dcoeffh = dcoeff
	reach   = reac

	if(n .le. minsize) then
		call gauss_seidel_smoothing(resh,b,X,timederivfactor,velh,dcoeffh,reach,&
				dirc_bc_flags,flux_bc_flags,dircvals,&
				fluxvals,dx,dt,n,1000)
	else

		!initial smoothing
		call gauss_seidel_smoothing(resh,b,X,timederivfactor,velh,dcoeffh,reach,&
				dirc_bc_flags,flux_bc_flags,dircvals,&
				fluxvals,dx,dt,n,4)

		!restriction of residual from fine to coarse grid
         	call restriction(resh,res2h,n/2+1)
         	
		call restriction(reach,reac2h,n/2+1)
         	call restriction(dcoeffh,dcoeff2h,n/2+1)
         	call restriction(velh,vel2h,n/2+1)

		e2h = 0.d0
		call dovcycle(e2h,res2h,timederivfactor,vel2h,dcoeff2h,reac2h,&
				dirc_bc_flags,flux_bc_flags,dircvals,&
				fluxvals,2*dx,dt,n/2+1)

		!prolong error from coarse to fine grid
		call prolong(eh,e2h,n/2+1)

		!update
		X(:) = X(:) + eh(:)

		!post smooth
		call gauss_seidel_smoothing(resh,b,X,timederivfactor,velh,dcoeffh,reach,&
				dirc_bc_flags,flux_bc_flags,dircvals,&
				fluxvals,dx,dt,n,4)
	endif

end subroutine dovcycle
!===============================================================
subroutine mgridprecond(MinvX,X,timederivfactor,vel,dcoeff,reac,dirc_bc_flags,&
				flux_bc_flags,dircvals,fluxvals,&
				dx,dt,n)
	
	integer, intent(in)    :: n
	real*8, intent(inout)    :: MinvX(n)
	real*8, intent(inout)  :: X(n)
	
	logical, intent(in) :: dirc_bc_flags(2),flux_bc_flags(2)
        real*8, intent(in) :: dircvals(2),fluxvals(2)
        real*8,intent(in)  :: vel(n),dcoeff(n),reac(n)
        real*8,intent(in)  :: dx,dt
	real*8,intent(in)  :: timederivfactor
	
	integer :: nvcycles,i
	!real*8 :: AX(n),res(n)
	!real*8 :: norm
	!integer :: j

	!norm=0.d0
	nvcycles=5

	do i=1,nvcycles
		call dovcycle(MinvX,X,timederivfactor,vel,dcoeff,reac,&
		dirc_bc_flags,flux_bc_flags,dircvals,&
				fluxvals,dx,dt,n)
	enddo
	

end subroutine mgridprecond
!===============================================================

end module convdiff
