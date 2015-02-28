program chem_driver

      use chem_module
	implicit none
	
      real*8 :: specdenvec(nspecies);
      real*8 :: specprod,inelterm

      print *,"nspecies:",nspecies

      specdenvec=1.d12

      specdenvec(1)=3.1d22
      specdenvec(2)=1.d14
      specdenvec(3)=1.d14

      call initializechemistry()
      call getspecproduction(2,11604.d0,300.d0,specdenvec,specprod,2.d7);
      call getelectroninelasticterm(11604.d0,300.d0,specdenvec,inelterm,2.d7);

      write(*,'(A E20.10)')"specprod:",specprod
      write(*,'(A E20.10)')"inelterm:",inelterm

end program chem_driver
