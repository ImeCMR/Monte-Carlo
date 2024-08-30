program armc
implicit none
!
! Defining variables
!
integer::natom,nsteps,nprint,i,j,k,nacc=0,ncont,ilts,istep,idd,istart
  real,dimension(500)::x,y,z 
  real::dseed,temp,epsln,sigma,Penergy,Tenergy,expTenergy,xbox,ybox,zbox,mar,kb,rcut,dx,dy,dz,dd,xx,yy, &
        zz,epsln4,sigma2,rcut2,itemp,scal,reqKenergy,avgPenergy,&
        avgPenergy2,sdPenergy,accr,potEnergy,dlt,dlte,oldPenergy,bf
  character(len=20)::gout,eout,xout,ipdb,fpdb
!
  write(6,*)'==================================='
  write(6,*)''
  write(6,*)'     Calculation is in progess     '
  write(6,*)''
  write(6,*)'==================================='
! Open input file
!
  open(unit=10,file='armc.dat',status='old')
!
! Read data from the file
!
  read(10,*) istart             ! CDT generate (0) or from the last MC run (1)
  read(10,*) ilts               ! rigid (0) OR random (1)
  read(10,*) natom              ! No. of atoms in the system
  read(10,*) nsteps             ! No. of MC steps
  read(10,*) nprint             ! printing interval
  read(10,*) epsln              ! epsilon
  read(10,*) sigma              ! sigma
  read(10,*) temp               ! Initial temperature
  read(10,*) expTenergy         ! Expected Total Energy (kJ/mol)
  read(10,*) mar                ! Mass of Argon (g/mol)
  read(10,*) kb                 ! Boltzmann Constant (kJ/mol/K)
  read(10,*) rcut               ! Cutoff distance
  read(10,*) dlt                ! MC step size
  read(10,*) xbox               !
  read(10,*) ybox               ! Simulation box dimensions
  read(10,*) zbox               !
  read(10,'(a)')gout            ! Name of the general output
  read(10,'(a)')eout            ! Energy out file
  read(10,'(a)')xout            ! Coordinate output file(trajectory file)
  read(10,'(a)')ipdb            ! Initial CDT file in PDB format)
  read(10,'(a)') fpdb               ! Fianl CDT file in PDB format
!
! Create output files
!
  open(unit=20,file=gout,status='unknown')
  open(unit=21,file=eout,status='unknown')
  open(unit=22,file=xout,status='unknown')
  open(unit=23,file=ipdb,status='unknown')
  open(unit=24,file=fpdb,status='unknown')
!
! write data
!
  write(20,500) istart,ilts,natom,nsteps,nprint,epsln,sigma,temp,expTenergy,mar,kb,rcut,xbox,ybox,zbox,gout,eout, &
                xout,ipdb,fpdb
 500 format(/,2x,' Start or restart option     ='i5,/,&
              2x,' Initial CDT option          ='i5,/,&
              2x,' Number of atoms             =',i5,/,&
              2x,' Number of MC steps          =',i7,/,&
              2x,' Printing interval           =',i5,/,&
              2x,' Ar Epsilon                  =',f8.4,' kJ/mol',/,&
              2x,' Ar Sigma                    =',f8.4,' nm',/,&
              2x,' Initial temperature         =',f8.4,' K',/,&
              2x,' Expected Total Energy       =',f9.2,' kJ/mol',/,&
              2x,' Mass of Ar                  =',f8.4,' g/mol',/,&
              2x,' Boltzmann Constant          =',f8.4,' kJ/mol/K',/,&
              2x,' Cutoff distance             =',f8.4,' nm',/,&
              2x,' xbox                        =',f8.4,' nm',/,&
              2x,' ybox                        =',f8.4,' nm',/,&
              2x,' zbox                        =',f8.4,' nm',/,&
              2x,' General output file         =',a20,/,&
              2x,' Energy output file          =',a20,/,&
              2x,' Coordinate out file         =',a20,/,&
              2x,' Initial Coordinates         =',a20,/,&
              2x,' Final Coordinates           =',a20,/)
!
 if (istart == 0) then   ! for start option
!
! -----------Generating initial coordinates
!
 if (ilts == 0) then
   ncont =0
   do i = 1,6
     xx = (i-1)*0.39
     do j = 1,6
       yy = (j-1)*0.39
       do k = 1,6
         zz = (k-1)*0.39
         ncont = ncont + 1
         x(ncont) = xx
         y(ncont) = yy
         z(ncont) = zz
       end do
     end do
   end do
!
 else
!
   ncont = 0
99 do i = 1,5000
     xx = (rand() - 0.5) * xbox
     yy = (rand() - 0.5) * ybox
     zz = (rand() - 0.5) * zbox
     if (ncont == 0) then
       ncont = ncont + 1
       x(ncont) = xx
       y(ncont) = yy
       z(ncont) = zz
     else
       do j = 1,ncont
         dx = x(j) - xx
         dy = y(j) - yy
         dz = z(j) - zz
         dd = sqrt(dx*dx + dy*dy + dz*dz)
         if (dd<0.38) then
           go to 99
         end if
       end do
       ncont = ncont + 1
       x(ncont) = xx
       y(ncont) = yy
       z(ncont) = zz
      !if (mod(ncont,25) == 0) write(6,'(i5)') ncont
       if (ncont == natom) go to 88
     end if
   end do
   88 continue
 end if
!
 else    ! for start option
!
   open(unit=11, file='finl.pdb',status='old')
   read(11,*)
   read(11,*)
   do i = 1,natom
     read(11,'(30x,3f8.3)') xx,yy,zz
     x(i) = xx*0.10 
     y(i) = yy*0.10  ! convert anstrons into nm 
     z(i) = zz*0.10 
   end do
!
 end if
!
!-------Initial CDT file
!
  do i = 1, natom
    xx = x(i)*10
    yy = y(i)*10  ! 10 is used to convert into anstrons
    zz = z(i)*10
    write(23,501) i,'Ar  ','Ar  ',i,xx,yy,zz,0.0,0.0
  end do
 501 format('ATOM',2x,i5,2x,a4,a4,1x,i4,4x,3f8.3,2f6.2)
!
!--------- Initial potential energy
!
  epsln4    = 4.0*epsln
  sigma2    = sigma*sigma
  rcut2     = rcut*rcut
  potEnergy = 0.0
!
  do i = 1,natom
    call eng(i,natom,mar,epsln4,sigma2,rcut2,x,y,z,xbox,ybox,zbox,Penergy)
    potEnergy = potEnergy + Penergy
  end do
  potEnergy = 0.5*potEnergy   ! Because same calculation is done twise (like 1,4 and 4,1 both are calculated)
!
!---------Initial potential energy and acceleration
!
  write(20,*)  'Initial Potential Energy =', potEnergy,' kJ/mol'
  write(20,*) ''
!
  idd         = 0
  avgPenergy  = 0.0
  avgPenergy2 = 0.0
!
!-------MC run starts here
!
  do istep = 1,nsteps
!
!--------- Pick an atom randomly
!
    idd = int(rand()*natom) + 1
!
!--------- Potential Energy before the movement
!
    call eng(idd,natom,mar,epsln4,sigma2,rcut2,x,y,z,xbox,ybox,zbox,Penergy)
    oldPenergy = Penergy
!
!--------- Move the atom 'idd' randomly
!
    xx = x(idd)
    yy = y(idd)   ! If the move is rejected we have oriinal coordinates are stored here for any usage
    zz = z(idd)
!
    x(idd) = x(idd) + dlt*(rand()-0.5)      ! getting positive or negative values
    y(idd) = y(idd) + dlt*(rand()-0.5) 
    z(idd) = z(idd) + dlt*(rand()-0.5) 
!
!--------- Periodic boundry conditions
!
     x(idd) = x(idd) - anint(x(i)/xbox)*xbox
     y(idd) = y(idd) - anint(y(i)/ybox)*ybox
     z(idd) = z(idd) - anint(z(i)/zbox)*zbox
!
!-----------Call energy
!
     call eng(idd,natom,mar,epsln4,sigma2,rcut2,x,y,z,xbox,ybox,zbox,Penergy)
     dlte = Penergy - oldPenergy
!
     bf = exp(-dlte/(kb*temp))
     if (bf > rand()) then
       potEnergy = potEnergy + dlte
       nacc      = nacc + 1
     else
       x(idd) = xx
       y(idd) = yy
       z(idd) = zz
     end if
!
! Calculating the averages
! 
     avgPenergy  = avgPenergy + potEnergy
     avgPenergy2 = avgPenergy2 + potEnergy*potEnergy
!        
     if(mod(istep,nprint) == 0) then
       accr = real(nacc)/istep
     ! write(6,'(a,i5)') 'Number of steps = ',istep
       idd = idd + 1
       write(21,'(i7,2f13.4)') istep,potEnergy,accr
       write(22,'(a,i8)')'MODEL',idd
       write(22,'(a,3f8.3,3f6.1)') 'REMARKS',xbox*10,ybox*10,zbox*10,90.0,90.0,90.0
       do i = 1,natom
         xx = x(i)*10.0
         yy = y(i)*10.0 ! write in anstron unit
         zz = z(i)*10.0
         write(22,501) i,'Ar  ','Ar  ',i,xx,yy,zz,0.0,0.0
       end do
       write(22,'(a)')'ENDMDL'
       write(20,503) istep,potEnergy,accr
     end if
     if(mod(istep,2000) == 0)write(6,'(a,i8)')' No. of steps = ',istep
  end do
  write(6,*)'=========================================='
  write(6,*)''
  write(6,*)'        Calculation is done !             '
  write(6,*)'    armc.out  ! General output file       '
  write(6,*)'    armc.eng  ! Energy output file        '
  write(6,*)'    armc.cdt  ! Coordinate output file    '
  write(6,*)'    init.pdb  ! Coordinate output file    '
  write(6,*)'    finl.pdb  ! Coordinate output file    '
  write(6,*)''
  write(6,*)'=========================================='
!
! -----------------MC run is over now
!
  503 format(/,2x,'No. of steps                  = ',i7,/,&
               2x,'Potential Energy              = ',f12.4,' kJ/mol',/,&
               2x,'Acc. Ratio                    = ',f12.4,/)
!
! ----------------Final configuration
!
    write(24,'(a)') 'Final configuration'
    write(24,'(a,3f8.3,3f6.1)') 'REMARKS',xbox*10,ybox*10,zbox*10,90.0,90.0,90.0
    do i = 1,natom
      xx = x(i)*10.0
      yy = y(i)*10.0 ! write in anstron unit
      zz = z(i)*10.0
      write(24,501) i,'Ar  ','Ar  ',i,xx,yy,zz,0.0,0.0
    end do
    write(24,'(a)')'ENDMDL'
    close(unit=24)        ! if write again it start from the top
!
! ----------------Calculating averages
!
     sdPenergy    = sqrt((nsteps*avgPenergy2 - avgPenergy*avgPenergy)/nsteps**2)
     avgPenergy   = avgPenergy/nsteps
     write(20,505) avgPenergy,sdPenergy,accr
505  format(//,2x,'Average properties -----------------',//,&
               2x,'Potential energy  = ',f12.4,' +/- ',f9.4,' kJ/mol',/,&
               2x,'Acc. Ratio        = ',f12.4,/)
!
end program armc
!
! ***************************************************************
!
! Subroutine
!
  subroutine eng(idd,natom,mar,epsln4,sigma2,rcut2,x,y,z,xbox,ybox,zbox,Penergy)
  implicit none
!
  integer::i,j,k,natom,ii,idd
  real,dimension(500)::x,y,z
  real::epsln4,sigma2,rcut2,xbox,ybox,zbox,Penergy,cc,c2,c6,c12,r2,xi,yi,zi,dx,dy,dz,mar
!
! ------Initialize
!
  Penergy = 0.0
!
    xi = x(idd)
    yi = y(idd)
    zi = z(idd)
    do j = 1,natom     ! To save computer time :: Parallel computation
      if (j /= idd) then
        dx = x(j) - xi
        dy = y(j) - yi
        dz = z(j) - zi
! 
!------Minimum image convention
!
        dx = dx - xbox*anint(dx/xbox)
        dy = dy - ybox*anint(dy/ybox)
        dz = dz - zbox*anint(dz/zbox)
!
        r2 = dx*dx + dy*dy + dz*dz   ! Calculate distance between two atoms
        if (r2 <= rcut2) then
          c2  = sigma2/r2
          c6  = c2*c2*c2
          c12 = c6*c6
          Penergy  = Penergy + epsln4*(c12 - c6)
        end if
      end if 
    end do
!
  end subroutine
