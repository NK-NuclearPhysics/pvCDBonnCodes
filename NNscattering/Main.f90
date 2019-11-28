Program Main
    ! Calculate the NN scattering phase shif with one-boson exchange potential
    ! Based on matrix inversion method providen by
    !       M. Haftel and F. Tabakin, Nucl. Phys. A, 158, 1 (1970)
    !   
    ! Codes written by Chencan Wang
    ! 2019.10.18
    use PhaseShift
    implicit none
    
    integer, parameter  :: N=11
    real*8              :: Elab(N), Klab, delta(6), Mn

    integer             :: ia, ib
    real*8,  parameter  :: Mpp=938.2720813d0,&  ! Reduced mass    
                           Mnp=938.9183019d0,&  ! corresponding to inn=1(pp)
                           Mnn=939.5654133d0,&  ! inn=2(np) & inn=3(nn)
                           Pi=3.1415926535898d0,&
                           deg=180.d0/Pi
    character           :: cSymbols(0:6)
    
    
    data Elab/1, 5, 10, 25, 50, 100, 150, 200, 250, 300, 350/
    data cSymbols/'S','P','D','F','G','H','I'/
    j = 1
    inn = 2
    heform = .false.
    sing   = .true.
    trip   = .true.
    coup   = .true.   
    ctype ='C'
 
    
    open(unit=1,file='PhaseShift.txt')
    write(1,10)
    selectcase(j)
    case(0)
        write(1,15) 
    case default
        write(1,12) cSymbols(j)  ,j,  cSymbols(j),  j, cSymbols(j-1),j,j,cSymbols(j+1),j 
    end select
    selectcase(inn)
    case(1)
        Mn = Mpp
    case(2)
        Mn = Mnp        
    case(3)
        Mn = Mnn
    end select 
    
    do ia = 1,N
        selectcase(inn)
        case(2)
                Klab = Mpp*sqrt(Elab(ia)*(Elab(ia)+2*Mnn)/((Mpp+Mnn)**2+2*Elab(ia)*Mpp))
        case default
                Klab = sqrt(Mn*Elab(ia)/2.0)
        end select 
 
        call GetPhase(Klab,delta)
        if(inn==2.and.Elab(ia)<=25.and.delta(4)<0) then
            delta(5) = -delta(5)
            delta(4) = delta(4)+180.d0
        end if
        write(1,20), Elab(ia), delta(1), delta(2), delta(4),&
                               delta(5), delta(3)
        
    end do
    
    
    
    
    close(1)
10  format(/'Phase Shifts[deg]'/)
12  format('Elab[MeV]',5x,'1',a,i1,9x,'3',a,i1,    &
                       9x,'3',a,i1,10x,'e',i1,9x,'3',a,i1) 
15  format('Elab[MeV]',5x,'1S0',9x,'Non',9x,'Non', &
              9x,'Non',9x,'3P0')   
20  format(1x, f5.1, 5x, f7.2, 5x, f7.2,  &
           5x, f7.2, 5x, f7.2, 5x, f7.2   )  
30  format(//)    
    end Program