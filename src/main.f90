program Laffer

    use Utilities
    use Model_Parameters
    use PolicyFunctions
    use Tauchen
    use hybrd_wrapper, only: setHybrParams
    use SolveRetired_EGM

    implicit none
    integer :: ik,tprint,it2,it3,it4,ium,iam,iuf,iaf,ix,j,iu2,ik2,ifc,counter
    real(8) :: dum,dum2,dum3
    integer :: i_k, i
    EXTERNAL labor2
    EXTERNAL labor1
    EXTERNAL labor3
    EXTERNAL labors
    !print *, "hello"
    !call TestLinInterp 
    !call OMP_SET_NUM_THREADS(80)

    call Initialize
    call Initialize_Retired()
    call setHybrParams(2)

    !Compute optimal policies in retirement

    do while(abs(epsilon4)>0.01d0)

        epsilon=1d0
        epsilon2=1d0
        epsilon3=1d0
        epsilon5=1d0

        iter=0

        do while((abs(epsilon)>0.002d0).OR.(abs(epsilon3)>0.001d0).OR.(abs(epsilon5)>0.001d0).OR.(abs(epsilon6)>0.01d0))

            iter=iter+1
            
            dum2= 1.5d0*wage(1,a(1,5),dble(T),u(1,5))/(1d0+t_employer)
            call MakeGrid(nw,wage_grid,0.01d0,dum2,2d0)

            !$OMP PARALLEL PRIVATE(ik)
            !$OMP DO SCHEDULE(DYNAMIC)
            do ik=1,nc
                call lsupply(ik)
            end do
            !$OMP END DO    
            !$OMP END PARALLEL

            do it=1,Tret

                !write(*, "(1x,' Solving retired,age ',i0,A,$)") it, CHAR(13)

                !$OMP PARALLEL PRIVATE(ik)
                !$OMP DO SCHEDULE(DYNAMIC)
                do ik=1,nk 
                    call SolveInRetirement(ik)
                end do
                !$OMP END DO    
                !$OMP END PARALLEL
                
                !$OMP PARALLEL PRIVATE(counter)
                !$OMP DO SCHEDULE(DYNAMIC)
                    do counter=1,nk*na*nu
                        call SolveInRetirement2(counter)
                    end do
                !$OMP END DO    
                !$OMP END PARALLEL

                !$OMP PARALLEL PRIVATE(counter)
                !$OMP DO SCHEDULE(DYNAMIC)
                    do counter=1,nk*na*nu
                        call SolveInRetirement3(counter)
                    end do
                !$OMP END DO    
                !$OMP END PARALLEL
                
                !$OMP PARALLEL PRIVATE(counter)
                !$OMP DO SCHEDULE(DYNAMIC)
                    do counter=1,nk*na*nu
                        call partest9(counter)
                    end do
                !$OMP END DO    
                !$OMP END PARALLEL

                ev_ret(1,1,:,:,:,:,:,:,:,Tret-it+1,:,:)=v_ret(1,1,:,:,:,:,:,:,:,Tret-it+1,:,:)
                evs_ret(1,:,:,:,:,:,Tret-it+1,:)=vs_ret(1,:,:,:,:,:,Tret-it+1,:)

                call Compute_EV_DEV(Tret-it+1)

                stop
                
            end do


            !Compute optimal policies at age 64

            !$OMP PARALLEL PRIVATE(counter)
            !$OMP DO SCHEDULE(DYNAMIC)
            do counter=1,nk*na*nu
                call Solvelastactive(counter)
            end do
            !$OMP END DO    
            !$OMP END PARALLEL

            
            it=0
            !Print *,'t is',T-it

            !$OMP PARALLEL PRIVATE(counter)
            !$OMP DO SCHEDULE(DYNAMIC)
            do counter=1,nk*na*nu
                call partest8(counter)
            end do
            !$OMP END DO    
            !$OMP END PARALLEL


            !$OMP PARALLEL PRIVATE(counter)
            !$OMP DO SCHEDULE(DYNAMIC)
            do counter=1,nk*na*nu
                call partest7(counter)
            end do
            !$OMP END DO    
            !$OMP END PARALLEL


            !! Compute spline coefficients:
            !
            !!$OMP PARALLEL PRIVATE(counter)
            !!$OMP DO SCHEDULE(DYNAMIC)
            !do counter = 1, nu*na*nfc
            !    call partest3(counter)
            !end do
            !!$OMP END DO    
            !!$OMP END PARALLEL

            ! Compute spline coefficients:


            !$OMP PARALLEL PRIVATE(counter)
            !$OMP DO SCHEDULE(DYNAMIC)
            do counter=1,nk*na*nu
                call partest(counter)
            end do
            !$OMP END DO    
            !$OMP END PARALLEL


            !$OMP PARALLEL PRIVATE(counter)
            !$OMP DO SCHEDULE(DYNAMIC)
            do counter=1,nk*na*nu
                call partest2(counter)
            end do
            !$OMP END DO    
            !$OMP END PARALLEL

            !Compute optimal policies for age 2-63
            do it=1,T-2

                !Print *,'t is',T-it

                !$OMP PARALLEL PRIVATE(counter)
                !$OMP DO SCHEDULE(DYNAMIC)
                do counter=1,nk*na*nu
                    call SolveActiveLife(counter)
                end do
                !$OMP END DO    
                !$OMP END PARALLEL


                !$OMP PARALLEL PRIVATE(counter)
                !$OMP DO SCHEDULE(DYNAMIC)
                do counter=1,nk*na*nu
                    call partest8(counter)
                end do
                !$OMP END DO    
                !$OMP END PARALLEL


                !$OMP PARALLEL PRIVATE(counter)
                !$OMP DO SCHEDULE(DYNAMIC)
                do counter = 1, nk*na*nu
                    call partest7(counter)
                end do
                !$OMP END DO    
                !$OMP END PARALLEL   

                !! Compute spline coefficients:
                !
                !!$OMP PARALLEL PRIVATE(counter)
                !!$OMP DO SCHEDULE(DYNAMIC)
                !do counter = 1, nfc*na*nu
                !    call partest3(counter)
                !end do
                !!$OMP END DO    
                !!$OMP END PARALLEL



                !$OMP PARALLEL PRIVATE(counter)
                !$OMP DO SCHEDULE(DYNAMIC)
                do counter=1,nk*na*nu
                    call partest(counter)
                end do
                !$OMP END DO    
                !$OMP END PARALLEL


                !$OMP PARALLEL PRIVATE(counter)
                !$OMP DO SCHEDULE(DYNAMIC)
                do counter=1,nk*na*nu
                    call partest2(counter)
                end do
                !$OMP END DO    
                !$OMP END PARALLEL

            end do


            it=T-1               

            !Print *,'t is',T-it

            !$OMP PARALLEL PRIVATE(counter)
            !$OMP DO SCHEDULE(DYNAMIC)
            do counter=1,nk*na*nu
                call Solvefirstactive(counter)
            end do
            !$OMP END DO    
            !$OMP END PARALLEL

            
               ! tprint=T
               ! 
               ! open(1,file='csingle.txt')
               !     open(2,file='ksingle.txt')
               !     open(3,file='nsingle.txt')
               !     open(4,file='vsingle.txt')
               !     open(5,file='vsinglem.txt')
               !     do it2 = 1, nk
               !         write(1,'(4f12.6)'), k_grid(it2), cs(1,it2,5,3,3,tprint,8),cs(2,it2,5,3,3,tprint,8),cs(2,it2,5,5,5,tprint,8)
               !         write(2,'(4f12.6)'), k_grid(it2), ks(1,it2,5,3,3,tprint,8),ks(2,it2,5,3,3,tprint,8),ks(2,it2,5,5,5,tprint,8)
               !         write(3,'(4f12.6)'), k_grid(it2), ns(1,it2,5,3,3,tprint,8),ns(2,it2,5,3,3,tprint,8),ns(2,it2,5,5,5,tprint,8)
               !         write(4,'(4f12.6)'), k_grid(it2), evs(1,it2,5,3,3,tprint,8),evs(2,it2,5,3,3,tprint,8),evs(2,it2,5,5,5,tprint,8)
               !         write(5,'(4f12.6)'), k_grid(it2), evm(1,it2,5,3,3,tprint,8),evm(2,it2,5,3,3,tprint,8),evm(2,it2,5,5,5,tprint,8)
               !     end do
               !     close(1)
               !     close(2)
               !     close(3)
               !     close(4)
               !     close(5)
               ! 
               !     open(1,file='ctest.txt')
               !         open(2,file='ktest.txt')
               !         open(3,file='vtest.txt')
               !         open(4,file='nmtest.txt')
               !         open(5,file='nftest.txt')
               !         !open(5,file='evtest.txt')
               !         !open(4,file='ve_vu.txt')
               !         !open(4,file='tbc.txt')
               !         do ik = 1, nk
               !             write (1,'(4f12.6)') k_grid(ik),c(ik,5,1,1,1,1,tprint,8),c(ik,5,1,1,5,5,tprint,8),c(ik,5,3,3,3,3,tprint,8)
               !             write (2,'(4f12.6)') k_grid(ik), k(ik,5,1,1,1,1,tprint,8), k(ik,5,1,1,5,5,tprint,8), k(ik,5,3,3,3,3,tprint,8)
               !             write (3,'(4f12.6)') k_grid(ik),eV(ik,5,1,1,1,1,tprint,8),eV(ik,5,1,1,5,5,tprint,8),eV(ik,5,3,3,3,3,tprint,8)
               !             write (4,'(4f12.6)') k_grid(ik),nm(ik,5,1,1,1,1,tprint,8),nm(ik,5,1,1,5,5,tprint,8),nm(ik,5,3,3,3,3,tprint,8)
               !             write (5,'(4f12.6)') k_grid(ik),nf(ik,5,1,1,1,1,tprint,8),nf(ik,5,1,1,5,5,tprint,8),nf(ik,5,3,3,3,3,tprint,8)
               !         end do
               !         close(1)
               !         close(2)
               !         close(3)   
               !         close(4)
               !         close(5)
               ! 
               !open(1,file='csingle2.txt')
               !     open(2,file='ksingle2.txt')
               !     open(3,file='nsingle2.txt')
               !     open(4,file='vsingle2.txt')
               !     open(5,file='vsinglem2.txt')
               !     do it2 = 1, nk
               !         write(1,'(4f12.6)'), k_grid(it2), cs(1,it2,5,3,3,tprint,20),cs(2,it2,5,3,3,tprint,20),cs(2,it2,5,5,5,tprint,20)
               !         write(2,'(4f12.6)'), k_grid(it2), ks(1,it2,5,3,3,tprint,20),ks(2,it2,5,3,3,tprint,20),ks(2,it2,5,5,5,tprint,20)
               !         write(3,'(4f12.6)'), k_grid(it2), ns(1,it2,5,3,3,tprint,20),ns(2,it2,5,3,3,tprint,20),ns(2,it2,5,5,5,tprint,20)
               !         write(4,'(4f12.6)'), k_grid(it2), evs(1,it2,5,3,3,tprint,20),evs(2,it2,5,3,3,tprint,20),evs(2,it2,5,5,5,tprint,20)
               !         write(5,'(4f12.6)'), k_grid(it2), evm(1,it2,5,3,3,tprint,20),evm(2,it2,5,3,3,tprint,20),evm(2,it2,5,5,5,tprint,20)
               !     end do
               !     close(1)
               !     close(2)
               !     close(3)
               !     close(4)
               !     close(5)
               ! 
               !     open(1,file='ctest2.txt')
               !         open(2,file='ktest2.txt')
               !         open(3,file='vtest2.txt')
               !         open(4,file='nmtest2.txt')
               !         open(5,file='nftest2.txt')
               !         !open(5,file='evtest.txt')
               !         !open(4,file='ve_vu.txt')
               !         !open(4,file='tbc.txt')
               !         do ik = 1, nk
               !             write (1,'(4f12.6)') k_grid(ik),c(ik,5,1,1,1,1,tprint,20),c(ik,5,1,1,5,5,tprint,20),c(ik,5,3,3,3,3,tprint,20)
               !             write (2,'(4f12.6)') k_grid(ik), k(ik,5,1,1,1,1,tprint,20), k(ik,5,1,1,5,5,tprint,20), k(ik,5,3,3,3,3,tprint,20)
               !             write (3,'(4f12.6)') k_grid(ik),eV(ik,5,1,1,1,1,tprint,20),eV(ik,5,1,1,5,5,tprint,20),eV(ik,5,3,3,3,3,tprint,20)
               !             write (4,'(4f12.6)') k_grid(ik),nm(ik,5,1,1,1,1,tprint,20),nm(ik,5,1,1,5,5,tprint,20),nm(ik,5,3,3,3,3,tprint,20)
               !             write (5,'(4f12.6)') k_grid(ik),nf(ik,5,1,1,1,1,tprint,20),nf(ik,5,1,1,5,5,tprint,20),nf(ik,5,3,3,3,3,tprint,20)
               !         end do
               !         close(1)
               !         close(2)
               !         close(3)   
               !         close(4)
               !         close(5)
               !          
               !         
               ! STOP
                
            
            
            
            !$OMP PARALLEL PRIVATE(ik)
            !$OMP DO SCHEDULE(DYNAMIC)
            do ik=1,nsim2
                call Simulation(ik)
            end do
            !$OMP END DO    
            !$OMP END PARALLEL

            call Statistics
!STOP
        end do

        epsilon4=ratio-ratiodum

        Print *,'epsilon4 is',epsilon4

        ratio=ratio-0.2d0*(ratio-ratiodum)

        w=(1d0-alpha)*ratio**alpha
        r=alpha*ratio**(alpha-1d0)-delta

    end do

    open(1, file='singledist.txt')

    write (1, *) fpartner,mpartner

    close(1)

    open(1, file='abilityprob.txt')

    write (1, *) ability_prob

    close(1)

    open(1,file='cpathm.txt')

    do it2=1,T
        dum2=0.0
        do it4=1,16
            do it3=1,10000
                dum2=dum2+Sim1m(it4,it3,it2,2)
            end do
        end do
        dum2=dum2/160000d0
        write (1,'(F8.3,F8.3)') it2*1d0, dum2
    end do

    do it2=1,Tret
        dum2=0.0
        do it4=1,16
            do it3=1,10000
                dum2=dum2+SimR1m(it4,it3,it2,2)
            end do
        end do
        dum2=dum2/160000d0
        write (1,'(F8.3,F8.3)') it2*1d0, dum2
    end do

    close(1)

    open(1,file='cpathf.txt')

    do it2=1,T
        dum2=0.0
        do it4=1,16
            do it3=1,10000
                dum2=dum2+Sim1f(it4,it3,it2,2)
            end do
        end do
        dum2=dum2/160000d0
        write (1,'(F8.3,F8.3)') it2*1d0, dum2
    end do

    do it2=1,Tret
        dum2=0.0
        do it4=1,16
            do it3=1,10000
                dum2=dum2+SimR1f(it4,it3,it2,2)
            end do
        end do
        dum2=dum2/160000d0
        write (1,'(F8.3,F8.3)') it2*1d0, dum2
    end do

    open(1,file='kpathm.txt')

    do it2=1,T
        dum2=0.0
        do it4=1,16
            do it3=1,10000
                dum2=dum2+Sim1m(it4,it3,it2,1)
            end do
        end do
        dum2=dum2/160000d0
        write (1,'(F8.3,F8.3)') it2*1d0, dum2
    end do

    do it2=1,Tret
        dum2=0.0
        do it4=1,16
            do it3=1,10000
                dum2=dum2+SimR1m(it4,it3,it2,1)
            end do
        end do
        dum2=dum2/160000d0
        write (1,'(F8.3,F8.3)') it2*1d0, dum2
    end do

    close(1)

    open(1,file='kpathf.txt')

    do it2=1,T
        dum2=0.0
        do it4=1,16
            do it3=1,10000
                dum2=dum2+Sim1f(it4,it3,it2,1)
            end do
        end do
        dum2=dum2/160000d0
        write (1,'(F8.3,F8.3)') it2*1d0, dum2
    end do

    do it2=1,Tret
        dum2=0.0
        do it4=1,16
            do it3=1,10000
                dum2=dum2+SimR1f(it4,it3,it2,1)
            end do
        end do
        dum2=dum2/160000d0
        write (1,'(F8.3,F8.3)') it2*1d0, dum2
    end do

    close(1)

    open(1,file='npathm.txt')

    do it2=1,T
        dum2=0.0
        do it4=1,16
            do it3=1,10000
                dum2=dum2+Sim1m(it4,it3,it2,4)
            end do
        end do
        dum2=dum2/160000d0
        write (1,'(F8.3,F8.3)') it2*1d0, dum2
    end do

    close(1)

    open(1,file='npathf.txt')

    do it2=1,T
        dum2=0.0
        do it4=1,16
            do it3=1,10000
                dum2=dum2+Sim1f(it4,it3,it2,4)
            end do
        end do
        dum2=dum2/160000d0
        write (1,'(F8.3,F8.3)') it2*1d0, dum2
    end do

    close(1)

    open(1,file='lfppathsingle.txt')

    do it2=1,T
        dum2=0.0d0
        dum3=0.0d0
        do it4=1,16
            do it3=1,10000
                if(Sim1f(it4,it3,it2,10)<0.5d0) then
                    dum3=dum3+1d0
                    if(Sim1f(it4,it3,it2,4)>0.001d0) then
                        dum2=dum2+1d0
                    end if
                end if
            end do
        end do
        dum2=dum2/dum3
        write (1,'(F8.3,F8.3)') it2*1d0, dum2
    end do

    !Single female labor force participation by age after 65


do it4=1,Tret

dum2=0d0
dum3=0d0
    
do it2=1,nsim2
do it=1,nsim
    if(Sim1f(it2,it,T,10)<0.5) then
        if(SimR1f(it2,it,it4,5)>1d-3) then
        dum2=dum2+(1d0)*WeightRet(it4)
    end if
    dum3=dum3+1d0*WeightRet(it4)
    end if
end do
end do



dum2=dum2/dum3

write (1,'(F8.3,F8.3)') (64+it4)*1d0, dum2

end do
    
    close(1)

    open(1,file='lfppathmarried.txt')

    do it2=1,T
        dum2=0.0d0
        dum3=0.0d0
        do it4=1,16
            do it3=1,10000
                if(Sim1f(it4,it3,it2,10)>0.5d0) then
                    dum3=dum3+1d0
                    if(Sim1f(it4,it3,it2,4)>0.001d0) then
                        dum2=dum2+1d0
                    end if
                end if
            end do
        end do
        dum2=dum2/dum3
        write (1,'(F8.3,F8.3)') it2*1d0, dum2
    end do
    
    !Married female labor force participation by age after 65

do it4=1,Tret

dum2=0d0
dum3=0d0    
    
do it2=1,nsim2
do it=1,nsim
    if(Sim1f(it2,it,T,10)>0.5) then
        if(SimR1f(it2,it,it4,5)>1d-3) then
        dum2=dum2+(1d0)*WeightRet(it4)
    end if
    dum3=dum3+1d0*WeightRet(it4)
    end if
end do
end do

dum2=dum2/dum3

write (1,'(F8.3,F8.3)') (64+it4)*1d0, dum2

end do

    close(1)

contains

    subroutine Initialize()

        !USE ANORDF_INT
        implicit none

        allocate(v(nk,nexp,nexp,na,nu,na,nu,T,nfc,nfc))
        allocate(ev(nk,nexp,nexp,na,nu,na,nu,T,nfc,nfc))
        allocate(c(nk,nexp,nexp,na,nu,na,nu,T,nfc,nfc))
        allocate(k(nk,nexp,nexp,na,nu,na,nu,T,nfc,nfc))
        allocate(nm(nk,nexp,nexp,na,nu,na,nu,T,nfc,nfc))
        allocate(nf(nk,nexp,nexp,na,nu,na,nu,T,nfc,nfc))
        allocate(ev_spln_coefs(nk,nexp,nexp,na,nu,na,nu,T,nfc,nfc))
        allocate(v_spln_coefs(nk,nexp,nexp,na,nu,na,nu,T,nfc,nfc))
        allocate(vdum(nk,nexp,nexp,na,nu,na,nu,T,nfc,nfc))
        allocate(cdum(nk,nexp,nexp,na,nu,na,nu,T,nfc,nfc))
        allocate(gkdum(nk,nexp,nexp,na,nu,na,nu,T,nfc,nfc))
        allocate(nmdum(nk,nexp,nexp,na,nu,na,nu,T,nfc,nfc))
        allocate(nfdum(nk,nexp,nexp,na,nu,na,nu,T,nfc,nfc))

        allocate(vs(2,nk,nexp,na,nu,T,nfc))
        allocate(evs(2,nk,nexp,na,nu,T,nfc))
        allocate(evm(2,nk,nexp,na,nu,T,nfc))
        allocate(cs(2,nk,nexp,na,nu,T,nfc))
        allocate(ks(2,nk,nexp,na,nu,T,nfc))
        allocate(ns(2,nk,nexp,na,nu,T,nfc))
        allocate(evs_spln_coefs(2,nk,nexp,na,nu,T,nfc))
        allocate(vs_spln_coefs(2,nk,nexp,na,nu,T,nfc))
        allocate(evm_spln_coefs(2,nk,nexp,na,nu,T,nfc))

        allocate(Sim1m(nsim2,nsim,46,11))
        allocate(Sim1f(nsim2,nsim,46,11))
        allocate(exp1m(nsim2,nsim,46,6))
        allocate(exp1f(nsim2,nsim,46,6))
        allocate(exp2m(nsim2,nsim,46,6))
        allocate(exp2f(nsim2,nsim,46,6))
        allocate(SimR1m(nsim2,nsim,37,11))
        allocate(SimR1f(nsim2,nsim,37,11))
        allocate(expR1m(nsim2,nsim,37,4))
        allocate(expR1f(nsim2,nsim,37,4))
        allocate(Random3m(nsim2,nsim,45+36))
        allocate(Random3f(nsim2,nsim,45+36))
        allocate(marstatm(nsim2,nsim,45))
        allocate(marstatf(nsim2,nsim,45))
        allocate(marstatm_init(nsim2,nsim))
        allocate(marstatf_init(nsim2,nsim))
        allocate(partshock(nsim2,nsim,1))
        allocate(partshock2(nsim2,nsim,1))
        allocate(Random1m(nsim2,nsim))
        allocate(Random1f(nsim2,nsim))
        allocate(Random2m(nsim2,nsim))
        allocate(Random2f(nsim2,nsim))

        allocate(a(2,na))
        allocate(Prob_a(2,na))
        allocate(fc(2,nfc))
        allocate(Prob_fc(2,nfc))
        allocate(fcm(2,nfc))
        allocate(Prob_fcm(2,nfc))
        allocate(u(2,nu))
        allocate(Prob_u(2,nu))
        allocate(trans_u(2,nu,nu))
        allocate(trans_a(2,na,na))
        allocate(trans_fc(2,nfc,nfc))
        allocate(trans_fcm(2,nfc,nfc))
        allocate(OmegaRet(Tret))
        allocate(OmegaRet2(Tret))
        allocate(OmegaActive(T))
        allocate(Probm(T))
        allocate(Probd(T))
        allocate(WeightRet(Tret))
        allocate(WeightActive(T))

        allocate(fpartner(nk,nexp,na,nu,T,nfc))
        allocate(mpartner(nk,nexp,na,nu,T,nfc))
        allocate(fpartnerdum(nk,nexp,na,nu,T,nfc))
        allocate(mpartnerdum(nk,nexp,na,nu,T,nfc))
        allocate(fpartnerdum2(nk,nexp,na,nu,T,nfc))
        allocate(mpartnerdum2(nk,nexp,na,nu,T,nfc))
        allocate(ability_prob(na,na))
        allocate(laborm(nc,nw,nw))
        allocate(laborf(nc,nw,nw))
        allocate(labormwork(nc,nw))
        allocate(laborfwork(nc,nw))
        allocate(laborsinglem(nc,nw))
        allocate(laborsinglef(nc,nw))

        allocate(c_grid(nc))
        allocate(wage_grid(nw))
        allocate(k_grid(nk))
        allocate(exp_grid(nexp,T+Tret))
        allocate(K_KNOT(nk+KORDER))
        allocate(EXP_KNOT(nexp+EXPORDER,T))
        allocate(ev_spln_coefs_ret(4,nk,Tret))
        ! for testing only
        allocate(p_ev_spln_coefs_ret(4,nk,Tret))

        allocate(evs_spln_coefs_ret(4,nk,Tret))
        ! for testing only
        allocate(p_evs_spln_coefs_ret(4,nk,Tret))

        allocate(c_ret(2,2,nk,nexp,nexp,na,nu,na,nu,Tret,nfc,nfc))
        allocate(v_ret(2,2,nk,nexp,nexp,na,nu,na,nu,Tret,nfc,nfc))
        allocate(ev_ret(2,2,nk,nexp,nexp,na,nu,na,nu,Tret,nfc,nfc))
        allocate(k_ret(2,2,nk,nexp,nexp,na,nu,na,nu,Tret,nfc,nfc))
        allocate(nm_ret(2,2,nk,nexp,nexp,na,nu,na,nu,Tret,nfc,nfc))
        allocate(nf_ret(2,2,nk,nexp,nexp,na,nu,na,nu,Tret,nfc,nfc))
        allocate(vs_ret(2,2,nk,nexp,na,nu,Tret,nfc))
        allocate(evs_ret(2,2,nk,nexp,na,nu,Tret,nfc))
        allocate(cs_ret(2,2,nk,nexp,na,nu,Tret,nfc))
        allocate(ks_ret(2,2,nk,nexp,na,nu,Tret,nfc))
        allocate(ns_ret(2,2,nk,nexp,na,nu,Tret,nfc))
        allocate(break(nk))

        open(1, file='singledist.txt')
        
        read (1, *) fpartner,mpartner
        
        close(1)
        
        open(1, file='abilityprob.txt')
        
        read (1, *) ability_prob
        
        close(1)


        !Print *,ability_prob(5,:)
        !STOP
        !ability_prob=1d0/5d0
        !
        !fpartner=0d0
        !mpartner=0d0
        !
        !do it=1,T
        !    do it2=1,nexp
        !        do ik=1,8
        !            mpartner(ik,it2,:,:,it,:)=1d0/(8*it*na*nu*nfc)
        !        end do
        !    end do
        !end do
        !
        !do it=1,T
        !    do it2=1,nexp
        !        do ik=1,8
        !            fpartner(ik,it2,:,:,it,:)=1d0/(8*it*na*nu*nfc)
        !        end do
        !    end do
        !end do

        trans_u = 0d0
        prob_u=0d0
        trans_a = 0d0
        prob_a=0d0
        gamma(1,:) = (/ 0.0605927d0, -0.0010648d0, 0.0000093d0 /)
        gamma(2,:) = (/ 0.0784408d0, -0.0025596d0, 0.0000256d0 /)
        gamma0=0.12d0

        theta(:) = (/ 0.93124354*tax_level_scale, 0.15002363*tax_prog_scale /)
        thetas(:) = (/ 0.81773322*tax_level_scale, 0.11060017*tax_prog_scale /)
        AE = 0.888039805730065d0
        !Unemp_benefit=0.201795*AE
        Unemp_benefit=0d0
        r=alpha*ratio**(alpha-1d0)-delta
        w=(1d0-alpha)*ratio**alpha
        call MakeGrid(nk,k_grid,0d0,100d0,3d0)
        call MakeGrid(nc,c_grid,0.01d0,100d0,3d0)

        exp_grid=0d0
        do it2=2,T+Tret
            call MakeGrid(nexp,exp_grid(:,it2),0d0,1d0*(it2-1),2d0)
        end do

        exp_grid(:,1)=exp_grid(:,2)
        EXP_KNOT(:,1)=EXP_KNOT(:,2)

        !Print *,exp_grid(:,2)

        !STOP

        !Filling in US divorce and marriage probabilities

        probd(1)=0.11908938
        probd(2)=0.09947980
        probd(3)=0.08335590
        probd(4)=0.07024181
        probd(5)=0.05970394
        probd(6)=0.05134912
        probd(7)=0.04482276
        probd(8)=0.03980704
        probd(9)=0.03601902
        probd(10)=0.03320882
        probd(11)=0.03115778
        probd(12)=0.02967661
        probd(13)=0.02860353
        probd(14)=0.02780249
        probd(15)=0.02716124
        probd(16)=0.02658955
        probd(17)=0.02601734
        probd(18)=0.02539286
        probd(19)=0.02468080
        probd(20)=0.02386051
        probd(21)=0.02292411
        probd(22)=0.02187467
        probd(23)=0.02072434
        probd(24)=0.01949255
        probd(25)=0.01820413
        probd(26)=0.01688749
        probd(27)=0.01557275
        probd(28)=0.01428994
        probd(29)=0.01306712
        probd(30)=0.01192853
        probd(31)=0.01089280
        probd(32)=0.00997105
        probd(33)=0.00916507
        probd(34)=0.00846549
        probd(35)=0.00784992
        probd(36)=0.00728109
        probd(37)=0.00670507
        probd(38)=0.00604934
        probd(39)=0.00522102
        probd(40)=0.00410500
        probd(41)=0.00256208
        probd(42)=0.00042715
        probd(43)=0.00000000
        probd(44)=0.00000000
        probd(45)=0.00000000

        probm(1)=0.08301004
        probm(2)=0.09522406
        probm(3)=0.10423554
        probm(4)=0.11041682
        probm(5)=0.11411954
        probm(6)=0.11567485
        probm(7)=0.11539366
        probm(8)=0.11356688
        probm(9)=0.11046570
        probm(10)=0.10634176
        probm(11)=0.10142749
        probm(12)=0.09593630
        probm(13)=0.09006280
        probm(14)=0.08398313
        probm(15)=0.07785512
        probm(16)=0.07181858
        probm(17)=0.06599554
        probm(18)=0.06049049
        probm(19)=0.05539063
        probm(20)=0.05076610
        probm(21)=0.04667024
        probm(22)=0.04313985
        probm(23)=0.04019541
        probm(24)=0.03784131
        probm(25)=0.03606616
        probm(26)=0.03484296
        probm(27)=0.03412940
        probm(28)=0.03386807
        probm(29)=0.03398674
        probm(30)=0.03439857
        probm(31)=0.03500238
        probm(32)=0.03568288
        probm(33)=0.03631093
        probm(34)=0.03674376
        probm(35)=0.03682526
        probm(36)=0.03638617
        probm(37)=0.03524438
        probm(38)=0.03320513
        probm(39)=0.03006129
        probm(40)=0.02559356
        probm(41)=0.01957079
        probm(42)=0.01175015
        probm(43)=0.00187740
        probm(44)=0.00000000
        probm(45)=0.00000000

        OmegaActive=1d0

        OmegaRet(1)=1d0-0.014319d0
        OmegaRet(2)=1d0-0.015540d0
        OmegaRet(3)=1d0-0.016920d0
        OmegaRet(4)=1d0-0.018448d0
        OmegaRet(5)=1d0-0.020170d0
        OmegaRet(6)=1d0-0.022022d0
        OmegaRet(7)=1d0-0.023973d0
        OmegaRet(8)=1d0-0.026203d0
        OmegaRet(9)=1d0-0.028771d0
        OmegaRet(10)=1d0-0.031629d0
        OmegaRet(11)=1d0-0.034611d0
        OmegaRet(12)=1d0-0.037710d0
        OmegaRet(13)=1d0-0.041264d0
        OmegaRet(14)=1d0-0.045405d0
        OmegaRet(15)=1d0-0.050128d0
        OmegaRet(16)=1d0-0.055339d0
        OmegaRet(17)=1d0-0.061005d0
        OmegaRet(18)=1d0-0.067396d0
        OmegaRet(19)=1d0-0.074476d0
        OmegaRet(20)=1d0-0.082272d0
        OmegaRet(21)=1d0-0.091816d0
        OmegaRet(22)=1d0-0.101898d0
        OmegaRet(23)=1d0-0.112870d0
        OmegaRet(24)=1d0-0.124763d0
        OmegaRet(25)=1d0-0.137597d0
        OmegaRet(26)=1d0-0.151383d0
        OmegaRet(27)=1d0-0.166117d0
        OmegaRet(28)=1d0-0.181778d0
        OmegaRet(29)=1d0-0.198331d0
        OmegaRet(30)=1d0-0.215721d0
        OmegaRet(31)=1d0-0.233874d0
        OmegaRet(32)=1d0-0.252699d0
        OmegaRet(33)=1d0-0.272086d0
        OmegaRet(34)=1d0-0.291912d0
        OmegaRet(35)=1d0-0.312040d0
        OmegaRet(36)=1d0-1d0

        OmegaRet2(1)=1d0
        do it2=1,Tret-1
            OmegaRet2(it2+1)=OmegaRet(it2)
        end do
        
        

        call tauchen_hans(sigma_am,rho_am,na,a(1,:),trans_a(1,:,:),prob_a(1,:))
        call tauchen_hans(sigma_um,rho_um,nu,u(1,:),trans_u(1,:,:),prob_u(1,:))
        
        call tauchen_hans(sigma_af,rho_af,na,a(2,:),trans_a(2,:,:),prob_a(2,:))
        call tauchen_hans(sigma_uf,rho_uf,nu,u(2,:),trans_u(2,:,:),prob_u(2,:))

        call tauchen_hans(sigma_fcm,rho_fcm,nfc,fc(1,:),trans_fc(1,:,:),prob_fc(1,:))
        call tauchen_hans(sigma_fcs,rho_fcs,nfc,fc(2,:),trans_fc(2,:,:),prob_fc(2,:))
        
        call tauchen_hans(sigma_fcmm,rho_fcmm,nfc,fcm(1,:),trans_fcm(1,:,:),prob_fcm(1,:))
        call tauchen_hans(sigma_fcsm,rho_fcsm,nfc,fcm(2,:),trans_fcm(2,:,:),prob_fcm(2,:))

        

        fc(1,:)=exp(fc(1,:)+mu_fcm)
        fc(2,:)=exp(fc(2,:)+mu_fcs)
        fcm(1,:)=exp(fcm(1,:)+mu_fcmm)
        fcm(2,:)=exp(fcm(2,:)+mu_fcsm)

        !open(1, file='transprob.txt')
        !    
        !    write (1, *) trans_u, Prob_u, trans_a, prob_a, trans_fc, prob_fc
        !    
        !close(1)
        !STOP
        !Print *,prob_u(1,:)
        !
        !Print *,trans_fc(1,1,1)
        !Print *,trans_fc(1,1,2)
        !Print *,trans_fc(1,1,3)
        !Print *,trans_fc(1,1,4)
        !Print *,trans_fc(1,1,5)
        !Print *,trans_fc(1,1,6)
        !Print *,trans_fc(1,1,7)
        !Print *,trans_fc(1,1,8)
        !Print *,trans_fc(1,1,9)
        !Print *,trans_fc(1,1,10)
        !Print *,trans_fc(1,1,11)
        !Print *,trans_fc(1,1,12)
        !Print *,trans_fc(1,1,13)
        !Print *,trans_fc(1,1,14)
        !Print *,trans_fc(1,1,15)

        !Print *,trans_u(1,4,:)
        !Print *,trans_u(1,5,:)
        !
        !dum=0D0
        !do it=1,15
        !    dum=dum+trans_fc(1,it,2)
        !end do

        !dum=0D0
        !do it=1,15
        !    dum=dum+prob_fc(2,it)
        !end do

        !dum=probfam(1,15)+probfam(2,15)+probfam(3,15)
        !
        !Print *,dum

        !Print *,fc(1,:)
        !
        !Print *,fc(2,:)
        !
        !Print *,prob_fc(1,:)
        !
        !Print *,prob_fc(2,:)
        
        !Print *,fcm(1,:)
        !
        !Print *,fcm(2,:)
        !
        !Print *,prob_fcm(1,:)
        !
        !Print *,prob_fcm(2,:)
        !
        !STOP

        call random_number(random1m)

        call random_number(random2m)

        call random_number(random3m)

        call random_number(random1f)

        call random_number(random2f)

        call random_number(random3f)

        call random_number(marstatm)

        call random_number(marstatf)

        call random_number(marstatm_init)

        call random_number(marstatf_init)

        call random_number(partshock)
        
        call random_number(partshock2)


        !    open(1, file='random1.txt')
        !    
        !        read (1, *) random1
        !        
        !    close(1)
        !    
        !    open(2, file='random2.txt')
        !    
        !        read (2, *) random2
        !        
        !    close(2)
        !    
        !    open(3, file='random3.txt')
        !    
        !        read (3, *) random3
        !        
        !    close(3)


    end subroutine Initialize  


end program Laffer
