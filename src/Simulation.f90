subroutine simulation(ik)
    !This subroutine simulates the lifecycle for nsim households
    use Model_Parameters
    use PolicyFunctions
    use Utilities

    implicit none

    integer, INTENT(IN) :: ik
    integer :: i,iam,ium,iaf,iuf,um,uf,it2,it3,it4,it5,j,count2,ifc,ifcm,iu2,ir,irm
    real(8) :: ix,ixm,dum2,dum3,dum4,dum5,dum6,&
        mixm(nsim,2), mixf(nsim,2),mixdum(nsim,2),mixdum2(nsim,2),pnt1(2),pnt2(3),exp_grid_dum(nexp),INTERP2D(nk,nexp),INTERP3D(nk,nexp,nexp)
    
    
    dum2=0.0
    !Print *,'Simulation',dum2


    !These vectors will hold various statistics

    Sim1m(ik,:,:,:)=0d0
    Sim1m(ik,:,:,:)=0d0

    !Assigning the asset level for 20-year olds

    Sim1m(ik,:,1,1)=0.440852*AE
    Sim1f(ik,:,1,1)=0.440852*AE

    !This vector will hold the level of experience for women, abilit2y and the idiosyncratic shock

    exp2m(ik,:,:,1)=0d0
    exp2f(ik,:,:,1)=0d0

    exp1m(ik,:,:,4)=0
    exp1f(ik,:,:,4)=0

    !Initialzing the idiosyncratic shocks to wages
    
    do it2=1,nsim
        if(Random1m(ik,it2)<Prob_u(1,1)) then
            exp1m(ik,it2,1,3)=1
        elseif((Prob_u(1,1)<Random1m(ik,it2)).AND.(Random1m(ik,it2)<(Prob_u(1,1)+Prob_u(1,2)))) then
            exp1m(ik,it2,1,3)=2
        elseif(((Prob_u(1,1)+Prob_u(1,2))<Random1m(ik,it2)).AND.(Random1m(ik,it2)<(Prob_u(1,1)+Prob_u(1,2)+Prob_u(1,3)))) then
            exp1m(ik,it2,1,3)=3
        elseif(((Prob_u(1,1)+Prob_u(1,2)+Prob_u(1,3))<Random1m(ik,it2)).AND.(Random1m(ik,it2)<(Prob_u(1,1)+Prob_u(1,2)+Prob_u(1,3)+Prob_u(1,4)))) then
            exp1m(ik,it2,1,3)=4
        else
            exp1m(ik,it2,1,3)=5
        end if
    end do
    
    do it2=1,nsim
        if(Random1f(ik,it2)<Prob_u(2,1)) then
            exp1f(ik,it2,1,3)=1
        elseif((Prob_u(2,1)<Random1f(ik,it2)).AND.(Random1f(ik,it2)<(Prob_u(2,1)+Prob_u(2,2)))) then
            exp1f(ik,it2,1,3)=2
        elseif(((Prob_u(2,1)+Prob_u(2,2))<Random1f(ik,it2)).AND.(Random1f(ik,it2)<(Prob_u(2,1)+Prob_u(2,2)+Prob_u(2,3)))) then
            exp1f(ik,it2,1,3)=3
        elseif(((Prob_u(2,1)+Prob_u(2,2)+Prob_u(2,3))<Random1f(ik,it2)).AND.(Random1f(ik,it2)<(Prob_u(2,1)+Prob_u(2,2)+Prob_u(2,3)+Prob_u(2,4)))) then
            exp1f(ik,it2,1,3)=4
        else
            exp1f(ik,it2,1,3)=5
        end if
    end do
    
    !Initialzing the distribution of abilit2ies
    
    do it2=1,nsim
        if(Random2m(ik,it2)<Prob_a(1,1)) then
            exp1m(ik,it2,:,2)=1
        elseif((Prob_a(1,1)<Random2m(ik,it2)).AND.(Random2m(ik,it2)<(Prob_a(1,1)+Prob_a(1,2)))) then
            exp1m(ik,it2,:,2)=2
        elseif(((Prob_a(1,1)+Prob_a(1,2))<Random2m(ik,it2)).AND.(Random2m(ik,it2)<(Prob_a(1,1)+Prob_a(1,2)+Prob_a(1,3)))) then
            exp1m(ik,it2,:,2)=3
        elseif(((Prob_a(1,1)+Prob_a(1,2)+Prob_a(1,3))<Random2m(ik,it2)).AND.(Random2m(ik,it2)<(Prob_a(1,1)+Prob_a(1,2)+Prob_a(1,3)+Prob_a(1,4)))) then
            exp1m(ik,it2,:,2)=4
        else
            exp1m(ik,it2,:,2)=5
        end if
    end do
    
    do it2=1,nsim
        if(Random2f(ik,it2)<Prob_a(2,1)) then
            exp1f(ik,it2,:,2)=1
        elseif((Prob_a(2,1)<Random2f(ik,it2)).AND.(Random2f(ik,it2)<(Prob_a(2,1)+Prob_a(2,2)))) then
            exp1f(ik,it2,:,2)=2
        elseif(((Prob_a(2,1)+Prob_a(2,2))<Random2f(ik,it2)).AND.(Random2f(ik,it2)<(Prob_a(2,1)+Prob_a(2,2)+Prob_a(2,3)))) then
            exp1f(ik,it2,:,2)=3
        elseif(((Prob_a(2,1)+Prob_a(2,2)+Prob_a(2,3))<Random2f(ik,it2)).AND.(Random2f(ik,it2)<(Prob_a(2,1)+Prob_a(2,2)+Prob_a(2,3)+Prob_a(2,4)))) then
            exp1f(ik,it2,:,2)=4
        else
            exp1f(ik,it2,:,2)=5
        end if
    end do
    
!Initializing marital status and partner number
    
!Initializing marital status and partner number

mixm=0d0
mixf=0d0
count2=0

do it2=1,nsim

!Single men

mixm(it2,2)=it2*1d0

exp1m(ik,it2,1,4)=0
Sim1m(ik,it2,1,10)=0d0

if(marstatm_init(ik,it2)<0.1259) then
    count2=count2+1
    mixm(it2,1)=marstatm_init(ik,it2)+match*A(1,exp1m(ik,it2,1,2))+1000d0
end if

end do

it3=0
do it2=1,nsim
    
!Single women
it3=it3+1
exp1f(ik,it2,1,4)=0
Sim1f(ik,it2,1,10)=0d0
mixf(it2,2)=it2*1d0

if(it3<count2+1) then
    mixf(it2,1)=marstatf_init(ik,it2)+match*A(2,exp1f(ik,it2,1,2))+1000d0
end if

end do

!Sorting single men and women by marriage shock

mixdum=0d0
mixdum2=0d0


do it2=1,nsim
    it3=0
    dum2=0d0
    do
        it3=it3+1
        if(mixm(it2,1)>mixdum(it3,1)) then
            do it4=it3,nsim
                mixdum2(it4,1)=mixdum(it4,1)
                mixdum2(it4,2)=mixdum(it4,2)
            end do
            mixdum(it3,1)=mixm(it2,1)
            mixdum(it3,2)=mixm(it2,2)
            
            do it4=(it3+1),nsim
                mixdum(it4,1)=mixdum2(it4-1,1)
                mixdum(it4,2)=mixdum2(it4-1,2)
            end do
            dum2=1d0
        end if
        if(it3==nsim) then
            mixdum(it3,1)=mixm(it2,1)
            mixdum(it3,2)=mixm(it2,2)
            dum2=1d0
        end if
        
        if(dum2>0.5d0) exit
        
    end do
    
end do

mixm=mixdum

!Print *,'Minval of mixm is',minval(mixm(:,2))
!Print *,'Minloc of mixm is',minloc(mixm(:,2))

mixdum=0d0
mixdum2=0d0

do it2=1,nsim
    it3=0
    dum2=0d0
    do
        it3=it3+1
        if(mixf(it2,1)>mixdum(it3,1)) then
            do it4=it3,nsim
                mixdum2(it4,1)=mixdum(it4,1)
                mixdum2(it4,2)=mixdum(it4,2)
            end do
            mixdum(it3,1)=mixf(it2,1)
            mixdum(it3,2)=mixf(it2,2)
            
            do it4=(it3+1),nsim
                mixdum(it4,1)=mixdum2(it4-1,1)
                mixdum(it4,2)=mixdum2(it4-1,2)
            end do
            dum2=1d0
        end if
        if(it3==nsim) then
            mixdum(it3,1)=mixf(it2,1)
            mixdum(it3,2)=mixf(it2,2)
            dum2=1d0
        end if
        if(dum2>0.5d0) exit
    end do
end do

mixf=mixdum    
    
!Arranging marriages

!count2=nsim
!it3=int(0.1259*nsim)

do it2=1,count2
        it4=int(mixm(it2,2)+0.01)
        it5=int(mixf(it2,2)+0.01)
        exp1m(ik,it4,1,4)=it5
        exp1f(ik,it5,1,4)=it4
        Sim1m(ik,it4,1,10)=1d0
        Sim1f(ik,it5,1,10)=1d0
        Sim1m(ik,it4,1,1)=Sim1m(ik,it4,1,1)+Sim1f(ik,it5,1,1)
        Sim1f(ik,it5,1,1)=Sim1m(ik,it4,1,1)
end do
    
    !Initialzing the participation costs
    
    do it2=1,nsim
        iu2=1
        dum5=Prob_fc(1,iu2)
        do while((dum5<Partshock(ik,it2,1)).AND.(iu2<nfc))
            iu2=iu2+1
            dum5=dum5+Prob_fc(1,iu2)
        end do
        exp1f(ik,it2,:,6)=iu2
    end do

   do it2=1,nsim
        iu2=1
        dum5=Prob_fcm(1,iu2)
        do while((dum5<Partshock2(ik,it2,1)).AND.(iu2<nfc))
            iu2=iu2+1
            dum5=dum5+Prob_fcm(1,iu2)
        end do
        exp1m(ik,it2,:,6)=iu2
    end do
    
    !do it2=1,100
    !    Print *,exp1(ik,it2,1,3)
    !end do
    !
    !Print *,prob_u(1,:)
    !
    !Print *,prob_u(2,:)
    !
    !STOP
    
    !Starting to simulate age 25-64, using piecewise linear interpolation of the policy functions
    
    do i=1,T
    
        exp_grid_dum=exp_grid(:,i)
        
        do it2=1,nsim
    
            !print *, it2, nsim
    
            if(Sim1m(ik,it2,i,10)>0.5d0) then
    
                it3=exp1m(ik,it2,i,4)
                dum2=Sim1m(ik,it2,i,1)
                iam=exp1m(ik,it2,i,2)
                ium=exp1m(ik,it2,i,3)
                ixm=exp2m(ik,it2,i,1)
                ifcm=exp1m(ik,it2,i,6)
                iaf=exp1f(ik,it3,i,2)
                iuf=exp1f(ik,it3,i,3)
                ix=exp2f(ik,it3,i,1)
                ifc=exp1f(ik,it3,i,6)
                pnt2 = (/dum2, ix, ixm/)
    
                !Next period's capital
                INTERP3D=k(:,:,:,iam,ium,iaf,iuf,i,ifc,ifcm)
                Sim1m(ik,it2,i+1,1) = trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                Sim1f(ik,it3,i+1,1) = Sim1m(ik,it2,i+1,1)
    
                !This period's consumption
                INTERP3D=c(:,:,:,iam,ium,iaf,iuf,i,ifc,ifcm)
                dum4 = trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                Sim1m(ik,it2,i,2) = dum4
                Sim1f(ik,it3,i,2) = dum4
    
                ! Male wage, work hours and earnings
                Sim1m(ik,it2,i,3) = wage(1,a(1,iam),dble(ixm),u(1,ium))/(1d0+t_employer)
                INTERP3D=nm(:,:,:,iam,ium,iaf,iuf,i,ifc,ifcm)
                dum5=trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                dum5=max(dum5,0d0)
                Sim1m(ik,it2,i,4) = dum5
                Sim1m(ik,it2,i,5) = dum5*wage(1,a(1,iam),dble(ixm),u(1,ium))/(1d0+t_employer)
    
                ! Female wage, work hours and earnings
    
                Sim1f(ik,it3,i,3) = wage(2,a(2,iaf),dble(ix),u(2,iuf))/(1d0+t_employer)
                INTERP3D=nf(:,:,:,iam,ium,iaf,iuf,i,ifc,ifcm)
                dum6=trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                dum6=max(dum6,0d0)
                Sim1f(ik,it3,i,4) = dum6
                Sim1f(ik,it3,i,5) = dum6*wage(2,a(2,iaf),dble(ix),u(2,iuf))/(1d0+t_employer)
    
                !Household income and taxes
    
                dum3=dum5*wage(1,a(1,iam),dble(ixm),u(1,ium))/(1d0+t_employer)+dum6*wage(2,a(2,iaf),dble(ix),u(2,iuf))/(1d0+t_employer)
                Sim1m(ik,it2,i,6) = dum3
                Sim1f(ik,it3,i,6) = dum3
                if(dum3>0d0) then
                    Sim1m(ik,it2,i,7)= tax_labor(dum3)*dum3
                    Sim1f(ik,it3,i,7)= tax_labor(dum3)*dum3
                else
                    Sim1m(ik,it2,i,7)= 0d0
                    Sim1f(ik,it3,i,7)= 0d0
                end if
    
                Sim1m(ik,it2,i,8)= dum4*tc
                Sim1f(ik,it3,i,8)= dum4*tc
                Sim1m(ik,it2,i,9)=dum3*t_employee+t_employer*dum3
                Sim1f(ik,it3,i,9)=dum3*t_employee+t_employer*dum3
    
                if(Sim1f(ik,it3,i,4)>1d-3) then
                    exp2f(ik,it3,i+1,1)=exp2f(ik,it3,i,1)+1d0
                else
                    exp2f(ik,it3,i+1,1)=exp2f(ik,it3,i,1)*(1d0-deltaexp)
                end if
                
                if(Sim1m(ik,it2,i,4)>1d-3) then
                    exp2m(ik,it2,i+1,1)=exp2m(ik,it2,i,1)+1d0
                else
                    exp2m(ik,it2,i+1,1)=exp2m(ik,it2,i,1)*(1d0-deltaexp)
                end if
    
                !Social Welfare
                INTERP3D=V(:,:,:,iam,ium,iaf,iuf,i,ifc,ifcm)
                dum3 = trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                Sim1m(ik,it2,i,11) = dum3
                Sim1f(ik,it3,i,11) = dum3
                
            else
    
                dum2=Sim1m(ik,it2,i,1)
                iam=exp1m(ik,it2,i,2)
                ium=exp1m(ik,it2,i,3)
                ixm=exp2m(ik,it2,i,1)
                ifcm=exp1m(ik,it2,i,6)
                j=1
                pnt1 = (/dum2, ixm/)
    
                !Next period's capital
                INTERP2D=ks(j,:,:,iam,ium,i,ifcm)
                Sim1m(ik,it2,i+1,1)=bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
    
                !This period's consumption
                INTERP2D=cs(j,:,:,iam,ium,i,ifcm)
                dum4 = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                Sim1m(ik,it2,i,2) = dum4
    
                ! Male wage, work hours and earnings
                Sim1m(ik,it2,i,3) = wage(1,a(1,iam),dble(ixm),u(1,ium))/(1d0+t_employer)
                INTERP2D=ns(j,:,:,iam,ium,i,ifcm)
                dum5=bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                dum5=max(dum5,0d0)
                Sim1m(ik,it2,i,4) = dum5
                Sim1m(ik,it2,i,5) = dum5*wage(1,a(1,iam),dble(ixm),u(1,ium))/(1d0+t_employer)
    
                !Household income and taxes
    
                dum3=dum5*wage(1,a(1,iam),dble(ixm),u(1,ium))/(1d0+t_employer)
                Sim1m(ik,it2,i,6) = dum3
                if(dum3>0d0) then
                    Sim1m(ik,it2,i,7)= tax_labors(dum3)*dum3
                else
                    Sim1m(ik,it2,i,7)=0d0
                end if
    
                Sim1m(ik,it2,i,8)= dum4*tc
                Sim1m(ik,it2,i,9)=dum3*t_employee+t_employer*dum3
                
                if(Sim1m(ik,it2,i,4)>1d-3) then
                    exp2m(ik,it2,i+1,1)=exp2m(ik,it2,i,1)+1d0
                else
                    exp2m(ik,it2,i+1,1)=exp2m(ik,it2,i,1)*(1d0-deltaexp)
                end if
                
                !Social welfare
                INTERP2D=Vs(j,:,:,iam,ium,i,ifcm)
                dum3 = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                Sim1m(ik,it2,i,11) = dum3
    
            end if
    
            !Single women
    
            if(Sim1f(ik,it2,i,10)<0.5d0) then
    
                dum2=Sim1f(ik,it2,i,1)
                iam=exp1f(ik,it2,i,2)
                ium=exp1f(ik,it2,i,3)
                ix=exp2f(ik,it2,i,1)
                ifc=exp1f(ik,it2,i,6)
                j=2
                pnt1 = (/dum2, ix/)
    
                !Next period's capital
                INTERP2D=ks(j,:,:,iam,ium,i,ifc)
                Sim1f(ik,it2,i+1,1)=bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
    
                !This period's consumption
                INTERP2D=cs(j,:,:,iam,ium,i,ifc)
                dum4=bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                Sim1f(ik,it2,i,2) = dum4    
    
                ! Female wage, work hours and earnings
    
                Sim1f(ik,it2,i,3) = wage(2,a(2,iam),dble(ix),u(2,ium))/(1d0+t_employer)
                INTERP2D=ns(j,:,:,iam,ium,i,ifc)
                dum5=bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                dum5=max(dum5,0d0)
                Sim1f(ik,it2,i,4) = dum5
                Sim1f(ik,it2,i,5) = dum5*wage(2,a(2,iam),dble(ix),u(2,ium))/(1d0+t_employer)
    
                !Household income and taxes
    
                dum3=dum5*wage(2,a(2,iam),dble(ix),u(2,ium))/(1d0+t_employer)
                Sim1f(ik,it2,i,6) = dum3
                if(dum3>0d0) then
                    Sim1f(ik,it2,i,7)= tax_labors(dum3)*dum3
                else
                    Sim1f(ik,it2,i,7)=0d0
                end if
    
                Sim1f(ik,it2,i,8)= dum4*tc
                Sim1f(ik,it2,i,9)=dum3*t_employee+t_employer*dum3
    
                if(Sim1f(ik,it2,i,4)>1d-3) then
                    exp2f(ik,it2,i+1,1)=exp2f(ik,it2,i,1)+1d0
                else
                    exp2f(ik,it2,i+1,1)=exp2f(ik,it2,i,1)*(1d0-deltaexp)
                end if
                
                !Social Welfare
                INTERP2D=Vs(j,:,:,iam,ium,i,ifc)
                dum3=bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                Sim1f(ik,it2,i,11) = dum3
    
            end if
    
        end do
    
    
        !print *, i
    
        !Uppdating the idiosyncratic wage shock
    
        do it2=1,nsim
            um=exp1m(ik,it2,i,3)
            if(Random3m(ik,it2,i)<trans_u(1,um,1)) then
                exp1m(ik,it2,i+1,3)=1
            elseif((trans_u(1,um,1)<Random3m(ik,it2,i)).AND.(Random3m(ik,it2,i)<(trans_u(1,um,1)+trans_u(1,um,2)))) then
                exp1m(ik,it2,i+1,3)=2
            elseif(((trans_u(1,um,1)+trans_u(1,um,2))<Random3m(ik,it2,i)).AND.(Random3m(ik,it2,i)<(trans_u(1,um,1)+trans_u(1,um,2)+trans_u(1,um,3)))) then
                exp1m(ik,it2,i+1,3)=3
            elseif(((trans_u(1,um,1)+trans_u(1,um,2)+trans_u(1,um,3))<Random3m(ik,it2,i)).AND.(Random3m(ik,it2,i)<(trans_u(1,um,1)+trans_u(1,um,2)+trans_u(1,um,3)+trans_u(1,um,4)))) then
                exp1m(ik,it2,i+1,3)=4
            else
                exp1m(ik,it2,i+1,3)=5
            end if
        end do
    
        do it2=1,nsim
            uf=exp1f(ik,it2,i,3)
            if(Random3f(ik,it2,i)<trans_u(2,uf,1)) then
                exp1f(ik,it2,i+1,3)=1
            elseif((trans_u(2,uf,1)<Random3f(ik,it2,i)).AND.(Random3f(ik,it2,i)<(trans_u(2,uf,1)+trans_u(2,uf,2)))) then
                exp1f(ik,it2,i+1,3)=2
            elseif(((trans_u(2,uf,1)+trans_u(2,uf,2))<Random3f(ik,it2,i)).AND.(Random3f(ik,it2,i)<(trans_u(2,uf,1)+trans_u(2,uf,2)+trans_u(2,uf,3)))) then
                exp1f(ik,it2,i+1,3)=3
            elseif(((trans_u(2,uf,1)+trans_u(2,uf,2)+trans_u(2,uf,3))<Random3f(ik,it2,i)).AND.(Random3f(ik,it2,i)<(trans_u(2,uf,1)+trans_u(2,uf,2)+trans_u(2,uf,3)+trans_u(2,uf,4)))) then
                exp1f(ik,it2,i+1,3)=4
            else
                exp1f(ik,it2,i+1,3)=5
            end if
        end do
    
        count2=0
        mixm=0d0
        mixf=0d0
    
        !dum5=0d0
        !do it2=1,nsim
        !    if((Sim1m(ik,it2,i,10)<0.5).AND.(Sim1m(ik,it2,i+1,1)>dum5)) then
        !        dum5=Sim1m(ik,it2,i+1,1)
        !    end if
        !end do
        !
        !dum5=0.0
        !
        !do it2=1,nsim
        !    if((Sim1f(ik,it2,i,10)<0.5).AND.(Sim1f(ik,it2,i+1,1)>dum5)) then
        !        dum5=Sim1f(ik,it2,i+1,1)
        !    end if
        !end do
    
        !dum5=maxval(Sim1m(ik,:,i+1,1))
        !
        !dum6=maxval(Sim1m(ik,:,i,2))
    
        do it2=1,nsim
    
            !Married staying married
    
            if((Sim1m(ik,it2,i,10)>0.5d0).AND.(marstatm(ik,it2,i)>Probd(i))) then
                it3=exp1m(ik,it2,i,4)
                exp1m(ik,it2,i+1,4)=it3
                exp1f(ik,it3,i+1,4)=it2
                Sim1m(ik,it2,i+1,10)=1d0
                Sim1f(ik,it3,i+1,10)=1d0
            end if
    
            !Married getting divorced
    
            if((Sim1m(ik,it2,i,10)>0.5d0).AND.(marstatm(ik,it2,i)<Probd(i))) then
                it3=exp1m(ik,it2,i,4)
                exp1m(ik,it2,i+1,4)=0
                exp1f(ik,it3,i+1,4)=0
                Sim1m(ik,it2,i+1,10)=0d0
                Sim1f(ik,it3,i+1,10)=0d0
                Sim1m(ik,it2,i+1,1)=Sim1m(ik,it2,i+1,1)/2d0
                Sim1f(ik,it3,i+1,1)=Sim1f(ik,it2,i+1,1)/2d0
    
                exp1f(ik,it3,i+1,5)=0
    
            end if
    
!Single men

if(Sim1m(ik,it2,i,10)<0.5d0) then
    exp1m(ik,it2,i+1,4)=0
    Sim1m(ik,it2,i+1,10)=0d0
    if(marstatm(ik,it2,i)<Probm(i)) then
        mixm(it2,2)=it2*1d0
        mixm(it2,1)=marstatm(ik,it2,i)+match*A(1,exp1m(ik,it2,i,2))+1000d0
        count2=count2+1
    end if
end if
    
        end do
    
it3=0
do it2=1,nsim
    !Single women
   
    if(Sim1f(ik,it2,i,10)<0.5d0) then
        it3=it3+1
        exp1f(ik,it2,i+1,4)=0
        Sim1f(ik,it2,i+1,10)=0d0
        if(it3<count2+1) then
            mixf(it2,2)=it2*1d0
            mixf(it2,1)=marstatf(ik,it2,i)+match*A(2,exp1f(ik,it2,i,2))+1000d0
        end if
    end if
end do
    
!Sorting single men and women by marriage shock Mn

mixdum=0d0
mixdum2=0d0

do it2=1,nsim
    it3=0
    dum2=0d0
    do
        it3=it3+1
        if(mixm(it2,1)>mixdum(it3,1)) then
            do it4=it3,nsim
                mixdum2(it4,1)=mixdum(it4,1)
                mixdum2(it4,2)=mixdum(it4,2)
            end do
            mixdum(it3,1)=mixm(it2,1)
            mixdum(it3,2)=mixm(it2,2)
            
            do it4=(it3+1),nsim
                mixdum(it4,1)=mixdum2(it4-1,1)
                mixdum(it4,2)=mixdum2(it4-1,2)
            end do
            dum2=1d0
        end if
        if(it3==nsim) then
            mixdum(it3,1)=mixm(it2,1)
            mixdum(it3,2)=mixm(it2,2)
            dum2=1d0
        end if
        
        if(dum2>0.5d0) exit
        
    end do
    
end do


mixm=mixdum


mixdum=0d0
mixdum2=0d0

do it2=1,nsim
    it3=0
    dum2=0d0
    do
        it3=it3+1
        if(mixf(it2,1)>mixdum(it3,1)) then
            do it4=it3,nsim
                mixdum2(it4,1)=mixdum(it4,1)
                mixdum2(it4,2)=mixdum(it4,2)
            end do
            mixdum(it3,1)=mixf(it2,1)
            mixdum(it3,2)=mixf(it2,2)
            
            do it4=(it3+1),nsim
                mixdum(it4,1)=mixdum2(it4-1,1)
                mixdum(it4,2)=mixdum2(it4-1,2)
            end do
            dum2=1d0
        end if
        if(it3==nsim) then
            mixdum(it3,1)=mixf(it2,1)
            mixdum(it3,2)=mixf(it2,2)
            dum2=1d0
        end if
        if(dum2>0.5d0) exit
    end do
end do

mixf=mixdum



!Arranging marriages (this means updating next period's marital status, assigning each spouse a partner number (who is he/she married too) and summng the assets of the couple)

!it3=int(count2*Probm(i))

do it2=1,count2
    !if(it2<it3+1) then
        it4=int(mixm(it2,2)+0.01)
        it5=int(mixf(it2,2)+0.01)
        exp1m(ik,it4,i+1,4)=it5
        exp1f(ik,it5,i+1,4)=it4
        Sim1m(ik,it4,i+1,10)=1d0
        Sim1f(ik,it5,i+1,10)=1d0
        Sim1m(ik,it4,i+1,1)=Sim1m(ik,it4,i+1,1)+Sim1f(ik,it5,i+1,1)
        Sim1f(ik,it5,i+1,1)=Sim1m(ik,it4,i+1,1)
    !else
        !it4=int(mixm(it2,2)+0.01)
        !it5=int(mixf(it2,2)+0.01)
        !exp1m(ik,it4,i+1,4)=0
        !exp1f(ik,it5,i+1,4)=0
        !Sim1m(ik,it4,i+1,10)=0d0
        !Sim1f(ik,it5,i+1,10)=0d0
    !end if
end do

dum5=maxval(Sim1m(ik,:,i+1,1))



end do    

    
    !Beginning the simulation for retired households
    
    SimR1m(ik,:,1,1)=Sim1m(ik,:,T+1,1)
    SimR1m(ik,:,1,5)=exp2m(ik,:,T+1,1)
    expR1m(ik,:,1,2)=exp1m(ik,:,T+1,3)
    expR1m(ik,:,1,3)=2
    SimR1f(ik,:,1,1)=Sim1f(ik,:,T+1,1)
    SimR1f(ik,:,1,5)=exp2f(ik,:,T+1,1)
    expR1f(ik,:,1,2)=exp1f(ik,:,T+1,3)
    expR1f(ik,:,1,3)=2
    
    do i=1,Tret
        expR1m(ik,:,i,4)=exp1m(ik,:,T+1,6)
        expR1m(ik,:,i,1)=exp1m(ik,:,T+1,2)
        expR1f(ik,:,i,4)=exp1f(ik,:,T+1,6)
        expR1f(ik,:,i,1)=exp1f(ik,:,T+1,2)
    end do
    
    
    do i=1,36
        
        exp_grid_dum=exp_grid(:,T+i)
    
        do it2=1,nsim
    
            if(Sim1m(ik,it2,T,10)>0.5) then
    
                it3=exp1m(ik,it2,T,4)
                dum2=SimR1m(ik,it2,i,1)
                ixm=SimR1m(ik,it2,i,5)
                iam=expR1m(ik,it2,i,1)
                ium=expR1m(ik,it2,i,2)
                irm=expR1m(ik,it2,i,3)
                ifcm=expR1m(ik,it2,i,4)
                ix=SimR1f(ik,it3,i,5)
                iaf=expR1f(ik,it3,i,1)
                iuf=expR1f(ik,it3,i,2)
                ir=expR1f(ik,it3,i,3)
                ifc=expR1f(ik,it3,i,4)
                pnt2 = (/dum2, ix, ixm/)
    
                !Next period's capital
                INTERP3D=k_ret(ir,irm,:,:,:,iam,ium,iaf,iuf,i,ifc,ifcm)
                SimR1m(ik,it2,i+1,1) = trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                SimR1f(ik,it3,i+1,1) = SimR1m(ik,it2,i+1,1)
                
                !Next periods's consumption
                INTERP3D=c_ret(ir,irm,:,:,:,iam,ium,iaf,iuf,i,ifc,ifcm)
                dum4 = trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                SimR1m(ik,it2,i,2) = dum4
                SimR1f(ik,it3,i,2) = dum4
                SimR1m(ik,it2,i,3)= dum4*tc
                SimR1f(ik,it3,i,3)= dum4*tc
                
                ! Male wage, work hours and earnings
                
                INTERP3D=nm_ret(ir,irm,:,:,:,iam,ium,iaf,iuf,i,ifc,ifcm)
                SimR1m(ik,it2,i,4) = wage(1,a(1,iam),dble(ixm),u(1,ium))/(1d0+t_employer)
                dum5=trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                dum5=max(dum5,0d0)
                SimR1m(ik,it2,i,5) = dum5
                SimR1m(ik,it2,i,6) = dum5*wage(1,a(1,iam),dble(ixm),u(1,ium))/(1d0+t_employer)
                
                ! Female wage, work hours and earnings
                
                INTERP3D=nf_ret(ir,irm,:,:,:,iam,ium,iaf,iuf,i,ifc,ifcm)
                SimR1f(ik,it3,i,4) = wage(2,a(2,iaf),dble(ix),u(2,iuf))/(1d0+t_employer)
                dum6=trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                dum6=max(dum6,0d0)
                SimR1f(ik,it3,i,5) = dum6
                SimR1f(ik,it3,i,6) = dum6*wage(2,a(2,iaf),dble(ix),u(2,iuf))/(1d0+t_employer)
    
                !Household income and taxes
    
                dum3=dum5*wage(1,a(1,iam),dble(ixm),u(1,ium))/(1d0+t_employer)+dum6*wage(2,a(2,iaf),dble(ix),u(2,iuf))/(1d0+t_employer)
                SimR1f(ik,it3,i,7) = dum3
                SimR1m(ik,it2,i,7) = dum3
                if(dum3>0d0) then
                    SimR1f(ik,it3,i,8) = tax_labor(dum3)*dum3
                    SimR1m(ik,it2,i,8) = tax_labor(dum3)*dum3
                else
                    SimR1f(ik,it3,i,8) = 0d0
                    SimR1m(ik,it2,i,8) = 0d0
                end if
                
                SimR1f(ik,it3,i,9)=dum3*t_employee+t_employer*dum3
                SimR1m(ik,it2,i,9)=dum3*t_employee+t_employer*dum3
                
                if(SimR1m(ik,it2,i,5)>1d-3) then
                    SimR1m(ik,it2,i+1,5)=SimR1m(ik,it2,i,5)+1d0
                    expR1m(ik,it2,i+1,3)=2
                else
                    SimR1m(ik,it2,i+1,5)=SimR1m(ik,it2,i,5)*(1d0-deltaexp)
                    expR1m(ik,it2,i+1,3)=1
                end if
                
                if(SimR1f(ik,it3,i,5)>1d-3) then
                    SimR1f(ik,it3,i+1,5)=SimR1f(ik,it3,i,5)+1d0
                    expR1f(ik,it3,i+1,3)=2
                else
                    SimR1f(ik,it3,i+1,5)=SimR1f(ik,it3,i,5)*(1d0-deltaexp)
                    expR1f(ik,it3,i+1,3)=1
                end if
    
                !Social welfare
                INTERP3D=v_ret(ir,irm,:,:,:,iam,ium,iaf,iuf,i,ifc,ifcm)
                dum3 = trilin_interp(k_grid, exp_grid_dum, exp_grid_dum, INTERP3D, nk, nexp, nexp, pnt2)
                SimR1m(ik,it2,i,11) = dum3
                SimR1f(ik,it3,i,11) = dum3
                
            else
    
                dum2=SimR1m(ik,it2,i,1)
                ixm=SimR1m(ik,it2,1,5)
                iam=expR1m(ik,it2,i,1)
                ium=expR1m(ik,it2,i,2)
                irm=expR1m(ik,it2,i,3)
                ifcm=expR1m(ik,it2,i,4)
                j=1
                pnt1 = (/dum2, ixm/)
                
                !Next period's capital
                INTERP2D=ks_ret(irm,j,:,:,iam,ium,i,ifcm)
                SimR1m(ik,it2,i+1,1) = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                
                !Next periods's consumption
                INTERP2D=cs_ret(irm,j,:,:,iam,ium,i,ifcm)
                dum4 = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                SimR1m(ik,it2,i,2) = dum4
                SimR1m(ik,it2,i,3)= dum4*tc
                
                ! Female wage, work hours and earnings
                
                INTERP2D=ns_ret(irm,j,:,:,iam,ium,i,ifcm)
                SimR1m(ik,it2,i,4) = wage(1,a(1,iam),dble(ixm),u(1,ium))/(1d0+t_employer)
                dum5=bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                dum5=max(dum5,0d0)
                SimR1m(ik,it2,i,5) = dum5
                SimR1m(ik,it2,i,6) = dum5*wage(1,a(1,iam),dble(ixm),u(1,ium))/(1d0+t_employer)
    
                !Household income and taxes
    
                dum3=dum5*wage(1,a(1,iam),dble(ixm),u(1,ium))/(1d0+t_employer)
                SimR1m(ik,it2,i,7) = dum3
                if(dum3>0d0) then
                    SimR1m(ik,it2,i,8) = tax_labors(dum3)*dum3
                else
                    SimR1m(ik,it2,i,8) = 0d0
                end if
                
                SimR1m(ik,it2,i,9)=dum3*t_employee+t_employer*dum3
                
                if(SimR1m(ik,it2,i,5)>1d-3) then
                    SimR1m(ik,it2,i+1,5)=SimR1m(ik,it2,i,5)+1d0
                    expR1m(ik,it2,i+1,3)=2
                else
                    SimR1m(ik,it2,i+1,5)=SimR1m(ik,it2,i,5)*(1d0-deltaexp)
                    expR1m(ik,it2,i+1,3)=1
                end if
    
                !Social welfare
                INTERP2D=vs_ret(ir,j,:,:,iaf,iuf,i,ifc)
                dum3 = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                SimR1m(ik,it2,i,11) = dum3
                
            end if
    
            
    
            if(Sim1f(ik,it2,T,10)<0.5) then
    
                dum2=SimR1f(ik,it2,i,1)
                ix=SimR1f(ik,it2,1,5)
                iaf=expR1f(ik,it2,i,1)
                iuf=expR1f(ik,it2,i,2)
                ir=expR1f(ik,it2,i,3)
                ifc=expR1f(ik,it2,i,4)
                j=2
                pnt1 = (/dum2, ix/)
                
                !Next period's capital
                INTERP2D=ks_ret(ir,j,:,:,iaf,iuf,i,ifc)
                SimR1f(ik,it2,i+1,1) = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                
                !Next periods's consumption
                INTERP2D=cs_ret(ir,j,:,:,iaf,iuf,i,ifc)
                dum4 = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                SimR1f(ik,it2,i,2) = dum4
                SimR1f(ik,it2,i,3)= dum4*tc
                
                ! Female wage, work hours and earnings
                
                INTERP2D=ns_ret(ir,j,:,:,iaf,iuf,i,ifc)
                SimR1f(ik,it2,i,4) = wage(2,a(2,iaf),dble(ix),u(2,iuf))/(1d0+t_employer)
                dum5=bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                dum5=max(dum5,0d0)
                SimR1f(ik,it2,i,5) = dum5
                SimR1f(ik,it2,i,6) = dum5*wage(2,a(2,iaf),dble(ix),u(2,iuf))/(1d0+t_employer)
    
                !Household income and taxes
    
                dum3=dum5*wage(2,a(2,iaf),dble(ix),u(2,iuf))/(1d0+t_employer)
                SimR1f(ik,it2,i,7) = dum3
                if(dum3>0d0) then
                    SimR1f(ik,it2,i,8) = tax_labors(dum3)*dum3
                else
                    SimR1f(ik,it2,i,8) = 0d0
                end if
                
                SimR1f(ik,it2,i,9)=dum3*t_employee+t_employer*dum3
                
                if(SimR1f(ik,it2,i,5)>1d-3) then
                    SimR1f(ik,it2,i+1,5)=SimR1f(ik,it2,i,5)+1d0
                    expR1f(ik,it2,i+1,3)=2
                else
                    SimR1f(ik,it2,i+1,5)=SimR1f(ik,it2,i,5)*(1d0-deltaexp)
                    expR1f(ik,it2,i+1,3)=1
                end if
                
                
                !Social welfare
                INTERP2D=vs_ret(ir,j,:,:,iaf,iuf,i,ifc)
                dum3 = bilin_interp(k_grid, exp_grid_dum, INTERP2D, nk, nexp, pnt1)
                SimR1f(ik,it2,i,11) = dum3
    
            end if
    
            
            
        end do
    
        do it2=1,nsim
            um=expR1m(ik,it2,i,2)
            if(Random3m(ik,it2,T+i)<trans_u(2,um,1)) then
                expR1m(ik,it2,i+1,2)=1
            elseif((trans_u(2,um,1)<Random3m(ik,it2,T+i)).AND.(Random3m(ik,it2,T+i)<(trans_u(2,um,1)+trans_u(2,um,2)))) then
                expR1m(ik,it2,i+1,2)=2
            elseif(((trans_u(2,um,1)+trans_u(2,um,2))<Random3m(ik,it2,T+i)).AND.(Random3m(ik,it2,T+i)<(trans_u(2,um,1)+trans_u(2,um,2)+trans_u(2,um,3)))) then
                expR1m(ik,it2,i+1,2)=3
            elseif(((trans_u(2,um,1)+trans_u(2,um,2)+trans_u(2,um,3))<Random3m(ik,it2,T+i)).AND.(Random3m(ik,it2,T+i)<(trans_u(2,um,1)+trans_u(2,um,2)+trans_u(2,um,3)+trans_u(2,um,4)))) then
                expR1m(ik,it2,i+1,2)=4
            else
                expR1m(ik,it2,i+1,2)=5
            end if
        end do
    
        do it2=1,nsim
            uf=expR1f(ik,it2,i,2)
            if(Random3f(ik,it2,T+i)<trans_u(2,uf,1)) then
                expR1f(ik,it2,i+1,2)=1
            elseif((trans_u(2,uf,1)<Random3f(ik,it2,T+i)).AND.(Random3f(ik,it2,T+i)<(trans_u(2,uf,1)+trans_u(2,uf,2)))) then
                expR1f(ik,it2,i+1,2)=2
            elseif(((trans_u(2,uf,1)+trans_u(2,uf,2))<Random3f(ik,it2,T+i)).AND.(Random3f(ik,it2,T+i)<(trans_u(2,uf,1)+trans_u(2,uf,2)+trans_u(2,uf,3)))) then
                expR1f(ik,it2,i+1,2)=3
            elseif(((trans_u(2,uf,1)+trans_u(2,uf,2)+trans_u(2,uf,3))<Random3f(ik,it2,T+i)).AND.(Random3f(ik,it2,T+i)<(trans_u(2,uf,1)+trans_u(2,uf,2)+trans_u(2,uf,3)+trans_u(2,uf,4)))) then
                expR1f(ik,it2,i+1,2)=4
            else
                expR1f(ik,it2,i+1,2)=5
            end if
        end do
        
    end do

end subroutine Simulation
