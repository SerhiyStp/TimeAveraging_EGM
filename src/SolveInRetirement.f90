subroutine SolveInRetirement(ik)
! Both men and women are retired (absorbing state)
    use Model_Parameters
    use PolicyFunctions
    use Utilities
    use MyInterpolation
    
    implicit none
    
    integer, INTENT(IN) :: ik
    real(8) :: c2, v2
    real(8) :: d1, d2,P1,P2,P3,P4,v3
    integer :: it2
    real(8) :: test1, test2
    real(8) :: vnext_p2, vnext_p3

    !Married    
    if(it==1) then
        ! Solve the very last period problem:           
        k_ret(1,1,ik,:,:,:,:,:,:,Tret,:,:) = 0d0
        c_ret(1,1,ik,:,:,:,:,:,:,Tret,:,:) = ((k_grid(ik) + Gamma_redistr)*(((1d0 + r)/OmegaRet2(Tret+1-it))*(1d0-tk)) + Psi_pension+lumpsum)/(1d0+tc)
        nf_ret(1,1,ik,:,:,:,:,:,:,Tret,:,:)=0d0
        nm_ret(1,1,ik,:,:,:,:,:,:,Tret,:,:)=0d0
        v_ret(1,1,ik,:,:,:,:,:,:,Tret,:,:) = Uc(c_ret(1,1,ik,1,1,1,1,1,1,Tret,1,1)) 
    else
        ! Solve the Tret-1 to 1st period of retirement:
        !Finding optimal capital by golden search
        P1=k_grid(1)
        P4=k_grid(ik)
        do
            P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
            P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)
            vnext_p2 = S_QDVAL(P2,k_grid,v_ret(1,1,:,1,1,1,1,1,1,Tret-it+2,1,1),nk)
            vnext_p3 = S_QDVAL(P3,k_grid,v_ret(1,1,:,1,1,1,1,1,1,Tret-it+2,1,1),nk)
            V2=Uc(((k_grid(ik) + Gamma_redistr)*(((1d0+r)/OmegaRet2(Tret+1-it))*(1d0-tk))+Psi_pension+lumpsum-P2*(1d0+mu))/(1.0+tc))                       
            V2=V2+beta*OmegaRet(Tret-it+1)*vnext_p2
            if(k_grid(ik)*(((1d0+r)/OmegaRet2(Tret+1-it))*(1d0-tk))+Psi_pension+lumpsum-P2*(1d0+mu)<0.01d0) then
                V2=-999999999.0
            end if

            V3=Uc(((k_grid(ik) + Gamma_redistr)*(((1d0+r)/OmegaRet2(Tret+1-it))*(1d0-tk))+Psi_pension+lumpsum-P3*(1d0+mu))/(1.0+tc))                       
            V3=V3+beta*OmegaRet(Tret-it+1)*vnext_p3
            if(k_grid(ik)*(((1d0+r)/OmegaRet2(Tret+1-it))*(1d0-tk))+Psi_pension+lumpsum-P3*(1d0+mu)<0.01d0) then
                V3=-999999999.0
            end if       
            if (V2 < V3) then
                P1=P2
            else
                P4=P3
            end if
            if((P4-P1)<1d-8) exit
        end do
        v_ret(1,1,ik,:,:,:,:,:,:,Tret-it+1,:,:)=V2
        k_ret(1,1,ik,:,:,:,:,:,:,Tret-it+1,:,:)=P2
        nf_ret(1,1,ik,:,:,:,:,:,:,Tret-it+1,:,:)=0d0
        nm_ret(1,1,ik,:,:,:,:,:,:,Tret-it+1,:,:)=0d0
        c_ret(1,1,ik,:,:,:,:,:,:,Tret-it+1,:,:)=((k_grid(ik) + Gamma_redistr)*(((1d0+r)/OmegaRet2(Tret+1-it))*(1d0-tk))+Psi_pension+lumpsum-P2*(1d0+mu))/(1d0+tc)
    end if

    !Single
    if(it==1) then
        ! Solve the very last period problem:           
        ks_ret(1,:,ik,:,:,:,Tret,:) = 0d0
        cs_ret(1,:,ik,:,:,:,Tret,:) = ((k_grid(ik) + Gamma_redistr*0.5d0)*(((1d0 + r)/OmegaRet2(Tret+1-it))*(1d0-tk)) + Psi_pension*0.5d0+lumpsum*0.5d0)/(1d0+tc)
        ns_ret(1,:,ik,:,:,:,Tret,:) = 0d0
        vs_ret(1,:,ik,:,:,:,Tret,:) = Uc(cs_ret(1,1,ik,1,1,1,Tret,1)) 
    else
        ! Solve the Tret-1 to 1st period of retirement:
        !Finding optimal capital by golden search
        P1=k_grid(1)
        P4=k_grid(ik)
        do
            P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
            P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)
            !vnext_p2 = D_QDVAL(P2,k_grid,vs_ret(1,1,:,1,1,1,Tret-it+2,1))
            vnext_p2 = S_QDVAL(P2,k_grid,vs_ret(1,1,:,1,1,1,Tret-it+2,1),nk)
            !vnext_p3 = D_QDVAL(P3,k_grid,vs_ret(1,1,:,1,1,1,Tret-it+2,1))
            vnext_p3 = S_QDVAL(P3,k_grid,vs_ret(1,1,:,1,1,1,Tret-it+2,1),nk)


            V2=Uc(((k_grid(ik) + Gamma_redistr*0.5d0)*(((1d0+r)/OmegaRet2(Tret+1-it))*(1d0-tk))+Psi_pension*0.5d0+lumpsum*0.5d0-P2*(1d0+mu))/(1.0+tc))                       
            !V2=V2+beta*OmegaRet(Tret-it+1)*LinInterp(P2,k_grid,v_ret(:,Tret-it+2),nk)
            !V2=V2+beta*OmegaRet(Tret-it+1)*D_CSVAL(P2,BREAK,evs_spln_coefs_ret(:,:,Tret-it+2))
            V2=V2+beta*OmegaRet(Tret-it+1)*vnext_p2
            if((k_grid(ik) + Gamma_redistr*0.5d0)*(((1d0+r)/OmegaRet2(Tret+1-it))*(1d0-tk))+Psi_pension*0.5d0+lumpsum*0.5d0-P2*(1d0+mu)<0.01d0) then
                V2=-999999999.0
            end if
            V3=Uc(((k_grid(ik) + Gamma_redistr*0.5d0)*(((1d0+r)/OmegaRet2(Tret+1-it))*(1d0-tk))+Psi_pension*0.5d0+lumpsum*0.5d0-P3*(1d0+mu))/(1.0+tc))                       
            !V3=V3+beta*OmegaRet(Tret-it+1)*LinInterp(P3,k_grid,v_ret(:,Tret-it+2),nk)
            !V3=V3+beta*OmegaRet(Tret-it+1)*D_CSVAL(P3,BREAK,evs_spln_coefs_ret(:,:,Tret-it+2))
            V3=V3+beta*OmegaRet(Tret-it+1)*vnext_p3
            if((k_grid(ik) + Gamma_redistr*0.5d0)*(((1d0+r)/OmegaRet2(Tret+1-it))*(1d0-tk))+Psi_pension*0.5d0+lumpsum*0.5d0-P3*(1d0+mu)<0.01d0) then
                V3=-999999999.0
            end if       
            if (V2 < V3) then
                P1=P2
            else
                P4=P3
            end if
            if((P4-P1)<1d-8) exit
        end do
        vs_ret(1,:,ik,:,:,:,Tret-it+1,:)=V2
        ks_ret(1,:,ik,:,:,:,Tret-it+1,:)=P2
        ns_ret(1,:,ik,:,:,:,Tret-it+1,:)=0d0
        cs_ret(1,:,ik,:,:,:,Tret-it+1,:)=((k_grid(ik) + Gamma_redistr*0.5d0)*(((1d0+r)/OmegaRet2(Tret+1-it))*(1d0-tk))+Psi_pension*0.5d0+lumpsum*0.5d0-P2*(1d0+mu))/(1d0+tc)
    end if

end subroutine SolveInRetirement
    
subroutine SolveInRetirement2(counter)
    ! Woman or Man can choose to work
    use Model_Parameters
    use PolicyFunctions
    use glob0
    use Utilities

    implicit none

    integer, INTENT(IN) :: counter
    integer :: ix,ixd,iam,ium,iaf,iuf,iu2,iu3,tprint,ikd,j,ifc,ik
    real(8) :: ce,cu,ke,ku,nem,nef,num,nuf,ve,vu
    real(8) :: ces,cus,kes,kus,nes,nus,ves,vus
    integer :: NEQ=0, IERSVR=0, IPACT=0, ISACT=0
    real(8) :: c2, MU2, d1, d2, vp(nu),dum3,dum4,dum5,dum6,y
    real(8) :: ACC=0.0001d0,ERREL=0.0001d0
    real(8) :: P1,P2,P3,P4,V2,V3,dum2,pnt2(3),pnt1(2)
    real(8) :: vnext, exp_grid_dum(nexp), INTERP2D(nk,nexp)

    dum3=((counter*1d0)/(nu*na*1d0))-0.00001d0
    ik=int(dum3)+1
    dum3=(((counter-(ik-1)*nu*na)*1d0)/(na*1d0))-0.00001d0
    ium=int(dum3)+1
    iam=counter-(ik-1)*na*nu-(ium-1)*nu

    if(it>1) then
        exp_grid_dum=exp_grid(:,T+Tret+2-it)
    end if

    ce=0d0
    cu=0d0
    ke=0d0
    ku=0d0
    nem=0d0
    nef=0d0
    num=0d0
    nuf=0d0
    ve=0d0
    vu=0d0


    ces=0d0
    cus=0d0
    kes=0d0
    kus=0d0
    nes=0d0
    nus=0d0
    ves=0d0
    vus=0d0

    !Married woman is working

    do ifc=1,nfc
        do ix = 1, nexp
            wagef = wage(2,a(2,iam),exp_grid(ix,T+Tret+1-it),u(2,ium))/(1d0+t_employer)

            if(it==1) then
                ! Solve the very last period problem:           
                ke = 0d0
                ce = ((k_grid(ik) + Gamma_redistr)*(((1d0 + r)/OmegaRet2(Tret+1-it))*(1d0-tk)) + Psi_pension+lumpsum)/(1d0+tc)
                nem = 0d0
                nef = 0d0
                ve = Uc(ce) 
            else
                ! Solve the Tret-1 to 1st period of retirement:
                !Finding optimal capital by golden search
                P1=0.01d0
                P4 = min((k_grid(nk)-0.001d0)/(1d0+tc), &
                    ((k_grid(ik) + Gamma_redistr)*(1d0+(r/OmegaRet2(Tret+1-it))*(1d0-tk))+lumpsum&
                    +Psi_pension*0.5d0+wagef*(1d0-tax_labor(wagef)-tSS_employee(wagef)))/(1d0+tc))
                do
                    P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                    P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)

                    pnt1 = (/P2, wagef/)
                    dum4 = min(max(bilin_interp(c_grid, wage_grid, laborfwork, nc, nw, pnt1),0d0),1d0)
                    y=dum4*wagef
                    dum2=((k_grid(ik) + Gamma_redistr)*(1d0+(r/OmegaRet2(Tret+1-it))*(1d0-tk)) &
                        +lumpsum+Psi_pension*0.5d0+y*(1d0-tax_labor(y)-tSS_employee(y))-P2*(1d0+tc))/(1d0+mu)
                    if(dum2<0.0001d0) then
                        V2=-999999999d0
                    else
                        V2=Uc(P2)-chif*(dum4**(1d0+etaf))/(1d0+etaf)-fc(1,ifc)

                        pnt1=[dum2, exp_grid(ix,T+Tret-it)+1d0]
                        vnext = bilin_interp( &
                                    k_grid, exp_grid_dum, &
                                    ev_ret(2,1,:,:,1,iam,ium,1,1,Tret+2-it,ifc,1), &
                                    nk, nexp, pnt1 )
                        V2=V2+beta*OmegaRet(Tret-it+1)*vnext
                    end if

                    pnt1 = [P3, wagef]
                    dum4 = min(max(bilin_interp(c_grid, wage_grid, laborfwork, nc, nw, pnt1),0d0),1d0)
                    y=dum4*wagef
                    dum2 = ((k_grid(ik) + Gamma_redistr)*(1d0+(r/OmegaRet2(Tret+1-it))*(1d0-tk)) &
                        +lumpsum+Psi_pension*0.5d0+y*(1d0-tax_labor(y)-tSS_employee(y))-P3*(1d0+tc))/(1d0+mu)
                    if(dum2<0.0001d0) then
                        V3=-999999999d0
                    else
                        V3=Uc(P3)-chif*(dum4**(1d0+etaf))/(1d0+etaf)-fc(1,ifc)
                        pnt1 = [dum2, exp_grid(ix,T+Tret-it)+1d0]
                        vnext = bilin_interp(&
                                    k_grid, exp_grid_dum,ev_ret(2,1,:,:,1,iam,ium,1,1,Tret+2-it,ifc,1), &
                                    nk, nexp, pnt1)
                        V3=V3+beta*OmegaRet(Tret-it+1)*vnext
                    end if


                    if (V2 < V3) then
                        P1=P2
                    else
                        P4=P3
                    end if
                    if((P4-P1)<1d-8) exit
                end do
                Ve=V2
                ke=dum2
                nef=dum4
                nem=0d0
                ce=P2
            end if


            vu=v_ret(1,1,ik,ix,1,iam,ium,1,1,Tret-it+1,ifc,1)
            ku=k_ret(1,1,ik,ix,1,iam,ium,1,1,Tret-it+1,ifc,1)
            nuf=0d0
            num=0d0
            cu=c_ret(1,1,ik,ix,1,iam,ium,1,1,Tret-it+1,ifc,1)

            if (ve >= vu) then
                v_ret(2,1,ik,ix,:,iam,ium,:,:,Tret-it+1,ifc,:)=ve
                c_ret(2,1,ik,ix,:,iam,ium,:,:,Tret-it+1,ifc,:)=ce
                k_ret(2,1,ik,ix,:,iam,ium,:,:,Tret-it+1,ifc,:)=ke
                nf_ret(2,1,ik,ix,:,iam,ium,:,:,Tret-it+1,ifc,:)=nef
                nm_ret(2,1,ik,ix,:,iam,ium,:,:,Tret-it+1,ifc,:)=nem
            else
                v_ret(2,1,ik,ix,:,iam,ium,:,:,Tret-it+1,ifc,:)=vu
                c_ret(2,1,ik,ix,:,iam,ium,:,:,Tret-it+1,ifc,:)=cu
                k_ret(2,1,ik,ix,:,iam,ium,:,:,Tret-it+1,ifc,:)=ku
                nf_ret(2,1,ik,ix,:,iam,ium,:,:,Tret-it+1,ifc,:)=nuf
                nm_ret(2,1,ik,ix,:,iam,ium,:,:,Tret-it+1,ifc,:)=num
            end if

        end do
    end do

    ce=0d0
    cu=0d0
    ke=0d0
    ku=0d0
    nem=0d0
    nef=0d0
    num=0d0
    nuf=0d0
    ve=0d0
    vu=0d0

    !Married man is working

    do ifc=1,nfc
        do ix = 1, nexp
            wagem = wage(1,a(1,iam),exp_grid(ix,T+Tret+1-it),u(1,ium))/(1d0+t_employer)


            if(it==1) then
                ! Solve the very last period problem:
                ke = 0d0
                ce = ((k_grid(ik) + Gamma_redistr)*(((1d0 + r)/OmegaRet2(Tret+1-it))*(1d0-tk)) + Psi_pension+lumpsum)/(1d0+tc)
                nem = 0d0
                nef = 0d0
                ve = Uc(ce)
            else
                ! Solve the Tret-1 to 1st period of retirement:
                !Finding optimal capital by golden search
                P1=0.01d0
                P4=min((k_grid(nk)-0.001d0)/(1d0+tc),((k_grid(ik) + Gamma_redistr)*(1d0+(r/OmegaRet2(Tret+1-it))*(1d0-tk))+lumpsum+Psi_pension*0.5d0+wagem*(1d0-tax_labor(wagem)-tSS_employee(wagem)))/(1d0+tc))
                do
                    P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                    P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)

                    pnt1 = (/P2, wagem/)
                    dum4 = min(max(bilin_interp(c_grid, wage_grid, labormwork, nc, nw, pnt1),0d0),1d0)
                    y=dum4*wagem
                    dum2=((k_grid(ik) + Gamma_redistr)*(1d0+(r/OmegaRet2(Tret+1-it))*(1d0-tk))+lumpsum+Psi_pension*0.5d0+y*(1d0-tax_labor(y)-tSS_employee(y))-P2*(1d0+tc))/(1d0+mu)
                    if(dum2<0.0001d0) then
                        V2=-999999999d0
                    else
                        V2=Uc(P2)-chim*(dum4**(1d0+etam))/(1d0+etam)-fcm(1,ifc)

                        pnt1 = [dum2, exp_grid(ix,T+Tret-it)+1d0]
                        vnext = bilin_interp( &
                                    k_grid, exp_grid(:,T+Tret-it+2), &
                                    ev_ret(1,2,:,1,:,iam,ium,1,1,Tret+2-it,1,ifc), nk, nexp, pnt1)
                        !vnext = D_QD2VL(dum2,exp_grid(ix,T+Tret+1-it)+1d0,k_grid,exp_grid(:,T+Tret-it+2),ev_ret(1,2,:,1,:,iam,ium,1,1,Tret+2-it,1,ifc))
                        V2=V2+beta*OmegaRet(Tret-it+1)*vnext
                    end if

                    pnt1 = [P3, wagem]
                    dum4 = min(max(bilin_interp(c_grid, wage_grid, labormwork, nc, nw, pnt1),0d0),1d0)
                    y=dum4*wagem
                    dum2=((k_grid(ik) + Gamma_redistr)*(1d0+(r/OmegaRet2(Tret+1-it))*(1d0-tk))+lumpsum+Psi_pension*0.5d0+y*(1d0-tax_labor(y)-tSS_employee(y))-P3*(1d0+tc))/(1d0+mu)
                    if(dum2<0.0001d0) then
                        V3=-999999999d0
                    else
                        V3=Uc(P3)-chim*(dum4**(1d0+etam))/(1d0+etam)-fcm(1,ifc)

                        pnt1 = [dum2, exp_grid(ix,T+Tret-it)+1d0]
                        vnext = bilin_interp( &
                                    k_grid, exp_grid(:,T+Tret-it+2), &
                                    ev_ret(1,2,:,1,:,iam,ium,1,1,Tret+2-it,1,ifc), nk, nexp, pnt1)
                        !vnext = D_QD2VL(dum2,exp_grid(ix,T+Tret+1-it)+1d0,k_grid,exp_grid(:,T+Tret-it+2),ev_ret(1,2,:,1,:,iam,ium,1,1,Tret+2-it,1,ifc))
                        V3=V3+beta*OmegaRet(Tret-it+1)*vnext
                    end if

                    if (V2 < V3) then
                        P1=P2
                    else
                        P4=P3
                    end if
                    if((P4-P1)<1d-8) exit
                end do
                Ve=V2
                ke=dum2
                nef=0d0
                nem=dum4
                ce=P2
            end if

            vu=v_ret(1,1,ik,ix,1,iam,ium,1,1,Tret-it+1,ifc,1)
            ku=k_ret(1,1,ik,ix,1,iam,ium,1,1,Tret-it+1,ifc,1)
            nuf=0d0
            num=0d0
            cu=c_ret(1,1,ik,ix,1,iam,ium,1,1,Tret-it+1,ifc,1)

            if (ve >= vu) then
                v_ret(1,2,ik,:,ix,iam,ium,:,:,Tret-it+1,:,ifc)=ve
                c_ret(1,2,ik,:,ix,iam,ium,:,:,Tret-it+1,:,ifc)=ce
                k_ret(1,2,ik,:,ix,iam,ium,:,:,Tret-it+1,:,ifc)=ke
                nf_ret(1,2,ik,:,ix,iam,ium,:,:,Tret-it+1,:,ifc)=nef
                nm_ret(1,2,ik,:,ix,iam,ium,:,:,Tret-it+1,:,ifc)=nem
            else
                v_ret(1,2,ik,:,ix,iam,ium,:,:,Tret-it+1,:,ifc)=vu
                c_ret(1,2,ik,:,ix,iam,ium,:,:,Tret-it+1,:,ifc)=cu
                k_ret(1,2,ik,:,ix,iam,ium,:,:,Tret-it+1,:,ifc)=ku
                nf_ret(1,2,ik,:,ix,iam,ium,:,:,Tret-it+1,:,ifc)=nuf
                nm_ret(1,2,ik,:,ix,iam,ium,:,:,Tret-it+1,:,ifc)=num
            end if
        end do
    end do


    !Single women

    j=2

    do ifc=1,nfc
        do ix = 1, nexp
            wagef = wage(2,a(2,iam),exp_grid(ix,T+Tret+1-it),u(2,ium))/(1d0+t_employer)

            if(it==1) then
                ! Solve the very last period problem:           
                kes = 0d0
                ces = ((k_grid(ik) + Gamma_redistr*0.5d0)*(((1d0 + r)/OmegaRet2(Tret+1-it))*(1d0-tk)) + Psi_pension*0.5d0+lumpsum*0.5d0)/(1d0+tc)
                nes = 0d0
                ves = Uc(ces) 
            else
                ! Solve the Tret-1 to 1st period of retirement:
                !Finding optimal capital by golden search
                P1=0.01d0
                P4=min((k_grid(nk)-0.001d0)/(1d0+tc),((k_grid(ik) + Gamma_redistr*0.5d0)*(1d0+(r/OmegaRet2(Tret+1-it))*(1d0-tk))+lumpsum*0.5d0+wagef*(1d0-tax_labors(wagef)-tSS_employee(wagef)))/(1d0+tc))
                do
                    P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                    P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)

                    pnt1 = [P2, wagef]
                    dum4 = min(max(bilin_interp(c_grid, wage_grid, laborsinglef, nc, nw, pnt1),0d0),1d0)
                    y=dum4*wagef
                    dum2=((k_grid(ik) + Gamma_redistr*0.5d0)*(1d0+(r/OmegaRet2(Tret+1-it))*(1d0-tk))+lumpsum*0.5d0+y*(1d0-tax_labors(y)-tSS_employee(y))-P2*(1d0+tc))/(1d0+mu)
                    if(dum2<0.0001d0) then
                        V2=-999999999d0
                    else
                        V2=Uc(P2)-chifs*(dum4**(1d0+etaf))/(1d0+etaf)-fc(2,ifc)

                        pnt1 = [dum2, exp_grid(ix,T+Tret+1-it)+1d0]
                        vnext = bilin_interp( &
                                    k_grid, exp_grid_dum, &
                                    evs_ret(2,j,:,:,iam,ium,Tret+2-it,ifc), &
                                    nk, nexp, pnt1)
                        V2=V2+beta*OmegaRet(Tret-it+1)*vnext

                    end if

                    pnt1 = [P3, wagef]
                    dum4 = min(max(bilin_interp(c_grid, wage_grid, laborsinglef, nc, nw, pnt1),0d0),1d0)
                    y=dum4*wagef
                    dum2=((k_grid(ik) + Gamma_redistr*0.5d0)*(1d0+(r/OmegaRet2(Tret+1-it))*(1d0-tk))+lumpsum*0.5d0+y*(1d0-tax_labors(y)-tSS_employee(y))-P3*(1d0+tc))/(1d0+mu)
                    if(dum2<0.0001d0) then
                        V3=-999999999d0
                    else
                        V3=Uc(P3)-chifs*(dum4**(1d0+etaf))/(1d0+etaf)-fc(2,ifc)
                        pnt1 = [dum2, exp_grid(ix,T+Tret+1-it)+1d0]
                        vnext = bilin_interp( &
                                    k_grid, exp_grid_dum, &
                                    evs_ret(2,j,:,:,iam,ium,Tret+2-it,ifc), &
                                    nk, nexp, pnt1)
                        V3=V3+beta*OmegaRet(Tret-it+1)*vnext
                    end if
                    if (V2 < V3) then
                        P1=P2
                    else
                        P4=P3
                    end if
                    if((P4-P1)<1d-8) exit
                end do
                Ves=V2
                kes=dum2
                nes=dum4
                ces=P2
            end if

            vus=vs_ret(1,j,ik,ix,iam,ium,Tret-it+1,ifc)
            kus=ks_ret(1,j,ik,ix,iam,ium,Tret-it+1,ifc)
            nus=ns_ret(1,j,ik,ix,iam,ium,Tret-it+1,ifc)
            cus=cs_ret(1,j,ik,ix,iam,ium,Tret-it+1,ifc)

            if (ves >= vus) then
                vs_ret(2,j,ik,ix,iam,ium,Tret-it+1,ifc)=ves
                cs_ret(2,j,ik,ix,iam,ium,Tret-it+1,ifc)=ces
                ks_ret(2,j,ik,ix,iam,ium,Tret-it+1,ifc)=kes
                ns_ret(2,j,ik,ix,iam,ium,Tret-it+1,ifc)=nes
            else
                vs_ret(2,j,ik,ix,iam,ium,Tret-it+1,ifc)=vus
                cs_ret(2,j,ik,ix,iam,ium,Tret-it+1,ifc)=cus
                ks_ret(2,j,ik,ix,iam,ium,Tret-it+1,ifc)=kus
                ns_ret(2,j,ik,ix,iam,ium,Tret-it+1,ifc)=nus
            end if

        end do
    end do

    !Single men

    j=1

    do ifc=1,nfc
        do ix = 1, nexp
            wagem = wage(1,a(1,iam),exp_grid(ix,T+Tret+1-it),u(1,ium))/(1d0+t_employer)

            if(it==1) then
                ! Solve the very last period problem:           
                kes = 0d0
                ces = ((k_grid(ik) + Gamma_redistr*0.5d0)*(((1d0 + r)/OmegaRet2(Tret+1-it))*(1d0-tk)) + Psi_pension*0.5d0+lumpsum*0.5d0)/(1d0+tc)
                nes = 0d0
                ves = Uc(ces) 
            else
                ! Solve the Tret-1 to 1st period of retirement:
                !Finding optimal capital by golden search
                P1=0.01d0
                P4=min((k_grid(nk)-0.001d0)/(1d0+tc),((k_grid(ik) + Gamma_redistr*0.5d0)*(1d0+(r/OmegaRet2(Tret+1-it))*(1d0-tk))+lumpsum*0.5d0+wagem*(1d0-tax_labors(wagem)-tSS_employee(wagem)))/(1d0+tc))
                do
                    P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                    P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)

                    pnt1 = [P2, wagem]
                    dum4 = min(max(bilin_interp(c_grid, wage_grid, laborsinglem, nc, nw, pnt1),0d0),1d0)
                    y = dum4*wagem
                    dum2 = ((k_grid(ik) + Gamma_redistr*0.5d0)*(1d0+(r/OmegaRet2(Tret+1-it))*(1d0-tk)) &
                        +lumpsum*0.5d0+y*(1d0-tax_labors(y)-tSS_employee(y))-P2*(1d0+tc))/(1d0+mu)
                    if(dum2<0.0001d0) then
                        V2=-999999999d0
                    else
                        V2=Uc(P2)-chims*(dum4**(1d0+etam))/(1d0+etam)-fcm(2,ifc)
                        pnt1 = [dum2, exp_grid(ix,T+Tret+1-it)+1d0]
                        vnext = bilin_interp( &
                                    k_grid, exp_grid_dum, &
                                    evs_ret(2,j,:,:,iam,ium,Tret+2-it,ifc), &
                                    nk, nexp, pnt1)
                        V2=V2+beta*OmegaRet(Tret-it+1)*vnext
                    end if

                    pnt1 = [P3, wagem]
                    dum4 = min(max(bilin_interp(c_grid, wage_grid, laborsinglem, nc, nw, pnt1),0d0),1d0)
                    y=dum4*wagem
                    dum2 = ((k_grid(ik) + Gamma_redistr*0.5d0)*(1d0+(r/OmegaRet2(Tret+1-it))*(1d0-tk)) &
                            +lumpsum*0.5d0+y*(1d0-tax_labors(y)-tSS_employee(y))-P3*(1d0+tc))/(1d0+mu)
                    if(dum2<0.0001d0) then
                        V3=-999999999d0
                    else
                        V3=Uc(P3)-chims*(dum4**(1d0+etam))/(1d0+etam)-fcm(2,ifc)
                        pnt1 = [dum2, exp_grid(ix,T+Tret+1-it)+1d0]
                        vnext = bilin_interp( &
                                k_grid, exp_grid_dum, &
                                evs_ret(2,j,:,:,iam,ium,Tret+2-it,ifc), &
                                nk, nexp, pnt1)
                        V3=V3+beta*OmegaRet(Tret-it+1)*vnext
                    end if
                    if (V2 < V3) then
                        P1=P2
                    else
                        P4=P3
                    end if
                    if((P4-P1)<1d-8) exit
                end do
                Ves=V2
                kes=dum2
                nes=dum4
                ces=P2

            end if



            vus=vs_ret(1,j,ik,ix,iam,ium,Tret-it+1,ifc)
            kus=ks_ret(1,j,ik,ix,iam,ium,Tret-it+1,ifc)
            nus=ns_ret(1,j,ik,ix,iam,ium,Tret-it+1,ifc)
            cus=cs_ret(1,j,ik,ix,iam,ium,Tret-it+1,ifc)

            if (ves >= vus) then
                vs_ret(2,j,ik,ix,iam,ium,Tret-it+1,ifc)=ves
                cs_ret(2,j,ik,ix,iam,ium,Tret-it+1,ifc)=ces
                ks_ret(2,j,ik,ix,iam,ium,Tret-it+1,ifc)=kes
                ns_ret(2,j,ik,ix,iam,ium,Tret-it+1,ifc)=nes
            else
                vs_ret(2,j,ik,ix,iam,ium,Tret-it+1,ifc)=vus
                cs_ret(2,j,ik,ix,iam,ium,Tret-it+1,ifc)=cus
                ks_ret(2,j,ik,ix,iam,ium,Tret-it+1,ifc)=kus
                ns_ret(2,j,ik,ix,iam,ium,Tret-it+1,ifc)=nus
            end if

        end do
    end do

end subroutine SolveInRetirement2
    
subroutine SolveInRetirement3(counter)
    ! Woman and Man can is working
    use Model_Parameters
    use PolicyFunctions
    use glob0
    use Utilities

    implicit none

    integer, INTENT(IN) :: counter
    integer :: ix,ixm,ixd,iam,ium,iaf,iuf,iu2,iu3,tprint,ikd,j,ifc,ifcm,ik
    real(8) :: ce,cu,ke,ku,nem,nef,num,nuf,ve,vu
    real(8) :: ces,cus,kes,kus,nes,nus,ves,vus
    integer :: NEQ=0, IERSVR=0, IPACT=0, ISACT=0
    real(8) :: c2, MU2, d1, d2, vp(nu),dum3,dum4,dum5,dum6,y
    real(8) :: ACC=0.0001d0,ERREL=0.0001d0
    real(8) :: P1,P2,P3,P4,V2,V3,dum2,pnt2(3),pnt1(2)
    real(8) :: vnext, exp_grid_dum(nexp), INTERP2D(nk,nexp)

    dum3=((counter*1d0)/(nu*na*1d0))-0.00001d0
    ik=int(dum3)+1
    dum3=(((counter-(ik-1)*nu*na)*1d0)/(na*1d0))-0.00001d0
    ium=int(dum3)+1
    iam=counter-(ik-1)*na*nu-(ium-1)*nu

    if(it>1) then
        exp_grid_dum=exp_grid(:,T+Tret+2-it)
    end if

    ce=0d0
    cu=0d0
    ke=0d0
    ku=0d0
    nem=0d0
    nef=0d0
    num=0d0
    nuf=0d0
    ve=0d0
    vu=0d0

    ces=0d0
    cus=0d0
    kes=0d0
    kus=0d0
    nes=0d0
    nus=0d0
    ves=0d0
    vus=0d0

    !Married man and woman is working

    do ifc=1,nfc
        do ifcm=1,nfc
            do iaf = 1, na
                do iuf = 1, nu
                    do ix = 1, nexp
                        do ixm = 1, nexp
                            wagem = wage(1,a(1,iam),exp_grid(ixm,T-it),u(1,ium))/(1d0+t_employer)
                            wagef = wage(2,a(2,iaf),exp_grid(ix,T-it),u(2,iuf))/(1d0+t_employer)

                            if(it==1) then
                                ! Solve the very last period problem:           
                                ke = 0d0
                                ce = ((k_grid(ik) + Gamma_redistr)*(((1d0 + r)/OmegaRet2(Tret+1-it))*(1d0-tk)) + Psi_pension+lumpsum)/(1d0+tc)
                                nef = 0d0
                                nem = 0d0
                                ve = Uc(ce)
                            else
                                ! Solve the Tret-1 to 1st period of retirement:
                                !Finding optimal capital by golden search
                                P1=0.01d0
                                P4=min((k_grid(nk)-0.001d0)/(1d0+tc),((k_grid(ik) + Gamma_redistr)*(1d0+(r/OmegaRet2(Tret+1-it))*(1d0-tk))+lumpsum+wagem*(1d0-tax_labor(wagem)-tSS_employee(wagem))+wagef*(1d0-tax_labor(wagef)-tSS_employee(wagef)))/(1d0+tc))
                                do
                                    P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                                    P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)

                                    pnt2 = [P2, wagem, wagef]
                                    dum4 = min(max(trilin_interp(c_grid, wage_grid, wage_grid, laborm, nc, nw, nw, pnt2),1d-10),1d0)
                                    dum5 = min(max(trilin_interp(c_grid, wage_grid, wage_grid, laborf, nc, nw, nw, pnt2),1d-10),1d0)
                                    y=dum4*wagem+dum5*wagef
                                    dum2=((k_grid(ik) + Gamma_redistr)*(1d0+(r/OmegaRet2(Tret+1-it))*(1d0-tk))+lumpsum+y*(1d0-tax_labor(y)-tSS_employee(y))-P2*(1d0+tc))/(1d0+mu)
                                    if(dum2<0.0001d0) then
                                        V2=-999999999d0
                                    else
                                        V2=Uc(P2)+Ul(dum4,dum5)-fc(1,ifc)-fcm(1,ifcm)

                                        pnt2 = [dum2,exp_grid(ix,T+Tret+1-it)+1d0,exp_grid(ixm,T+Tret+1-it)+1d0]
                                        vnext = trilin_interp( &
                                                    k_grid,exp_grid(:,T+Tret-it+2),exp_grid(:,T+Tret-it+2), &
                                                    ev_ret(2,2,:,:,:,iam,ium,iaf,iuf,Tret+2-it,ifc,ifcm), &
                                                    nk, nexp, nexp, pnt2 )

                                        !vnext = D_QD3VL(dum2,exp_grid(ix,T+Tret+1-it)+1d0,exp_grid(ixm,T+Tret+1-it)+1d0,k_grid,exp_grid(:,T+Tret-it+2),exp_grid(:,T+Tret-it+2),ev_ret(2,2,:,:,:,iam,ium,iaf,iuf,Tret+2-it,ifc,ifcm))
                                        V2=V2+beta*OmegaRet(Tret-it+1)*vnext
                                    end if

                                    pnt2 = [P2, wagem, wagef]
                                    dum4 = min(max(trilin_interp(c_grid, wage_grid, wage_grid, laborm, nc, nw, nw, pnt2),1d-10),1d0)
                                    dum5 = min(max(trilin_interp(c_grid, wage_grid, wage_grid, laborf, nc, nw, nw, pnt2),1d-10),1d0)
                                    y=dum4*wagem+dum5*wagef
                                    dum2=((k_grid(ik) + Gamma_redistr)*(1d0+(r/OmegaRet2(Tret+1-it))*(1d0-tk))+lumpsum+y*(1d0-tax_labor(y)-tSS_employee(y))-P3*(1d0+tc))/(1d0+mu)
                                    if(dum2<0.0001d0) then
                                        V3=-999999999d0
                                    else
                                        V3=Uc(P3)+Ul(dum4,dum5)-fc(1,ifc)-fcm(1,ifcm)

                                        pnt2 = [dum2,exp_grid(ix,T+Tret+1-it)+1d0,exp_grid(ixm,T+Tret+1-it)+1d0]
                                        vnext = trilin_interp( &
                                                    k_grid,exp_grid(:,T+Tret-it+2),exp_grid(:,T+Tret-it+2), &
                                                    ev_ret(2,2,:,:,:,iam,ium,iaf,iuf,Tret+2-it,ifc,ifcm), &
                                                    nk, nexp, nexp, pnt2 )

                                        !vnext = D_QD3VL(dum2,exp_grid(ix,T+Tret+1-it)+1d0,exp_grid(ixm,T+Tret+1-it)+1d0,k_grid,exp_grid(:,T+Tret-it+2),exp_grid(:,T+Tret-it+2),ev_ret(2,2,:,:,:,iam,ium,iaf,iuf,Tret+2-it,ifc,ifcm))
                                        V3=V3+beta*OmegaRet(Tret-it+1)*vnext
                                    end if

                                    if (V2 < V3) then
                                        P1=P2
                                    else
                                        P4=P3
                                    end if
                                    if((P4-P1)<1d-8) exit
                                end do
                                Ve=V2
                                ke=dum2
                                nem=dum4
                                nef=dum5
                                ce=P2
                            end if

                            vu=v_ret(1,1,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)
                            ku=k_ret(1,1,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)
                            nuf=0d0
                            num=0d0
                            cu=c_ret(1,1,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)

                            if (ve >= vu) then
                                v_ret(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)=ve
                                c_ret(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)=ce
                                k_ret(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)=ke
                                nf_ret(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)=nef
                                nm_ret(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)=nem
                            else
                                v_ret(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)=vu
                                c_ret(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)=cu
                                k_ret(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)=ku
                                nf_ret(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)=nuf
                                nm_ret(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)=num
                            end if

                            if (v_ret(2,1,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm) >= v_ret(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)) then
                                v_ret(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)=v_ret(2,1,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)
                                c_ret(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)=c_ret(2,1,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)
                                k_ret(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)=k_ret(2,1,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)
                                nf_ret(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)=nf_ret(2,1,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)
                                nm_ret(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)=nm_ret(2,1,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)
                            end if

                            if (v_ret(1,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm) >= v_ret(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)) then
                                v_ret(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)=v_ret(1,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)
                                c_ret(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)=c_ret(1,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)
                                k_ret(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)=k_ret(1,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)
                                nf_ret(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)=nf_ret(1,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)
                                nm_ret(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)= &
                                    nm_ret(1,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm)
                            end if

                        end do
                    end do
                end do
            end do
        end do
    end do


end subroutine SolveInRetirement3
