subroutine Solvelastactive(counter)

    !This subroutine computes optimal policies at age 64
    use Model_Parameters
    use PolicyFunctions
    use glob0
    use Utilities

    implicit none

    integer, INTENT(IN) :: counter
    integer :: ik,ix,ixm,ixd,iam,ium,iaf,iuf,iu2,iu3,j,ifc,ifcm
    real(8) :: ce,cu,ke,ku,nem,nef,num,nuf,ve,vu
    real(8) :: ces,cus,kes,kus,nes,nus,ves,vus
    real(8), dimension (:), allocatable :: Expdum
    integer :: NEQ=0, IERSVR=0, IPACT=0, ISACT=0
    real(8) :: c2, MU2, d1, d2, vp(nu),dum3,dum4,dum5,dum6,y
    real(8) :: ACC=0.0001d0,ERREL=0.0001d0
    real(8) :: P1,P2,P3,P4,V2,V3,dum2,pnt2(3),pnt1(2)
    real(8) :: vnext, exp_grid_dum(nexp), INTERP2D(nk,nexp)

    real(8) :: dd1, dd2, dd3, dd4, dd5


    !Assigning the grid points
    dum3=((counter*1d0)/(nu*na*1d0))-0.00001d0
    ik=int(dum3)+1
    dum3=(((counter-(ik-1)*nu*na)*1d0)/(na*1d0))-0.00001d0
    ium=int(dum3)+1
    iam=counter-(ik-1)*na*nu-(ium-1)*nu

    exp_grid_dum=exp_grid(:,T+1)

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


    ! Last period of active life
    !if employed
    do ifc=1,nfc
        do ifcm=1,nfc
            do iaf = 1, na
                do iuf = 1, nu
                    do ix = 1, nexp
                        do ixm = 1, nexp
                            !Print *,'ix is',ix
                            !Finding optimal capital by golden search
                            wagem = wage(1,a(1,iam),exp_grid(ixm,T),u(1,ium))/(1d0+t_employer)
                            wagef = wage(2,a(2,iaf),exp_grid(ix,T),u(2,iuf))/(1d0+t_employer)
                            P1=0.01d0
                            P4=min((k_grid(nk)-0.001d0)/(1d0+tc),((k_grid(ik) + Gamma_redistr)*(1d0+r*(1d0-tk))+lumpsum+(wagem+wagef)*(1d0-tax_labor(wagem+wagef)-tSS_employee(wagem+wagef)))/(1d0+tc))
                            dd1 = (k_grid(ik) + Gamma_redistr)*(1d0+r*(1d0-tk))
                            dd2 = tax_labor(wagem+wagef)
                            dd3 = wagem+wagef
                            dd4 = tSS_employee(wagem+wagef)
                            do
                                P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                                P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)

                                pnt2 = (/P2, wagem, wagef/)
                                dum4 = min(max(trilin_interp(c_grid, wage_grid, wage_grid, laborm, nc, nw, nw, pnt2),1d-10),1d0)
                                dum5 = min(max(trilin_interp(c_grid, wage_grid, wage_grid, laborf, nc, nw, nw, pnt2),1d-10),1d0)
                                y=dum4*wagem+dum5*wagef
                                dum2=((k_grid(ik) + Gamma_redistr)*(1d0+r*(1d0-tk))+lumpsum+y*(1d0-tax_labor(y)-tSS_employee(y))-P2*(1d0+tc))/(1d0+mu)
                                if(dum2<0.0001d0) then
                                    V2=-999999999d0
                                else
                                    V2=Uc(P2)+Ul(dum4,dum5)-fc(1,ifc)-fcm(1,ifcm)

                                    pnt2 = [dum2,exp_grid(ix,T)+1d0,exp_grid(ixm,T)+1d0]
                                    vnext = trilin_interp(k_grid,exp_grid(:,T+1),exp_grid(:,T+1), &
                                                        ev_ret(2,2,:,:,:,iam,ium,iaf,iuf,1,ifc,ifcm), &
                                                        nk, nexp, nexp, pnt2)

                                    !vnext = D_QD3VL(dum2,exp_grid(ix,T)+1d0,exp_grid(ixm,T)+1d0,k_grid,exp_grid(:,T+1),exp_grid(:,T+1),ev_ret(2,2,:,:,:,iam,ium,iaf,iuf,1,ifc,ifcm))
                                    V2=V2+beta*OmegaActive(T)*vnext
                                end if

                                pnt2 = [P3, wagem, wagef]
                                dum4 = min(max(trilin_interp(c_grid, wage_grid, wage_grid, laborm, nc, nw, nw, pnt2),1d-10),1d0)
                                dum5 = min(max(trilin_interp(c_grid, wage_grid, wage_grid, laborf, nc, nw, nw, pnt2),1d-10),1d0)
                                y=dum4*wagem+dum5*wagef
                                dum2=((k_grid(ik) + Gamma_redistr)*(1d0+r*(1d0-tk))+lumpsum+y*(1d0-tax_labor(y)-tSS_employee(y))-P3*(1d0+tc))/(1d0+mu)
                                if(dum2<0.0001d0) then
                                    V3=-999999999d0
                                else
                                    V3=Uc(P3)+Ul(dum4,dum5)-fc(1,ifc)-fcm(1,ifcm)

                                    pnt2 = [dum2,exp_grid(ix,T)+1d0,exp_grid(ixm,T)+1d0]
                                    vnext = trilin_interp(k_grid,exp_grid(:,T+1),exp_grid(:,T+1), &
                                                        ev_ret(2,2,:,:,:,iam,ium,iaf,iuf,1,ifc,ifcm), &
                                                        nk, nexp, nexp, pnt2)

                                    !vnext = D_QD3VL(dum2,exp_grid(ix,T)+1d0,exp_grid(ixm,T)+1d0,k_grid,exp_grid(:,T+1),exp_grid(:,T+1),ev_ret(2,2,:,:,:,iam,ium,iaf,iuf,1,ifc,ifcm))
                                    V3=V3+beta*OmegaActive(T)*vnext
                                end if

                                if (V2 < V3) then
                                    P1=P2
                                else
                                    P4=P3
                                end if
                                if((P4-P1)<1d-6) exit
                            end do

                            Ve=V2
                            ke=dum2
                            nem=dum4
                            nef=dum5
                            ce=P2

                            if(nem<0d0) then
                                Print *,'ik is',ik
                                Print *,'ix is',ix
                                Print *,'iam is',iam
                                Print *,'ium is',ium
                                Print *,'iaf is',iaf
                                Print *,'iuf is',iuf
                                Print *,'num is',num
                                Print *,'nuf is',nuf
                                read(*,*)
                            end if

                            if(nef<0d0) then
                                Print *,'ik is',ik
                                Print *,'ix is',ix
                                Print *,'iam is',iam
                                Print *,'ium is',ium
                                Print *,'iaf is',iaf
                                Print *,'iuf is',iuf
                                Print *,'nem is',nem
                                Print *,'nef is',nef
                                read(*,*)
                            end if


                            !If female unemployed

                            P1=0.01d0
                            P4=min((k_grid(nk)-0.001d0)/(1d0+tc),((k_grid(ik) + Gamma_redistr)*(1d0+r*(1d0-tk))+lumpsum+Unemp_benefit+(wagem)*(1d0-tax_labor(wagem)-tSS_employee(wagem)))/(1d0+tc))
                            do
                                P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                                P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)

                                pnt1 = (/P2, wagem/)
                                dum4 = min(max(bilin_interp(c_grid, wage_grid, labormwork, nc, nw, pnt1),0d0),1d0)
                                y=dum4*wagem
                                dum2=((k_grid(ik) + Gamma_redistr)*(1d0+r*(1d0-tk))+lumpsum+Unemp_benefit+y*(1d0-tax_labor(y)-tSS_employee(y))-P2*(1d0+tc))/(1d0+mu)
                                if(dum2<0.0001d0) then
                                    V2=-999999999d0
                                else
                                    V2=Uc(P2)+Ul(dum4,0d0)-fcm(1,ifcm)
                                    pnt2 = [dum2,exp_grid(ix,T)*(1d0-deltaexp),exp_grid(ixm,T)+1d0]
                                    vnext = trilin_interp(k_grid,exp_grid(:,T+1),exp_grid(:,T+1), &
                                                        ev_ret(2,2,:,:,:,iam,ium,iaf,iuf,1,ifc,ifcm), &
                                                        nk, nexp, nexp, pnt2)

                                    !vnext = D_QD3VL(dum2,exp_grid(ix,T)*(1d0-deltaexp),exp_grid(ixm,T)+1d0,k_grid,exp_grid(:,T+1),exp_grid(:,T+1),ev_ret(2,2,:,:,:,iam,ium,iaf,iuf,1,ifc,ifcm))
                                    V2=V2+beta*OmegaActive(T)*vnext
                                end if

                                pnt1 = (/P3, wagem/)
                                dum4 = min(max(bilin_interp(c_grid, wage_grid, labormwork, nc, nw, pnt1),0d0),1d0)
                                y=dum4*wagem
                                dum2=((k_grid(ik) + Gamma_redistr)*(1d0+r*(1d0-tk))+lumpsum+Unemp_benefit+y*(1d0-tax_labor(y)-tSS_employee(y))-P3*(1d0+tc))/(1d0+mu)
                                if(dum2<0.0001d0) then
                                    V3=-999999999d0
                                else
                                    V3=Uc(P3)+Ul(dum4,0d0)-fcm(1,ifcm)
                                    pnt2 = [dum2,exp_grid(ix,T)*(1d0-deltaexp),exp_grid(ixm,T)+1d0]
                                    vnext = trilin_interp(k_grid,exp_grid(:,T+1),exp_grid(:,T+1), &
                                                        ev_ret(2,2,:,:,:,iam,ium,iaf,iuf,1,ifc,ifcm), &
                                                        nk, nexp, nexp, pnt2)

                                    !vnext = D_QD3VL(dum2,exp_grid(ix,T)*(1d0-deltaexp),exp_grid(ixm,T)+1d0,k_grid,exp_grid(:,T+1),exp_grid(:,T+1),ev_ret(2,2,:,:,:,iam,ium,iaf,iuf,1,ifc,ifcm))
                                    V3=V3+beta*OmegaActive(T)*vnext
                                end if

                                if (V2 < V3) then
                                    P1=P2
                                else
                                    P4=P3
                                end if
                                if((P4-P1)<1d-6) exit
                            end do

                            Vu=V2
                            ku=dum2
                            num=dum4
                            nuf=0d0
                            cu=P2

                            if(num<0d0) then
                                Print *,'ik is',ik
                                Print *,'ix is',ix
                                Print *,'iam is',iam
                                Print *,'ium is',ium
                                Print *,'iaf is',iaf
                                Print *,'iuf is',iuf
                                Print *,'num is',num
                                Print *,'nuf is',nuf
                                read(*,*)
                            end if

                            if(nuf<0d0) then
                                Print *,'ik is',ik
                                Print *,'ix is',ix
                                Print *,'iam is',iam
                                Print *,'ium is',ium
                                Print *,'iaf is',iaf
                                Print *,'iuf is',iuf
                                Print *,'num is',num
                                Print *,'nuf is',nuf
                                read(*,*)
                            end if

                            !If male unemployed

                            P1=0.01d0
                            P4=min((k_grid(nk)-0.001d0)/(1d0+tc),((k_grid(ik) + Gamma_redistr)*(1d0+r*(1d0-tk))+lumpsum+Unemp_benefit+(wagef)*(1d0-tax_labor(wagef)-tSS_employee(wagef)))/(1d0+tc))
                            do
                                P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                                P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)

                                pnt1 = (/P2, wagef/)
                                dum4 = min(max(bilin_interp(c_grid, wage_grid, laborfwork, nc, nw, pnt1),0d0),1d0)
                                y=dum4*wagef
                                dum2=((k_grid(ik) + Gamma_redistr)*(1d0+r*(1d0-tk))+lumpsum+Unemp_benefit+y*(1d0-tax_labor(y)-tSS_employee(y))-P2*(1d0+tc))/(1d0+mu)
                                if(dum2<0.0001d0) then
                                    V2=-999999999d0
                                else
                                    V2=Uc(P2)+Ul(0d0,dum4)-fc(1,ifc)
                                    pnt2 = [dum2,exp_grid(ix,T)+1d0,exp_grid(ixm,T)*(1d0-deltaexp)]
                                    vnext = trilin_interp(k_grid,exp_grid(:,T+1),exp_grid(:,T+1), &
                                                        ev_ret(2,2,:,:,:,iam,ium,iaf,iuf,1,ifc,ifcm), &
                                                        nk, nexp, nexp, pnt2)
                                    !vnext = D_QD3VL(dum2,exp_grid(ix,T)+1d0,exp_grid(ixm,T)*(1d0-deltaexp),k_grid,exp_grid(:,T+1),exp_grid(:,T+1),ev_ret(2,2,:,:,:,iam,ium,iaf,iuf,1,ifc,ifcm))
                                    V2=V2+beta*OmegaActive(T)*vnext
                                end if

                                pnt1 = (/P3, wagef/)
                                dum4 = min(max(bilin_interp(c_grid, wage_grid, laborfwork, nc, nw, pnt1),0d0),1d0)
                                y=dum4*wagef
                                dum2=((k_grid(ik) + Gamma_redistr)*(1d0+r*(1d0-tk))+lumpsum+Unemp_benefit+y*(1d0-tax_labor(y)-tSS_employee(y))-P3*(1d0+tc))/(1d0+mu)
                                if(dum2<0.0001d0) then
                                    V3=-999999999d0
                                else
                                    V3=Uc(P3)+Ul(0d0,dum4)-fc(1,ifc)
                                    pnt2 = [dum2,exp_grid(ix,T)+1d0,exp_grid(ixm,T)*(1d0-deltaexp)]
                                    vnext = trilin_interp(k_grid,exp_grid(:,T+1),exp_grid(:,T+1), &
                                                        ev_ret(2,2,:,:,:,iam,ium,iaf,iuf,1,ifc,ifcm), &
                                                        nk, nexp, nexp, pnt2)
                                    !vnext = D_QD3VL(dum2,exp_grid(ix,T)+1d0,exp_grid(ixm,T)*(1d0-deltaexp),k_grid,exp_grid(:,T+1),exp_grid(:,T+1),ev_ret(2,2,:,:,:,iam,ium,iaf,iuf,1,ifc,ifcm))
                                    V3=V3+beta*OmegaActive(T)*vnext
                                end if

                                if (V2 < V3) then
                                    P1=P2
                                else
                                    P4=P3
                                end if
                                if((P4-P1)<1d-6) exit
                            end do

                            if (V2 >= vu) then
                                Vu=V2
                                ku=dum2
                                num=0d0
                                nuf=dum4
                                cu=P2
                            end if

                            if(num<0d0) then
                                Print *,'ik is',ik
                                Print *,'ix is',ix
                                Print *,'iam is',iam
                                Print *,'ium is',ium
                                Print *,'iaf is',iaf
                                Print *,'iuf is',iuf
                                Print *,'num is',num
                                Print *,'nuf is',nuf
                                read(*,*)
                            end if

                            if(nuf<0d0) then
                                Print *,'ik is',ik
                                Print *,'ix is',ix
                                Print *,'iam is',iam
                                Print *,'ium is',ium
                                Print *,'iaf is',iaf
                                Print *,'iuf is',iuf
                                Print *,'num is',num
                                Print *,'nuf is',nuf
                                read(*,*)
                            end if


                            !If both spouses unemployed

                            P1=0.01d0
                            P4=min((k_grid(nk)-0.001d0)/(1d0+tc),((k_grid(ik) + Gamma_redistr)*(1d0+r*(1d0-tk))+lumpsum+2d0*Unemp_benefit)/(1d0+tc))
                            do
                                P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                                P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)

                                dum2=((k_grid(ik) + Gamma_redistr)*(1d0+r*(1d0-tk))+lumpsum+2d0*Unemp_benefit-P2*(1d0+tc))/(1d0+mu)
                                if(dum2<0.0001d0) then
                                    V2=-999999999d0
                                else
                                    V2=Uc(P2)+Ul(0d0,0d0)
                                    pnt2 = [dum2,exp_grid(ix,T)*(1d0-deltaexp),exp_grid(ixm,T)*(1d0-deltaexp)]
                                    vnext = trilin_interp(k_grid,exp_grid(:,T+1),exp_grid(:,T+1), &
                                                        ev_ret(2,2,:,:,:,iam,ium,iaf,iuf,1,ifc,ifcm), &
                                                        nk, nexp, nexp, pnt2)
                                    !vnext = D_QD3VL(dum2,exp_grid(ix,T)*(1d0-deltaexp),exp_grid(ixm,T)*(1d0-deltaexp),k_grid,exp_grid(:,T+1),exp_grid(:,T+1),ev_ret(2,2,:,:,:,iam,ium,iaf,iuf,1,ifc,ifcm))
                                    V2=V2+beta*OmegaActive(T)*vnext
                                end if

                                dum2=((k_grid(ik) + Gamma_redistr)*(1d0+r*(1d0-tk))+lumpsum+2d0*Unemp_benefit-P3*(1d0+tc))/(1d0+mu)
                                if(dum2<0.0001d0) then
                                    V3=-999999999d0
                                else
                                    V3=Uc(P3)+Ul(0d0,0d0)
                                    pnt2 = [dum2,exp_grid(ix,T)*(1d0-deltaexp),exp_grid(ixm,T)*(1d0-deltaexp)]
                                    vnext = trilin_interp(k_grid,exp_grid(:,T+1),exp_grid(:,T+1), &
                                                        ev_ret(2,2,:,:,:,iam,ium,iaf,iuf,1,ifc,ifcm), &
                                                        nk, nexp, nexp, pnt2)

                                    !vnext = D_QD3VL(dum2,exp_grid(ix,T)*(1d0-deltaexp),exp_grid(ixm,T)*(1d0-deltaexp),k_grid,exp_grid(:,T+1),exp_grid(:,T+1),ev_ret(2,2,:,:,:,iam,ium,iaf,iuf,1,ifc,ifcm))
                                    V3=V3+beta*OmegaActive(T)*vnext
                                end if

                                if (V2 < V3) then
                                    P1=P2
                                else
                                    P4=P3
                                end if
                                if((P4-P1)<1d-6) exit
                            end do

                            if (V2 >= vu) then
                                Vu=V2
                                ku=dum2
                                num=0d0
                                nuf=0d0
                                cu=P2
                            end if


                            if (ve >= vu) then
                                v(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm)=ve
                                c(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm)=ce
                                k(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm)=ke
                                nm(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm)=nem
                                nf(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm)=nef
                            else
                                v(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm)=vu
                                c(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm)=cu
                                k(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm)=ku
                                nm(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm)=num
                                nf(ik,ix,ixm,iam,ium,iaf,iuf,T,ifc,ifcm)=nuf
                            end if



                        end do
                    end do
                end do
            end do
        end do
    end do

    !Singles    
    j=2
    do ifc=1,nfc
        do ix = 1, nexp
            !Print *,'ix is',ix
            !Finding optimal capital by golden search
            wagef = wage(2,a(2,iam),exp_grid(ix,T),u(2,ium))/(1d0+t_employer)
            P1=0.01d0
            P4=min((k_grid(nk)-0.001d0)/(1d0+tc),((k_grid(ik) + Gamma_redistr*0.5d0)*(1d0+r*(1d0-tk))+lumpsum*0.5d0+wagef*(1d0-tax_labors(wagef)-tSS_employee(wagef)))/(1d0+tc))
            do
                P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)

                pnt1 = (/P2, wagef/)
                dum4 = min(max(bilin_interp(c_grid, wage_grid, laborsinglef, nc, nw, pnt1),0d0),1d0)
                y=dum4*wagef
                dum2=((k_grid(ik) + Gamma_redistr*0.5d0)*(1d0+r*(1d0-tk))+lumpsum*0.5d0+y*(1d0-tax_labors(y)-tSS_employee(y))-P2*(1d0+tc))/(1d0+mu)
                if(dum2<0.0001d0) then
                    V2=-999999999d0
                else
                    V2=Uc(P2)-chifs*(dum4**(1d0+etaf))/(1d0+etaf)-fc(2,ifc)
                    pnt1 = [dum2,exp_grid(ix,T)+1d0]
                    vnext = bilin_interp(k_grid,exp_grid(:,T+1),evs_ret(2,j,:,:,iam,ium,1,ifc),nk,nexp, pnt1)
                    !vnext = D_QD2VL(dum2,exp_grid(ix,T)+1d0,k_grid,exp_grid(:,T+1),evs_ret(2,j,:,:,iam,ium,1,ifc))
                    V2=V2+beta*OmegaActive(T)*vnext
                end if

                pnt1 = (/P3, wagef/)
                dum4 = min(max(bilin_interp(c_grid, wage_grid, laborsinglef, nc, nw, pnt1),0d0),1d0)
                y=dum4*wagef
                dum2=((k_grid(ik) + Gamma_redistr*0.5d0)*(1d0+r*(1d0-tk))+lumpsum*0.5d0+y*(1d0-tax_labors(y)-tSS_employee(y))-P3*(1d0+tc))/(1d0+mu)
                if(dum2<0.0001d0) then
                    V3=-999999999d0
                else
                    pnt1 = [dum2,exp_grid(ix,T)+1d0]
                    vnext = bilin_interp(k_grid,exp_grid(:,T+1),evs_ret(2,j,:,:,iam,ium,1,ifc),nk,nexp, pnt1)
                    !vnext = D_QD2VL(dum2,exp_grid(ix,T)+1d0,k_grid,exp_grid(:,T+1),evs_ret(2,j,:,:,iam,ium,1,ifc))
                    V3=Uc(P3)-chifs*(dum4**(1d0+etaf))/(1d0+etaf)-fc(2,ifc)
                    V3=V3+beta*OmegaActive(T)*vnext
                end if

                if (V2 < V3) then
                    P1=P2
                else
                    P4=P3
                end if
                if((P4-P1)<1d-6) exit
            end do

            Ves=V2
            kes=dum2
            nes=dum4
            ces=P2

            if(nes<0d0) then
                Print *,'ik is',ik
                Print *,'ix is',ix
                Print *,'iam is',iam
                Print *,'ium is',ium
                Print *,'nes is',nes
                read(*,*)
            end if

            !If single female unemployed
            !Print *,'ix is',ix
            !Finding optimal capital by golden search
            P1=0.01d0
            P4=min((k_grid(nk)-0.001d0)/(1d0+tc),((k_grid(ik) + Gamma_redistr*0.5d0)*(1d0+r*(1d0-tk))+lumpsum*0.5d0+Unemp_benefit)/(1d0+tc))
            do
                P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)

                dum2=((k_grid(ik) + Gamma_redistr*0.5d0)*(1d0+r*(1d0-tk))+lumpsum*0.5d0+Unemp_benefit-P2*(1d0+tc))/(1d0+mu)
                if(dum2<0.0001d0) then
                    V2=-999999999d0
                else
                    pnt1 = [dum2,exp_grid(ix,T)*(1d0-deltaexp)]
                    vnext = bilin_interp(k_grid,exp_grid(:,T+1),evs_ret(2,j,:,:,iam,ium,1,ifc),nk,nexp, pnt1)
                    !vnext = D_QD2VL(dum2,exp_grid(ix,T)*(1d0-deltaexp),k_grid,exp_grid(:,T+1),evs_ret(2,j,:,:,iam,ium,1,ifc))
                    V2=Uc(P2)
                    V2=V2+beta*OmegaActive(T)*vnext
                end if

                dum2=((k_grid(ik) + Gamma_redistr*0.5d0)*(1d0+r*(1d0-tk))+lumpsum*0.5d0+Unemp_benefit-P3*(1d0+tc))/(1d0+mu)
                if(dum2<0.0001d0) then
                    V3=-999999999d0
                else
                    pnt1 = [dum2,exp_grid(ix,T)*(1d0-deltaexp)]
                    vnext = bilin_interp(k_grid,exp_grid(:,T+1),evs_ret(2,j,:,:,iam,ium,1,ifc),nk,nexp, pnt1)
                    !vnext = D_QD2VL(dum2,exp_grid(ix,T)*(1d0-deltaexp),k_grid,exp_grid(:,T+1),evs_ret(2,j,:,:,iam,ium,1,ifc))
                    V3=Uc(P3)
                    V3=V3+beta*OmegaActive(T)*vnext
                end if

                if (V2 < V3) then
                    P1=P2
                else
                    P4=P3
                end if
                if((P4-P1)<1d-6) exit
            end do

            Vus=V2
            kus=dum2
            nus=0d0
            cus=P2

            if (ves >= vus) then
                vs(j,ik,ix,iam,ium,T,ifc)=ves
                cs(j,ik,ix,iam,ium,T,ifc)=ces
                ks(j,ik,ix,iam,ium,T,ifc)=kes
                ns(j,ik,ix,iam,ium,T,ifc)=nes
            else
                vs(j,ik,ix,iam,ium,T,ifc)=vus
                cs(j,ik,ix,iam,ium,T,ifc)=cus
                ks(j,ik,ix,iam,ium,T,ifc)=kus
                ns(j,ik,ix,iam,ium,T,ifc)=nus
            end if
        end do
    end do


    !Men   
    j=1
    do ifc=1,nfc
        do ix = 1, nexp
            !Print *,'ix is',ix
            !Finding optimal capital by golden search
            wagem = wage(1,a(1,iam),exp_grid(ix,T),u(1,ium))/(1d0+t_employer)
            P1=0.01d0
            P4=min((k_grid(nk)-0.001d0)/(1d0+tc),((k_grid(ik) + Gamma_redistr*0.5d0)*(1d0+r*(1d0-tk))+lumpsum*0.5d0+wagem*(1d0-tax_labors(wagem)-tSS_employee(wagem)))/(1d0+tc))
            do
                P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)

                pnt1 = (/P2, wagem/)
                dum4 = min(max(bilin_interp(c_grid, wage_grid, laborsinglem, nc, nw, pnt1),0d0),1d0)
                y=dum4*wagem
                dum2=((k_grid(ik) + Gamma_redistr*0.5d0)*(1d0+r*(1d0-tk))+lumpsum*0.5d0+y*(1d0-tax_labors(y)-tSS_employee(y))-P2*(1d0+tc))/(1d0+mu)
                if(dum2<0.0001d0) then
                    V2=-999999999d0
                else
                    V2=Uc(P2)-chims*(dum4**(1d0+etam))/(1d0+etam)-fcm(2,ifc)
                    pnt1 = [dum2,exp_grid(ix,T)+1d0]
                    vnext = bilin_interp(k_grid,exp_grid(:,T+1),evs_ret(2,j,:,:,iam,ium,1,ifc),nk,nexp, pnt1)
                    !vnext = D_QD2VL(dum2,exp_grid(ix,T)+1d0,k_grid,exp_grid(:,T+1),evs_ret(2,j,:,:,iam,ium,1,ifc))
                    V2=V2+beta*OmegaActive(T)*vnext
                end if

                pnt1 = (/P3, wagem/)
                dum4 = min(max(bilin_interp(c_grid, wage_grid, laborsinglem, nc, nw, pnt1),0d0),1d0)
                y=dum4*wagem
                dum2=((k_grid(ik) + Gamma_redistr*0.5d0)*(1d0+r*(1d0-tk))+lumpsum*0.5d0+y*(1d0-tax_labors(y)-tSS_employee(y))-P3*(1d0+tc))/(1d0+mu)
                if(dum2<0.0001d0) then
                    V3=-999999999d0
                else
                    pnt1 = [dum2,exp_grid(ix,T)+1d0]
                    vnext = bilin_interp(k_grid,exp_grid(:,T+1),evs_ret(2,j,:,:,iam,ium,1,ifc),nk,nexp, pnt1)
                    !vnext = D_QD2VL(dum2,exp_grid(ix,T)+1d0,k_grid,exp_grid(:,T+1),evs_ret(2,j,:,:,iam,ium,1,ifc))
                    V3=Uc(P3)-chims*(dum4**(1d0+etam))/(1d0+etam)-fcm(2,ifc)
                    V3=V3+beta*OmegaActive(T)*vnext
                end if

                if (V2 < V3) then
                    P1=P2
                else
                    P4=P3
                end if
                if((P4-P1)<1d-6) exit
            end do

            Ves=V2
            kes=dum2
            nes=dum4
            ces=P2

            if(nes<0d0) then
                Print *,'ik is',ik
                Print *,'ix is',ix
                Print *,'iam is',iam
                Print *,'ium is',ium
                Print *,'nes is',nes
                read(*,*)
            end if

            !If single male unemployed
            !Print *,'ix is',ix
            !Finding optimal capital by golden search
            P1=0.01d0
            P4=min((k_grid(nk)-0.001d0)/(1d0+tc),((k_grid(ik) + Gamma_redistr*0.5d0)*(1d0+r*(1d0-tk))+lumpsum*0.5d0+Unemp_benefit)/(1d0+tc))
            do
                P2 = P1 + ((3.0-sqrt(5.0))/2.0)*(P4-P1)
                P3 = P1 + ((sqrt(5.0)-1.0)/2.0)*(P4-P1)

                dum2=((k_grid(ik) + Gamma_redistr*0.5d0)*(1d0+r*(1d0-tk))+lumpsum*0.5d0+Unemp_benefit-P2*(1d0+tc))/(1d0+mu)
                if(dum2<0.0001d0) then
                    V2=-999999999d0
                else
                    pnt1 = [dum2,exp_grid(ix,T)*(1d0-deltaexp)]
                    vnext = bilin_interp(k_grid,exp_grid(:,T+1),evs_ret(2,j,:,:,iam,ium,1,ifc),nk,nexp, pnt1)
                    !vnext = D_QD2VL(dum2,exp_grid(ix,T)*(1d0-deltaexp),k_grid,exp_grid(:,T+1),evs_ret(2,j,:,:,iam,ium,1,ifc))
                    V2=Uc(P2)
                    V2=V2+beta*OmegaActive(T)*vnext
                end if

                dum2=((k_grid(ik) + Gamma_redistr*0.5d0)*(1d0+r*(1d0-tk))+lumpsum*0.5d0+Unemp_benefit-P3*(1d0+tc))/(1d0+mu)
                if(dum2<0.0001d0) then
                    V3=-999999999d0
                else
                    pnt1 = [dum2,exp_grid(ix,T)*(1d0-deltaexp)]
                    vnext = bilin_interp(k_grid,exp_grid(:,T+1),evs_ret(2,j,:,:,iam,ium,1,ifc),nk,nexp, pnt1)
                    !vnext = D_QD2VL(dum2,exp_grid(ix,T)*(1d0-deltaexp),k_grid,exp_grid(:,T+1),evs_ret(2,j,:,:,iam,ium,1,ifc))
                    V3=Uc(P3)
                    V3=V3+beta*OmegaActive(T)*vnext
                end if

                if (V2 < V3) then
                    P1=P2
                else
                    P4=P3
                end if
                if((P4-P1)<1d-6) exit
            end do

            Vus=V2
            kus=dum2
            nus=0d0
            cus=P2

            if (ves >= vus) then
                vs(j,ik,ix,iam,ium,T,ifc)=ves
                cs(j,ik,ix,iam,ium,T,ifc)=ces
                ks(j,ik,ix,iam,ium,T,ifc)=kes
                ns(j,ik,ix,iam,ium,T,ifc)=nes
            else
                vs(j,ik,ix,iam,ium,T,ifc)=vus
                cs(j,ik,ix,iam,ium,T,ifc)=cus
                ks(j,ik,ix,iam,ium,T,ifc)=kus
                ns(j,ik,ix,iam,ium,T,ifc)=nus
            end if
        end do
    end do


end subroutine Solvelastactive
