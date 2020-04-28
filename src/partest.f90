subroutine partest(counter)
    !This subroutine computes optimal policies at age 64
    use Model_Parameters
    use PolicyFunctions
    use glob0
    use Utilities
    use MyInterpolation

    implicit none

    integer, INTENT(IN) :: counter
    integer :: ix,ixm,ixd,iam,ium,iaf,iuf,iu2,iu3,j,ik2,j2,ifc,ifcm,ik
    real(8), dimension (:,:,:,:,:,:), allocatable :: ce,cu,ke,ku,nem,nef,num,nuf,ve,vu
    real(8), dimension (:,:,:,:,:), allocatable :: ces,cus,kes,kus,nes,nus,ves,vus
    real(8), dimension (:), allocatable :: Expdum
    integer :: NEQ=0, IERSVR=0, IPACT=0, ISACT=0
    real(8) :: c2, MU2, d1, d2, vp(nu),dum3,dum4,dum5,dum6,y
    real(8) :: ACC=0.0001d0,ERREL=0.0001d0
    real(8) :: P1,P2,P3,P4,V2,V3,dum2
    real(8) :: vnext

    !Assigning the grid points
    dum3=((counter*1d0)/(nu*na*1d0))-0.00001d0
    ik=int(dum3)+1
    dum3=(((counter-(ik-1)*nu*na)*1d0)/(na*1d0))-0.00001d0
    ium=int(dum3)+1
    iam=counter-(ik-1)*na*nu-(ium-1)*nu  

    j=2
    evm(j,ik,:,iam,ium,T-it,:)=0d0
    do ifc=1,nfc
    do ix = 1, nexp
    do iu2=1,nu
    do ik2=1,nk
    do ixm = 1, nexp
    do iaf=1,na
    do iuf=1,nu
    do ifcm=1,nfc
        dum=k_grid(ik)+k_grid(ik2)
        if(dum<k_grid(nk)-0.001d0) then
            !vnext = D_QDVAL(dum, k_grid, v(:,ix,ixm,iaf,iuf,iam,iu2,T-it,ifc,ifcm))
            vnext = S_QDVAL(dum, k_grid, v(:,ix,ixm,iaf,iuf,iam,iu2,T-it,ifc,ifcm), nk)
            evm(j,ik,ix,iam,ium,T-it,ifc) = &
                evm(j,ik,ix,iam,ium,T-it,ifc) + &
                trans_u(2,ium,iu2)*ability_prob(iam,iaf)*mpartner(ik2,ixm,iaf,iuf,T-it,ifcm)*vnext
        else
            evm(j,ik,ix,iam,ium,T-it,ifc) = &
                evm(j,ik,ix,iam,ium,T-it,ifc) + &
                trans_u(2,ium,iu2)*ability_prob(iam,iaf)*mpartner(ik2,ixm,iaf,iuf,T-it,ifcm)*&
                LinInterp(dum,k_grid,v(:,ix,ixm,iaf,iuf,iam,iu2,T-it,ifc,ifcm),nk)
        end if
    end do
    end do
    end do
    end do
    end do
    end do
    end do
    end do

end subroutine partest
    
subroutine partest2(counter)
    !This subroutine computes optimal policies at age 64
    use Model_Parameters
    use PolicyFunctions
    use glob0
    use Utilities
    use MyInterpolation

    implicit none

    integer, INTENT(IN) :: counter
    integer :: ix,ixm,ixd,iam,ium,iaf,iuf,iu2,iu3,j,ik2,j2,ifc,ifcm,ik
    real(8), dimension (:,:,:,:,:,:), allocatable :: ce,cu,ke,ku,nem,nef,num,nuf,ve,vu
    real(8), dimension (:,:,:,:,:), allocatable :: ces,cus,kes,kus,nes,nus,ves,vus
    real(8), dimension (:), allocatable :: Expdum
    integer :: NEQ=0, IERSVR=0, IPACT=0, ISACT=0
    real(8) :: c2, MU2, d1, d2, vp(nu),dum3,dum4,dum5,dum6,y
    real(8) :: ACC=0.0001d0,ERREL=0.0001d0
    real(8) :: P1,P2,P3,P4,V2,V3,dum2
    real(8) :: vnext

    !Assigning the grid points
    dum3=((counter*1d0)/(nu*na*1d0))-0.00001d0
    ik=int(dum3)+1
    dum3=(((counter-(ik-1)*nu*na)*1d0)/(na*1d0))-0.00001d0
    ium=int(dum3)+1
    iam=counter-(ik-1)*na*nu-(ium-1)*nu

    j=1
    evm(j,ik,:,iam,ium,T-it,:)=0d0
    do ixm=1,nexp
        do ifcm=1,nfc
            do iu2=1,nu
                do ik2=1,nk
                    do ix=1,nexp
                        do iaf=1,na
                            do iuf=1,nu
                                do ifc=1,nfc
                                    dum=k_grid(ik)+k_grid(ik2)
                                    if(dum<k_grid(nk)-0.001d0) then
                                        vnext = S_QDVAL(dum, k_grid, v(:,ix,ixm,iam,iu2,iaf,iuf,T-it,ifc,ifcm), nk)
                                        evm(j,ik,ixm,iam,ium,T-it,ifcm) = &
                                            evm(j,ik,ixm,iam,ium,T-it,ifcm) + &
                                            trans_u(1,ium,iu2)*ability_prob(iam,iaf)*fpartner(ik2,ix,iaf,iuf,T-it,ifc)*vnext
                                    else
                                        evm(j,ik,ixm,iam,ium,T-it,ifcm) = &
                                            evm(j,ik,ixm,iam,ium,T-it,ifcm) + &
                                            trans_u(1,ium,iu2)*ability_prob(iam,iaf)*fpartner(ik2,ix,iaf,iuf,T-it,ifc)*&
                                            LinInterp(dum,k_grid,v(:,ix,ixm,iam,iu2,iaf,iuf,T-it,ifc,ifcm),nk)
                                    end if
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do

end subroutine partest2
    
subroutine partest3(counter)
    !This subroutine computes optimal policies at age 64
    use Model_Parameters
    use PolicyFunctions
    use glob0
    use Utilities

    implicit none

    integer, INTENT(IN) :: counter
    integer :: ix,ixd,iam,iaf,iuf,iu2,iu3,j,ik2,j2,ifc,ium
    real(8), dimension (:,:,:,:,:,:), allocatable :: ce,cu,ke,ku,nem,nef,num,nuf,ve,vu
    real(8), dimension (:,:,:,:,:), allocatable :: ces,cus,kes,kus,nes,nus,ves,vus
    real(8), dimension (:), allocatable :: Expdum
    integer :: NEQ=0, IERSVR=0, IPACT=0, ISACT=0
    real(8) :: c2, MU2, d1, d2, vp(nu),dum3,dum4,dum5,dum6,y
    real(8) :: ACC=0.0001d0,ERREL=0.0001d0
    real(8) :: P1,P2,P3,P4,V2,V3,dum2

    !Assigning the grid points
    dum3=((counter*1d0)/(na*nfc*1d0))-0.00001d0
    ium=int(dum3)+1
    dum3=(((counter-(ium-1)*nfc*na)*1d0)/(na*1d0))-0.00001d0
    ifc=int(dum3)+1
    iam=counter-(ium-1)*nfc*na-(ifc-1)*na


    !do iaf=1,na
    !    do iuf=1,nu
    !        !call d_csint(k_grid, ev(:,ix,iam,ium,iaf,iuf,T-it,j2), BREAK, ev_spln_coefs(:,:,ix,iam,ium,iaf,iuf,T-it,j2))
    !        CALL D_BS2IN(k_grid, exp_grid(:,T-it), ev(:,:,iam,ium,iaf,iuf,T-it,ifc), KORDER, EXPORDER, K_KNOT, EXP_KNOT(:,T-it), ev_spln_coefs(:,:,iam,ium,iaf,iuf,T-it,ifc) , nk)
    !        CALL D_BS2IN(k_grid, exp_grid(:,T-it), v(:,:,iam,ium,iaf,iuf,T-it,ifc), KORDER, EXPORDER, K_KNOT, EXP_KNOT(:,T-it), v_spln_coefs(:,:,iam,ium,iaf,iuf,T-it,ifc) , nk)
    !    end do
    !end do


end subroutine partest3

subroutine partest4(counter)
    !This subroutine computes optimal policies at age 64
    use Model_Parameters
    use PolicyFunctions
    use glob0
    use Utilities

    implicit none

    integer, INTENT(IN) :: counter
    integer :: ix,ixd,iam,iaf,iuf,iu2,iu3,j,ik2,j2,ifc,ium
    real(8), dimension (:,:,:,:,:,:), allocatable :: ce,cu,ke,ku,nem,nef,num,nuf,ve,vu
    real(8), dimension (:,:,:,:,:), allocatable :: ces,cus,kes,kus,nes,nus,ves,vus
    real(8), dimension (:), allocatable :: Expdum
    integer :: NEQ=0, IERSVR=0, IPACT=0, ISACT=0
    real(8) :: c2, MU2, d1, d2, vp(nu),dum3,dum4,dum5,dum6,y
    real(8) :: ACC=0.0001d0,ERREL=0.0001d0
    real(8) :: P1,P2,P3,P4,V2,V3,dum2

    !Assigning the grid points
    dum3=((counter*1d0)/(na*nfc*1d0))-0.00001d0
    ium=int(dum3)+1
    dum3=(((counter-(ium-1)*nfc*na)*1d0)/(na*1d0))-0.00001d0
    ifc=int(dum3)+1
    iam=counter-(ium-1)*nfc*na-(ifc-1)*na

    do ix = 1, nexp
        do iaf=1,na
            do iuf=1,nu
                !call d_csint(k_grid, v(:,ix,iam,ium,iaf,iuf,T-it,j2), BREAK, v_spln_coefs(:,:,ix,iam,ium,iaf,iuf,T-it,j2))
                !call ppp_csint(k_grid, v(:,ix,iam,ium,iaf,iuf,T-it,ifc), v_spln_coefs(:,:,ix,iam,ium,iaf,iuf,T-it,ifc), nk)
            end do
        end do
    end do

end subroutine partest4
    
subroutine partest7(counter)
    !This subroutine computes optimal policies at age 64
    use Model_Parameters
    use PolicyFunctions
    use glob0
    use Utilities
    use MyInterpolation

    implicit none

    integer, INTENT(IN) :: counter
    integer :: ix,ixm,ixd,iam,ium,iaf,iuf,iu2,iu3,j,ik2,j2,ifc,ifcm,ik
    real(8), dimension (:,:,:,:,:,:), allocatable :: ce,cu,ke,ku,nem,nef,num,nuf,ve,vu
    real(8), dimension (:,:,:,:,:), allocatable :: ces,cus,kes,kus,nes,nus,ves,vus
    real(8), dimension (:), allocatable :: Expdum
    integer :: NEQ=0, IERSVR=0, IPACT=0, ISACT=0
    real(8) :: c2, MU2, d1, d2, vp(nu),dum3,dum4,dum5,dum6,y
    real(8) :: ACC=0.0001d0,ERREL=0.0001d0
    real(8) :: P1,P2,P3,P4,V2,V3,dum2
    real(8) :: vnext

    !Assigning the grid points
    dum3=((counter*1d0)/(nu*na*1d0))-0.00001d0
    ik=int(dum3)+1
    dum3=(((counter-(ik-1)*nu*na)*1d0)/(na*1d0))-0.00001d0
    ium=int(dum3)+1
    iam=counter-(ik-1)*na*nu-(ium-1)*nu

    if(T-it>1) then
        ev(ik,:,:,iam,ium,:,:,T-it,:,:)=0d0
        dum3=K_grid(ik)/2d0
        do ifc=1,nfc
        do ifcm=1,nfc
        do ix = 1, nexp
        do ixm = 1, nexp
        do iuf = 1, nu
        do iaf = 1, na
            do iu2=1,nu
            do iu3=1,nu
                ev(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm) = &
                    ev(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm) &
                    + (1d0-Probd(T-it-1))*trans_u(1,ium,iu2)*trans_u(2,iuf,iu3)*V(ik,ix,ixm,iam,iu2,iaf,iu3,T-it,ifc,ifcm)
            end do
            end do
            do iu3=1,nu
                vnext = S_QDVAL(dum3,k_grid,vs(2,:,ix,iaf,iu3,T-it,ifc), nk)
                ev(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm) = &
                    ev(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm) + &
                    Probd(T-it-1)*trans_u(2,iuf,iu3)*0.5d0*vnext
            end do
            do iu2=1,nu
                vnext = S_QDVAL(dum3,k_grid,vs(1,:,ixm,iam,iu2,T-it,ifcm), nk)
                ev(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm) = &
                    ev(ik,ix,ixm,iam,ium,iaf,iuf,T-it,ifc,ifcm) + &
                    Probd(T-it-1)*trans_u(1,ium,iu2)*0.5d0*vnext
            end do
        end do
        end do
        end do
        end do
        end do
        end do
    end if

end subroutine partest7
    
subroutine partest8(counter)
    !This subroutine computes optimal policies at age 64
    use Model_Parameters
    use PolicyFunctions
    use glob0
    use Utilities

    implicit none

    integer, INTENT(IN) :: counter
    integer :: ix,ixd,iam,ium,iaf,iuf,iu2,iu3,j,ik2,j2,ifc,ik
    real(8), dimension (:,:,:,:,:,:), allocatable :: ce,cu,ke,ku,nem,nef,num,nuf,ve,vu
    real(8), dimension (:,:,:,:,:), allocatable :: ces,cus,kes,kus,nes,nus,ves,vus
    real(8), dimension (:), allocatable :: Expdum
    integer :: NEQ=0, IERSVR=0, IPACT=0, ISACT=0
    real(8) :: c2, MU2, d1, d2, vp(nu),dum3,dum4,dum5,dum6,y
    real(8) :: ACC=0.0001d0,ERREL=0.0001d0
    real(8) :: P1,P2,P3,P4,V2,V3,dum2
    real(8) :: vnext

    !Assigning the grid points
    dum3=((counter*1d0)/(nu*na*1d0))-0.00001d0
    ik=int(dum3)+1
    dum3=(((counter-(ik-1)*nu*na)*1d0)/(na*1d0))-0.00001d0
    ium=int(dum3)+1
    iam=counter-(ik-1)*na*nu-(ium-1)*nu

    evs(:,ik,:,iam,ium,T-it,:)=0d0

    do j=1,2
        do ifc=1,nfc
            do ix = 1, nexp
                do iu2=1,nu
                    evs(j,ik,ix,iam,ium,T-it,ifc)=evs(j,ik,ix,iam,ium,T-it,ifc)+trans_u(2,ium,iu2)*Vs(j,ik,ix,iam,iu2,T-it,ifc)
                end do
            end do
        end do
    end do        

end subroutine partest8
    
subroutine partest9(counter)
    !This subroutine computes optimal policies at age 64
    use Model_Parameters
    use PolicyFunctions
    use glob0
    use Utilities

    implicit none

    integer, INTENT(IN) :: counter
    integer :: ix,ixm,ixd,iam,ium,iaf,iuf,iu2,iu3,j,ik2,j2,ifc,ifcm,ik
    real(8), dimension (:,:,:,:,:,:), allocatable :: ce,cu,ke,ku,nem,nef,num,nuf,ve,vu
    real(8), dimension (:,:,:,:,:), allocatable :: ces,cus,kes,kus,nes,nus,ves,vus
    real(8), dimension (:), allocatable :: Expdum
    integer :: NEQ=0, IERSVR=0, IPACT=0, ISACT=0
    real(8) :: c2, MU2, d1, d2, vp(nu),dum3,dum4,dum5,dum6,y
    real(8) :: ACC=0.0001d0,ERREL=0.0001d0
    real(8) :: P1,P2,P3,P4,V2,V3,dum2
    real(8) :: vnext

    !Assigning the grid points
    dum3=((counter*1d0)/(nu*na*1d0))-0.00001d0
    ik=int(dum3)+1
    dum3=(((counter-(ik-1)*nu*na)*1d0)/(na*1d0))-0.00001d0
    ium=int(dum3)+1
    iam=counter-(ik-1)*na*nu-(ium-1)*nu

    !Singles
    evs_ret(2,:,ik,:,iam,ium,Tret-it+1,:)=0d0
    do j=1,2
    do ifc=1,nfc
    do ix = 1, nexp
    do iu2=1,nu
        evs_ret(2,j,ik,ix,iam,ium,Tret-it+1,ifc) &
            =evs_ret(2,j,ik,ix,iam,ium,Tret-it+1,ifc)&
            +trans_u(j,ium,iu2)*vs_ret(2,j,ik,ix,iam,iu2,Tret-it+1,ifc)
    end do
    end do
    end do
    end do


    !Woman works
    ev_ret(2,1,ik,:,:,:,:,iam,ium,Tret-it+1,:,:)=0d0
    do ifc=1,nfc
    do ix = 1, nexp
    do iu2=1,nu
        ev_ret(2,1,ik,ix,:,:,:,iam,ium,Tret-it+1,ifc,:) = &
            ev_ret(2,1,ik,ix,:,:,:,iam,ium,Tret-it+1,ifc,:) &
            +trans_u(2,ium,iu2)*v_ret(2,1,ik,ix,:,:,:,iam,iu2,Tret-it+1,ifc,:)
    end do
    end do
    end do

    !Man works
    ev_ret(1,2,ik,:,:,iam,ium,:,:,Tret-it+1,:,:)=0d0
    do ifc=1,nfc
    do ix = 1, nexp
    do iu2=1,nu
        ev_ret(1,2,ik,:,ix,iam,ium,:,:,Tret-it+1,:,ifc) = &
            ev_ret(1,2,ik,:,ix,iam,ium,:,:,Tret-it+1,:,ifc) &
            +trans_u(1,ium,iu2)*v_ret(1,2,ik,:,ix,iam,iu2,:,:,Tret-it+1,:,ifc)
    end do
    end do
    end do

    !Both spouses work
    ev_ret(2,2,ik,:,:,iam,ium,:,:,Tret-it+1,:,:)=0d0
    do ifc=1,nfc
    do ifcm=1,nfc
    do ix = 1, nexp
    do ixm = 1, nexp
    do iuf = 1, nu
    do iaf = 1, na
    do iu2=1,nu
    do iu3=1,nu
        ev_ret(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm) &
            =ev_ret(2,2,ik,ix,ixm,iam,ium,iaf,iuf,Tret-it+1,ifc,ifcm) &
            +trans_u(1,ium,iu2)*trans_u(2,iuf,iu3)*v_ret(2,2,ik,ix,ixm,iam,iu2,iaf,iu3,Tret-it+1,ifc,ifcm)
    end do
    end do
    end do
    end do
    end do
    end do
    end do
    end do

end subroutine partest9
