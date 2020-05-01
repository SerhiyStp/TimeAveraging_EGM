module SolveRetired_EGM
use Model_Parameters
    implicit none

    real(8), dimension(:,:,:,:,:,:,:,:,:,:,:), allocatable, target :: evv_ret
    real(8), dimension (:,:,:,:,:,:,:), allocatable, target :: evvs_ret

    real(8) :: z_a, z_b, z_c
    real(8) :: z_eta

    type dev
        real(8) :: f(nk)
        integer :: nc_ilo, nc_ihi
    contains
        procedure :: FindNonConvexRegion
    end type dev

    type(dev) :: dev_ret(2,2,nexp,nexp,na,nu,na,nu,nfc,nfc)
    type(dev) :: devs_ret(2,2,nexp,na,nu,nfc)

contains

    subroutine SolveRetired_All()

        call OptimizeRetired_LastPeriod()
        call Compute_EV_DEV(Tret)
        call OptimizeRetired(Tret-1)
        
    end subroutine SolveRetired_All

    function hrs_foc(h) result(fval)
        implicit none
        real(8), intent(in) :: h
        real(8) :: fval

        fval = z_a*h**z_eta - z_b*h**(-theta(2)) + z_c
        
    end function hrs_foc

    subroutine Initialize_Retired()
        integer :: j, ic, ia, iu, ix
        integer :: ia2, iu2, ix2

        allocate(evv_ret(2,2,nk,nexp,nexp,na,nu,na,nu,nfc,nfc))
        allocate(evvs_ret(2,2,nk,nexp,na,nu,nfc))

    end subroutine Initialize_Retired

    subroutine Compute_EV_DEV(i_t)
    use Utilities, only: LinInterp
    use PolicyFunctions, only: evs_ret, ev_ret, k_grid, vs_ret, v_ret
    use pyplot_module
        integer :: i_t
        integer :: j, ik, ix, ia, iu, ifc
        real(8) :: v0, v1, k0, k1
        integer :: i1, i2
        integer :: ixm, ixf, iam, iaf, ium, iuf
        integer :: ifcm, ifcf
        real(8) :: eps
        real(8), dimension(:), pointer :: ev_k_ptr
        type(pyplot) :: plt
        real(8) :: dev_plot(nk)
        integer :: ix_plt, ia_plt, iu_plt, ifc_plt
        integer :: ixm_plt, iam_plt, ium_plt, ixf_plt, iaf_plt, iuf_plt
        integer :: ifcm_plt, ifcf_plt
        real(8) :: evvs
        integer :: iu2, iu3
        real(8) :: df

        !Singles
        do j = 1, 2
            do iu = 1, nu
                evvs_ret(2,j,:,:,:,iu,:) = 0.0d0
                do iu2 = 1, nu
                    evvs_ret(2,j,:,:,:,iu,:) = evvs_ret(2,j,:,:,:,iu,:) &
                        + trans_u(j,iu,iu2)*vs_ret(2,j,:,:,:,iu2,i_t,:)
                end do
            end do
        end do

        evvs_ret(1,:,:,:,:,:,:) = vs_ret(1,:,:,:,:,:,i_t,:)

        eps = 1.0d-5

        do ix = 1, nexp
        do ia = 1, na
        do iu = 1, nu
        do ifc = 1, nfc
        do i1 = 1, 2
        do i2 = 1, 2
            ev_k_ptr => &
                evvs_ret(i1,i2,:,ix,ia,iu,ifc)
                !evs_ret(i1,i2,:,ix,ia,iu,i_t,ifc)
            do ik = 1, nk
                k0 = k_grid(ik) - eps
                k1 = k_grid(ik) + eps
                v0 = LinInterp(k0,k_grid,ev_k_ptr,nk)
                v1 = LinInterp(k1,k_grid,ev_k_ptr,nk)
                df = (v1 - v0)/(2d0*eps)
                devs_ret(i1,i2,ix,ia,iu,ifc)%f(ik) = df
            end do
        end do
        end do
        end do
        end do
        end do
        end do

        !Married (only man or woman works)
        !(2,2,nk,nexp,nexp,na,nu,na,nu,nfc,nfc)
        do iu = 1, nu
            evv_ret(2,1,:,:,:,:,:,:,iu,:,:) = 0.0d0
            do iu2 = 1, nu
                evv_ret(2,1,:,:,:,:,:,:,iu,:,:) = &
                    evv_ret(2,1,:,:,:,:,:,:,iu,:,:) &
                    + trans_u(2,iu,iu2)*v_ret(2,1,:,:,:,:,:,:,iu2,i_t,:,:)
            end do
            evv_ret(1,2,:,:,:,:,iu,:,:,:,:) = 0.0d0
            do iu2 = 1, nu
                evv_ret(1,2,:,:,:,:,iu,:,:,:,:) = &
                    evv_ret(1,2,:,:,:,:,iu,:,:,:,:) &
                    + trans_u(1,iu,iu2)*v_ret(1,2,:,:,:,:,iu,:,:,i_t,:,:)
            end do
        end do
        !Married (both man and woman work)
        do ium = 1, nu
            do iuf = 1, nu
                evv_ret(2,2,:,:,:,:,ium,:,iuf,:,:) = 0.0d0
                do iu2 = 1, nu
                    do iu3 = 1, nu
                        evv_ret(2,2,:,:,:,:,ium,:,iuf,:,:) &
                            = evv_ret(2,2,:,:,:,:,ium,:,iuf,:,:) &
                            + trans_u(1,ium,iu2)*trans_u(2,iuf,iu3) &
                            * v_ret(2,2,:,:,:,:,iu2,:,iu3,i_t,:,:)
                    end do
                end do
            end do
        end do
        !Married (both retired)
        evv_ret(1,1,:,:,:,:,:,:,:,:,:) = v_ret(1,1,:,:,:,:,:,:,:,i_t,:,:)

        do i1 = 1, 2
        do i2 = 1, 2
        do ixm = 1, nexp
        do ixf = 1, nexp
        do iam = 1, na
        do iaf = 1, na
        do ium = 1, nu
        do iuf = 1, nu
        do ifcm = 1, nfc
        do ifcf = 1, nfc
            ev_k_ptr => &
                evv_ret(i1,i2,:,ixm,ixf,iam,ium,iaf,iuf,ifcm,ifcf)
            do ik = 1, nk
                k0 = k_grid(ik) - eps
                k1 = k_grid(ik) + eps
                v0 = LinInterp(k0,k_grid,ev_k_ptr,nk)
                v1 = LinInterp(k1,k_grid,ev_k_ptr,nk)
                df = (v1 - v0)/(2d0*eps)
                dev_ret(i1,i2,ixm,ixf,iam,ium,iaf,iuf,ifcm,ifcf)%f(ik) = df
            end do
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        end do

        ix_plt = 3
        ia_plt = 3
        iu_plt = 3
        ifc_plt = 2
        dev_plot = devs_ret(1,1,ix_plt,ia_plt,iu_plt,ifc_plt)%f

        call plt%initialize(grid=.true.,xlabel='k',&
            title='MU',legend=.true.)
        call plt%add_plot(k_grid,dev_plot,label='MU (single)',&
            linestyle='b-o',markersize=5,linewidth=2)

        ixm_plt = 3
        ixf_plt = 2
        iam_plt = 2
        iaf_plt = 3
        ium_plt = 2
        iuf_plt = 2
        ifcm_plt = 3
        ifcf_plt = 2

        dev_plot = &
            dev_ret(1,1,ixm_plt,ixf_plt,iam_plt,ium_plt,iaf_plt,iuf_plt,&
                ifcm_plt,ifcf_plt)%f

        call plt%add_plot(k_grid,dev_plot,label='MU (married)',&
            linestyle='r--',markersize=5,linewidth=2)
        !call plt%savefig('dev3.png', pyfile='dev2.py')
        call plt%savefig('dev3.png')

    end subroutine Compute_EV_DEV

    subroutine OptimizeRetired_LastPeriod()
    use PolicyFunctions
    use Utilities, only: Uc
        integer :: ik
        real(8) :: cons
        real(8) :: rr, kk

        rr = ((1d0 + r)/OmegaRet2(Tret))*(1d0-tk)
        do ik = 1, nk
            kk = k_grid(ik)+Gamma_redistr
            k_ret(:,:,ik,:,:,:,:,:,:,Tret,:,:) = 0d0
            cons = (kk*rr+Psi_pension+lumpsum)/(1d0+tc)
            c_ret(:,:,ik,:,:,:,:,:,:,Tret,:,:) = cons
            v_ret(:,:,ik,:,:,:,:,:,:,Tret,:,:) = Uc(cons)

            kk = k_grid(ik)+Gamma_redistr*0.5d0
            ks_ret(:,:,ik,:,:,:,Tret,:) = 0.0d0
            cons = (kk*rr+(Psi_pension+lumpsum)*0.5d0)/(1d0+tc)
            cs_ret(:,:,ik,:,:,:,Tret,:) = cons
            ns_ret(:,:,ik,:,:,:,Tret,:) = 0.0d0
            vs_ret(:,:,ik,:,:,:,Tret,:) = Uc(cons)
        end do
        
    end subroutine OptimizeRetired_LastPeriod

    subroutine OptimizeRetired(i_t) 
        integer :: i_t
        integer :: i1, i2, ixm, ixf, iam, ium, iaf, iuf, ifm, iff
        integer :: j, i, ix, ia, iu, ifc
        integer :: ik
        real(8) :: d_ev
        real(8) :: ap_a_map_s(nk, nk, 2)
        real(8) :: cc, hh
        real(8) :: hlo, hhi, aerr, rerr
        real(8) :: zeroin
        external :: zeroin

        do i1 = 1, 2
        do i2 = 1, 2
        do ixm = 1, nexp
        do ixf = 1, nexp
        do iam = 1, na
        do ium = 1, nu
        do iaf = 1, na
        do iuf = 1, nu
        do ifm = 1, nfc
        do iff = 2, nfc
            call dev_ret(i1,i2,ixm,ixf,&
                iam,ium,iaf,iuf,ifm,iff)%FindNonConvexRegion() 
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        end do

        do j = 1, 2
        do i = 1, 2
        do ix = 1, nexp
        do ia = 1, na
        do iu = 1, nu
        do ifc = 1, nfc
            call devs_ret(j,i,ix,ia,iu,ifc)%FindNonConvexRegion()
            do ik = 1, nk
                d_ev = devs_ret(j,i,ix,ia,iu,ifc)%f(ik)
                cc = 1d0/(beta*OmegaRet(i_t+1)*d_ev*(1d0+tc))
                hlo = 1d-6
                hhi = 100d0
                aerr = 1d-5
                rerr = 1d-5
                hh = zeroin(hrs_foc, hlo, hhi, aerr, rerr)
            end do 
        end do
        end do
        end do
        end do
        end do
        end do

    end subroutine OptimizeRetired

    subroutine FindNonConvexRegion(this)
        class(dev) :: this
        integer, dimension(nk) :: INCR_IDX_H, INCR_IDX_L, ILO_IDX, IHI_IDX
        real(8) :: vmax, vmin
        logical :: ldv_ncr, ldv_gv, ldv_sv
        integer :: ia
        integer :: ilo, ihi

        do ia = 2, nk
            ldv_ncr          = this%f(ia) .gt. this%f(ia-1)
            INCR_IDX_H(ia-1) = logicfunc(ldv_ncr)
            INCR_IDX_L(ia)   = logicfunc(ldv_ncr)
        end do

        if( maxval(INCR_IDX_H) .ne. 0) then
            vmin = minval(this%f, mask = INCR_IDX_H .eq. 1)
            vmax = maxval(this%f, mask = INCR_IDX_L .eq. 1)
            do ia = 1, nk
                ldv_gv      = this%f(ia) .gt. vmax
                ILO_IDX(ia) = logicfunc(ldv_gv)*ia
            end do
            if( maxval(ILO_IDX) .ne. 0) then
                ilo = maxval(ILO_IDX, mask = ILO_IDX .gt. 0) !sum(ILO_IDX)-1
            else
                ilo = 1
            end if
            do ia = 1, nk
                ldv_sv      = this%f(ia) .lt. vmin
                IHI_IDX(ia) = logicfunc(ldv_sv)*ia
            end do
            if( maxval(IHI_IDX) .ne. 0) then
                ihi = minval(IHI_IDX, mask = IHI_IDX .gt. 0) !N-sum(IHI_IDX)+1
            else
                ihi = 1
            end if
        else
            ilo = nk + 1
            ihi = 0
        end if
        this%nc_ilo = ilo
        this%nc_ihi = ihi

    end subroutine FindNonConvexRegion

    elemental pure double precision function logic2dbl(a)
        logical, intent(in) :: a
        if (a) then
            logic2dbl = 1.d0
        else
            logic2dbl = 0.d0
        end if
    end function logic2dbl
    
    elemental pure double precision function logicfunc(a)
        logical, intent(in) :: a
        if (a) then
            logicfunc = 1
        else
            logicfunc = 0
        end if
    end function logicfunc

end module SolveRetired_EGM
