module SolveRetired_EGM
use Model_Parameters
    implicit none
    real(8), dimension(:,:,:,:,:,:,:,:,:,:,:), allocatable, target :: dev_ret
    real(8), dimension (:,:,:,:,:,:,:), allocatable, target :: devs_ret
    real(8), dimension(:,:,:,:,:,:,:,:,:,:,:), allocatable, target :: evv_ret
    real(8), dimension (:,:,:,:,:,:,:), allocatable, target :: evvs_ret

contains

    subroutine Initialize_Retired()

        allocate(dev_ret(2,2,nk,nexp,nexp,na,nu,na,nu,nfc,nfc))
        allocate(devs_ret(2,2,nk,nexp,na,nu,nfc))
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
                devs_ret(i1,i2,ik,ix,ia,iu,ifc) = &
                    (v1 - v0)/(2d0*eps)
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
                !ev_ret(i1,i2,:,ixm,ixf,iam,ium,iaf,iuf,i_t,ifcm,ifcf)
            do ik = 1, nk
                k0 = k_grid(ik) - eps
                k1 = k_grid(ik) + eps
                v0 = LinInterp(k0,k_grid,ev_k_ptr,nk)
                v1 = LinInterp(k1,k_grid,ev_k_ptr,nk)
                dev_ret(i1,i2,ik,ixm,ixf,iam,ium,iaf,iuf,ifcm,ifcf) = &
                    (v1 - v0)/(2d0*eps)
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
        dev_plot = devs_ret(1,1,:,ix_plt,ia_plt,iu_plt,ifc_plt)

        call plt%initialize(grid=.true.,xlabel='k',&
            title='MU',legend=.true.)
        call plt%add_plot(k_grid,dev_plot,label='MU (single)',linestyle='b-o',markersize=5,linewidth=2)

        ixm_plt = 3
        ixf_plt = 2
        iam_plt = 2
        iaf_plt = 3
        ium_plt = 2
        iuf_plt = 2
        ifcm_plt = 3
        ifcf_plt = 2
        dev_plot = dev_ret(1,1,:,ixm_plt,ixf_plt,iam_plt,ium_plt,iaf_plt,iuf_plt,ifcm_plt,ifcf_plt)

        call plt%add_plot(k_grid,dev_plot,label='MU (married)',linestyle='r--',markersize=5,linewidth=2)
        call plt%savefig('dev2.png', pyfile='dev.py')

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

    subroutine SolveRetired_Mode1(i_t)
        integer :: i_t
        
    end subroutine SolveRetired_Mode1

end module SolveRetired_EGM
