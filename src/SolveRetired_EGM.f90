module SolveRetired_EGM
use Model_Parameters
    implicit none

    !real(8), dimension(:,:,:,:,:,:,:,:,:,:,:), allocatable :: dev_ret
    !real(8), dimension (:,:,:,:,:,:,:), allocatable :: devs_ret
    real(8), dimension(:,:,:,:,:,:,:,:,:,:,:), allocatable, target :: evv_ret
    real(8), dimension (:,:,:,:,:,:,:), allocatable, target :: evvs_ret

    type dev
        real(8) :: f(nk)
        integer :: ilo, ihi
    end type dev

    type(dev) :: dev_ret(2,2,nexp,nexp,na,nu,na,nu,nfc,nfc)
    type(dev) :: devs_ret(2,2,nexp,na,nu,nfc)

contains

    subroutine Initialize_Retired()

        !allocate(dev_ret(2,2,nk,nexp,nexp,na,nu,na,nu,nfc,nfc))
        !allocate(devs_ret(2,2,nk,nexp,na,nu,nfc))
        allocate(evv_ret(2,2,nk,nexp,nexp,na,nu,na,nu,nfc,nfc))
        allocate(evvs_ret(2,2,nk,nexp,na,nu,nfc))

    end subroutine Initialize_Retired

    subroutine Compute_EV_DEV(i_t)
    use Utilities, only: LinInterp
    use PolicyFunctions, only: evs_ret, ev_ret, k_grid, vs_ret
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
        real(8) :: df
        integer :: iu2

        do j = 1, 2
            do iu = 1, nu
                evvs_ret(2,j,:,:,:,iu,:) = 0.0d0
                do iu2 = 1, nu
                    evvs_ret(2,j,:,:,:,iu,:) = evvs_ret(2,j,:,:,:,iu,:) + trans_u(j,iu,iu2)*vs_ret(2,j,:,:,:,iu2,i_t,:)
                    !evvs_ret(2,j,:,:,:,iu,:) = vs_ret(2,j,:,:,:,iu2,i_t,:)
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
                evs_ret(i1,i2,:,ix,ia,iu,i_t,ifc)
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
                ev_ret(i1,i2,:,ixm,ixf,iam,ium,iaf,iuf,i_t,ifcm,ifcf)
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

        ix_plt = 1
        ia_plt = 1
        iu_plt = 1
        ifc_plt = 1
        dev_plot = devs_ret(1,1,ix_plt,ia_plt,iu_plt,ifc_plt)%f

        call plt%initialize(grid=.true.,xlabel='k',&
            title='MU',legend=.true.)
        call plt%add_plot(k_grid,dev_plot,label='MU (single)',&
            linestyle='b-o',markersize=5,linewidth=2)

        ixm_plt = 1
        ixf_plt = 1
        iam_plt = 1
        iaf_plt = 1
        ium_plt = 1
        iuf_plt = 1
        ifcm_plt = 1
        ifcf_plt = 1
        dev_plot = &
            dev_ret(1,1,ixm_plt,ixf_plt,iam_plt,ium_plt,iaf_plt,iuf_plt,&
                ifcm_plt,ifcf_plt)%f

        call plt%add_plot(k_grid,dev_plot,label='MU (married)',&
            linestyle='r--',markersize=5,linewidth=2)
        call plt%savefig('dev.png', pyfile='dev.py')

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
