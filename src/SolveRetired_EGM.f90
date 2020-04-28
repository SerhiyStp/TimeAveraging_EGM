module SolveRetired_EGM
use Model_Parameters
    implicit none
    real(8), dimension(:,:,:,:,:,:,:,:,:,:,:,:), allocatable, target :: dev_ret
    real(8), dimension (:,:,:,:,:,:,:,:), allocatable, target :: devs_ret

contains

    subroutine Initialize_Retired()

        allocate(dev_ret(2,2,nk,nexp,nexp,na,nu,na,nu,Tret,nfc,nfc))
        allocate(devs_ret(2,2,nk,nexp,na,nu,Tret,nfc))

    end subroutine Initialize_Retired

    subroutine Compute_DEV(i_t)
    use Utilities, only: LinInterp
    use PolicyFunctions, only: evs_ret, ev_ret, k_grid
        integer :: i_t
        integer :: j, ik, ix, ia, iu, ifc
        real(8) :: v0, v1, k0, k1
        integer :: i1, i2
        integer :: ixm, ixf, iam, iaf, ium, iuf
        integer :: ifcm, ifcf
        real(8) :: eps
        real(8), dimension(:), pointer :: ev_k_ptr

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
                devs_ret(i1,i2,ik,ix,ia,iu,i_t,ifc ) = &
                    (v1 - v0)/(2d0*eps)
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
                devs_ret(i1,i2,ik,ixm,ixf,iam,ium,iaf,iuf,i_t,ifcm,infcf) = &
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

    end subroutine Compute_DEV

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