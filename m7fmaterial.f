C
C ADJUSTING PARAMETERS TO FIT GIVEN MATERIAL TEST DATA
C
C Using the free parameters to fit given data may sometimes lead to less
C than perfect fits of test data. In that case, it is necessary to
C adjust some of the fixed, or hard-to-change, parameters c_i. This is
C quite challenging. Improving the fit of one type of test by changing
C one parameter, c_j, may spoil the fit of another type of test. It may
C also necessitate changing several parameters c_j, which is a rather
C difficult nonlinear optimization problem. Attempts by several
C researchers to automatically adjust these parameters using
C optimization software has been unsuccessful in the past because of the
C complexity of the response. To succeed requires considerable
C experience. In that case, you can obtain assistance from
C Prof. F.C. Caner ( <ferhun.caner@upc.edu>, Tel. +34659816715, Skype
C ferhun.caner). He has ample experience in the fitting of test data
C with M7 using the finite element method.
C
C The program is based on the following papers, which can be freely
C downloaded as 527.pdf, 528.pdf, 519.pdf and 547.pdf, respectively,
C from ZP Bazant's website:
C http://www.civil.northwestern.edu/people/bazant/PDFs/Papers/
C
C [1] Caner, F.C., and Ba\v zant, Z.P. (2013). ``Microplane model M7 for
C     plain concrete: I. formulation." {\em ASCE J. of Engrg. Mechanics}
C     139 (12), Dec., 1714--1723.
C
C [2] Caner, F.C., and Ba\v zant, Z.P. (2013). ``Microplane model M7 for
C     plain concrete: II. calibration and verification." {\em ASCE J. of
C     Engrg. Mechanics} 139 (12), Dec., 1724--1735.
C
C [3] Caner, F.C., Ba\v zant, Z.P., and Wendner, R. (2013). ``Microplane
C     model M7f for fiber reinforced concrete." {\em Engrg. Fracture
C     Mechanics} 105, 41--57.
C
C [4] Kirane, K., and Ba\v zant, Z.P. (2014). ``Microplane damage model
C     for fatigue of quasibrittle materials: Sub-critical crack growth,
C     lifetime and residual strength." {\em International Journal of
C     Fatigue} 70, 93--105.
C

C +--------------------------------------------------------------------+
C |        SUBROUTINE C_NORM_ELASTIC           |
C +--------------------------------------------------------------------+
      subroutine c_norm_elastic(ef,def,sn0,eps_N0_pos,eps_N0_neg,sv0,
     $snf_no_split,E_N,zeta,damageOld,damageNew,
     . !not modified
     . young, poisson, C0, c_18, c_19, c_20, c_21
     . !modified
     . )
      use, intrinsic :: iso_fortran_env, only: RK => real64
      implicit none
!================================================input
      !not modified
      double precision :: young, poisson, C0, c_18, 
     . c_19, c_20, c_21
      !modified

!====================================================
      integer, parameter ::  np=37, nvhm=5, nvhf=2, nvhi=np*nvhm+nvhf
      double precision, parameter :: PI=3.1415926535897932384626433832d0
      double precision, dimension(1:np) :: ef, def, sn0
      double precision, dimension(1:np) :: eps_N0_pos, eps_N0_neg
      double precision :: sv0,zeta
      double precision, dimension(1:np) :: snf_no_split
      double precision, dimension(1:np) :: E_N
      double precision, dimension(1:np) :: damageOld, damageNew

      double precision :: E_N_0
      double precision :: weight_factor, t0, t1, t2, t3, t4, tsum, fzeta
      integer :: isize, i
      E_N_0= young / (1.0d0 - 2.0d0*poisson)

      isize=size(ef)

      do i=1,isize
      if ((sn0(i) > 0.d0)) then
      weight_factor=4.0d1
      t0=1.d0
      t1=(weight_factor*zeta)**2.d0
      t2=(weight_factor*zeta)**4.d0
      t3=(weight_factor*zeta)**6.d0
      t4=(weight_factor*zeta)**8.d0
      tsum=t0+t1+t2
      fzeta=1.d0/tsum
      E_N(i) = E_N_0*fzeta*exp(-c_19*eps_N0_pos(i))
      if (sn0(i) > E_N_0*ef(i) .and. sn0(i)*def(i)<0.0d0) then
      E_N(i)=E_N_0
      end if
      else
      E_N(i) = E_N_0*(exp(-c_20*abs(eps_N0_neg(i))/
     $(1.d0+c_18*max(-sv0,0.d0)/E_N_0))+
     $c_21*max(-sv0,0.d0)/E_N_0)
      end if
      damageNew(i) = max(1.0d0 - E_N(i)/E_N_0, damageOld(i))
      end do

      snf_no_split = sn0*C0 + E_N*def

      return
      end subroutine c_norm_elastic

C +--------------------------------------------------------------------+
C |          SUBROUTINE C_NORM           |
C +--------------------------------------------------------------------+
      subroutine c_norm(ef,sv0,R_N,snb,
     . !not modified
     . young, poisson, k_1, c_2, 
     . c_3, c_4, C_R2, V_f,d_1, d_2, d_3, d_4, d_5,
     . d_6,
     . !modified
     . c_1)
      use, intrinsic :: iso_fortran_env, only: RK => real64
      implicit none
!================================================input
      !not modified
      double precision :: young, poisson, k_1, c_2, 
     . c_3, c_4, C_R2, V_f,d_1, d_2, d_3, d_4, d_5,
     . d_6,
     . !modified
     . c_1

!====================================================
      double precision, parameter :: PI=3.1415926535897932384626433832d0
      integer, parameter ::  np=37, nvhm=5, nvhf=2, nvhi=np*nvhm+nvhf

      double precision, dimension(1:np) :: ef
      double precision, dimension(1:np) :: snb
      double precision :: sv0, R_N
      double precision, dimension(1:np)  :: aux1
      double precision :: C_N0, RATE_FACTOR, fstar, beta_N
      double precision :: E_V, eb_N0, eb_NK, eb_N
      integer :: isize

      isize = size(ef)

      E_V = young/(1.0d0-2.0d0*poisson)
      eb_N0 = c_3*k_1
      eb_NK = c_4
!Beta_1
      c_1 = d_1*tanh(d_2*V_f - d_3)+d_4*exp(-max(-sv0-d_6,0.d0)/E_V*d_5)
      if (sv0.lt.0.0d0) then
      eb_N = eb_N0 - eb_NK / E_V * sv0
      else
      eb_N = eb_N0
      end if

      fstar = k_1*young*c_1
      beta_N = c_2*c_1*k_1

      C_N0=C_R2

      aux1 = max(ef-beta_N,0.0d0)/eb_N
      snb = fstar*exp(-aux1)

      RATE_FACTOR = C_N0*R_N
      snb = snb * (1.d0 + RATE_FACTOR)

      return
      end subroutine c_norm

! +--------------------------------------------------------------------+
! |           SUBROUTINE C_NORM_FIB         |
! +--------------------------------------------------------------------+
      subroutine c_norm_fib(ef,R_N,snb_fib,
     . !not Modified
     . young, k_1, C_R2, sig_fib_0,
     . p_1, p_2, xp_1, xp_2, xp_3
     . !modified
     . )
      use, intrinsic :: iso_fortran_env, only: RK => real64
      implicit none
!================================================input
      !not modified
      double precision :: young, k_1, C_R2, sig_fib_0,
     . p_1, p_2, xp_1, xp_2, xp_3
     . !modified
     .

!====================================================
      double precision, parameter :: PI=3.1415926535897932384626433832d0
      integer, parameter ::  np=37, nvhm=5, nvhf=2, nvhi=np*nvhm+nvhf

      double precision, dimension(1:np) :: ef
      double precision, dimension(1:np) :: snb_fib

      double precision :: R_N

      double precision :: C_N0, RATE_FACTOR

      double precision :: e_N_max, e_N_0, e_N_hat, eps_shift, pow
      integer :: isize, i

      isize = size(ef)

      eps_shift = xp_2
      pow = xp_3

      e_N_max = pow/p_2 + eps_shift

      e_N_0 = xp_1

      e_N_hat = sig_fib_0

      C_N0=C_R2

      do i=1, isize
      if ( ef(i)/k_1 < e_N_max ) then
      snb_fib(i) = young*p_1*k_1*
     $max(ef(i)/k_1 - eps_shift,0.0d0)**pow*
     $exp(-p_2*max(ef(i)/k_1 - eps_shift,0.0d0))
      else if ((ef(i)/k_1 >= e_N_max) .and. (ef(i)/k_1 < e_N_0)) then
      snb_fib(i) = young*p_1*k_1*(pow/p_2)**pow*exp(-1.0d0*pow)
      else if (ef(i)/k_1 >= e_N_0) then
      snb_fib(i) = young*p_1*k_1*
     $(ef(i)/k_1 - e_N_0  + e_N_max - eps_shift)**pow*
     $exp(-p_2*(ef(i)/k_1 - e_N_0 + e_N_max - eps_shift))

      end if
      end do

      RATE_FACTOR = C_N0*R_N
      snb_fib = snb_fib * (1.d0 + RATE_FACTOR)
      return
      end subroutine c_norm_fib

! +--------------------------------------------------------------------+
! |          SUBROUTINE C_DEV            |
! +--------------------------------------------------------------------+
      subroutine c_dev(ded, ed0, dev, ev0, C_d, R_D,
     $  sdneg, sdpos,
     . !not modified
     . young, poisson, k_1, c_20, C_R2, 
     . c_5_0, c_5_1, c_5_M4, c_6_0, c_6_1, c_6_2, c_7_0, 
     . c_7_1, c_7_M4, c_8_0, c_8_1, c_8_M4, c_9_0, c_9_1, 
     . c_9_2, c_9_M4, c_6_M4, p_5,p_6,V_f,cf_1,cf_2,cf_3, 
     . cf_1_1, cf_2_2,cf_3_3,f_c0,E_0,f_cp,
     . !modified
     . c_5,c_6,c_7, c_8, c_9)
      use, intrinsic :: iso_fortran_env, only: RK => real64
      implicit none
!================================================input
      !not modified
      double precision :: young, poisson, k_1, c_20, C_R2, 
     . c_5_0, c_5_1, c_5_M4, c_6_0, c_6_1, c_6_2, c_7_0, 
     . c_7_1, c_7_M4, c_8_0, c_8_1, c_8_M4, c_9_0, c_9_1, 
     . c_9_2, c_9_M4, c_6_M4, p_5,p_6,V_f,cf_1,cf_2,cf_3, 
     . cf_1_1, cf_2_2,cf_3_3,
     . !modified
     . c_5,c_6,c_7, c_8, c_9

!====================================================
      double precision, parameter :: PI=3.1415926535897932384626433832d0
      integer, parameter ::  np=37, nvhm=5, nvhf=2, nvhi=np*nvhm+nvhf

      double precision, dimension(1:np) :: ded, ed0
      double precision, dimension(1:np) :: C_d

      double precision, dimension(1:np) :: sdneg, sdpos

      double precision :: dev, ev0, R_D
      double precision :: par5, par6, Cd0, Ev
      double precision, dimension(1:np) :: edf, sdlower, sdupper
      double precision, dimension(1:np) :: sdupper_fib
      double precision :: RATE_FACTOR_C, RATE_FACTOR_T, C_DC0, C_DT0
      double precision :: evf, E_v

      double precision :: f_c0, E_0, f_cp, c_40, beta_15, beta_16
      double precision :: beta_17, beta_18,
     $beta_19
      integer :: isize

      isize = size(ded)

      evf = ev0 + dev

      Ev=young/(1.0d0 - 2.0d0*poisson)

      c_40=1.0d+0
      beta_15=c_5_1*exp(-c_40*(f_cp/young-f_c0/E_0))
      beta_16=c_8_1*exp(-c_40*(f_cp/young-f_c0/E_0))
      beta_17=c_7_1*exp(-c_40*(f_cp/young-f_c0/E_0))
      beta_18=c_6_0*exp(-c_40*(f_cp/young-f_c0/E_0))
      beta_19=c_9_0*exp(-c_40*(f_cp/young-f_c0/E_0))

      c_5 = (beta_15*tanh(c_5_0*max(-evf,0.0d0)/k_1))*
     $(1.0d0 + cf_1*tanh(cf_1_1*V_f))+cf_1*tanh(cf_1_1*V_f)+c_5_M4
!Beta_3
      c_8 = (beta_16*tanh(c_8_0*max(-evf,0.0d0)/k_1))*
     $(1.0d0 + cf_3*tanh(cf_3_3*V_f))+cf_3*tanh(cf_3_3*V_f)+c_8_M4
!Beta_2
      c_7 = beta_17*tanh(c_7_0*max(-evf,0.0d0)/k_1)+
     $cf_2*tanh(cf_2_2*V_f) + c_7_M4
!Question in the paper c_9 c_6 are zeros, why is it calculated like this here?
      if (beta_18*max(-evf/k_1-c_6_1,0.0d0)>=log(c_6_2)) then
      c_6 = c_6_M4*c_6_2
      else
      c_6 = c_6_M4*exp(beta_18*max(-evf/k_1-c_6_1,0.0d0))
      end if
      if (beta_19*max(-evf/k_1-c_9_1,0.0d0)>=log(c_9_2)) then
      c_9 = c_9_M4*c_9_2
      else
      c_9 = c_9_M4*exp(beta_19*max(-evf/k_1-c_9_1,0.0d0))
      end if

      par5 = k_1*young*c_8
      par6 = young*c_5*k_1
      Cd0=young/(1.d0+poisson)*(1.d0-4.d0*poisson)/(1.d0-2.d0*poisson)
      E_v = young/(1.d0-2.d0*poisson)

      C_DC0 = C_R2

      C_DT0 = C_R2

      C_d=Cd0

      edf = ed0+ded

      call dev_comp(edf, sdlower,
     . ! not modified
     . young, k_1, c_7, c_8, c_9
     . !modified
     . )

      call dev_tens(edf, sdupper,
     . !not modified
     . young,k_1,c_5,c_6,c_7,c_20
     . !modified
     . )
      call dev_tens_fib(edf, sdupper_fib,
     . !not modified
     . young, k_1,p_5,p_6
     . !modified
     . )

      RATE_FACTOR_C = C_DC0*R_D
      RATE_FACTOR_T = C_DT0*R_D

      sdlower = sdlower * (1.d0 + RATE_FACTOR_C)

      sdupper = sdupper * (1.d0 + RATE_FACTOR_T)

      sdneg=sdlower
      sdpos=sdupper

      return
      end subroutine c_dev

! +--------------------------------------------------------------------+
! |         SUBROUTINE DEV_COMP          |
! +--------------------------------------------------------------------+
      subroutine dev_comp(edf, sdlower,
     . ! not modified
     . young, k_1, c_7, c_8, c_9
     . !modified
     . )
      use, intrinsic :: iso_fortran_env, only: RK => real64
      implicit none
!================================================input
      !not modified
      double precision :: young, k_1, c_7, c_8, c_9
     . !modified
!====================================================
      integer, parameter ::  np=37
      double precision, dimension(1:np) :: edf
      double precision, dimension(1:np) :: sdlower
      double precision :: par1, par2, par5

      par1 = c_9*c_8*k_1
      par2 = c_7*k_1
      par5 = k_1*young*c_8

      sdlower=-par5/(1.d0+(max(-edf-par1,0.0d0)/par2)**2.0d0)

      return
      end subroutine dev_comp
C
C
C
C +--------------------------------------------------------------------+
C |          SUBROUTINE DEV_TENS         |
C +--------------------------------------------------------------------+

      subroutine dev_tens(edf, sdupper,
     . !not modified
     . young,k_1,c_5,c_6,c_7,c_20
     . !modified
     . )
      use, intrinsic :: iso_fortran_env, only: RK => real64
      implicit none
!================================================input
      !not modified
      double precision :: young,k_1,c_5,c_6,c_7,c_20
     . !modified
!====================================================
      integer, parameter ::  np=37
      double precision, dimension(1:np) :: edf
      double precision, dimension(1:np) :: sdupper
      double precision :: par2, par3, par4, par6
      par2 = c_7*k_1
      par3 = c_20
      par4 = c_6*c_5*k_1
      par6 = young*c_5*k_1
      sdupper=par6/(1.0d0+(max(edf-par4,0.0d0)/par2/par3)**2.0d0)

      return
      end subroutine dev_tens

! +--------------------------------------------------------------------+
! |          SUBROUTINE DEV_TENS_FIB        |
! +--------------------------------------------------------------------+
      subroutine dev_tens_fib(edf, sdupper_fib,
     . !not modified
     . young, k_1,p_5,p_6
     . !modified
     . )
      use, intrinsic :: iso_fortran_env, only: RK => real64
      implicit none
!================================================input
      !not modified
      double precision:: young, k_1,p_5,p_6
     . !modified
!====================================================
      integer, parameter ::  np=37
      double precision, dimension(1:np) :: edf
      double precision, dimension(1:np) :: sdupper_fib
      sdupper_fib=young*p_5*k_1*max(edf/k_1,0.0d0)*
     $exp(-p_6*max(edf/k_1,0.0d0))

      return
      end subroutine dev_tens_fib

      subroutine c_shear_tens(sig_o, s_0, fsp_0,
     . !not modified
     . c_10
     . !modified
     . )
      use, intrinsic :: iso_fortran_env, only: RK => real64
      implicit none
!================================================input
      !not modified
      double precision :: c_10
     . !modified
!====================================================
      double precision :: sig_o, s_0
      double precision :: fsp_0
      double precision :: c_6_p
      c_6_p = c_10
      fsp_0 = c_6_p * max(sig_o, 0.0d0) / (1.d0 + c_6_p / s_0 *
     $max(sig_o, 0.0d0))
      return
      end subroutine c_shear_tens

! +--------------------------------------------------------------------+
! |           SUBROUTINE C_SHEAR2           |
! +--------------------------------------------------------------------+
      subroutine c_shear2(eps_L, eps_M, C_d, snf, R_S, del, dem, sl0,
     $  sm0, slf, smf, eps_V,
     . !not modified
     . young, poisson, k_1, C0, k_2, c_10, c_11, c_12,C_R2
     . !modified
     . )
      use, intrinsic :: iso_fortran_env, only: RK => real64
      implicit none
!================================================input
      !not modified
      double precision :: young, poisson, k_1, C0, k_2, 
     . c_10, c_11, c_12,C_R2
     . !modified
!====================================================
      integer, parameter ::  np=37
      double precision, dimension(1:np)  :: eps_L, eps_M, C_d, snf, del,
     $dem, sl0, sm0
      double precision  :: R_S, eps_V
      double precision, dimension(1:np) :: slf, smf
      double precision :: E_T, c_12_p, c_13_p, c_6_p, s_0
      double precision, dimension(1:np) :: C_t, fsp
      double precision, dimension(1:np) :: stf, fsp_0, epsf, depsf
      double precision :: RATE_FACTOR_1, C_S0, sig_o
      integer :: isize, j

      isize = size(snf)

      E_T=young/(1.d0+poisson)*(1.d0-4.d0*poisson)/(1.d0-2.d0*poisson)
      c_6_p = c_10
      c_12_p = c_11
      c_13_p = c_12*3.371d-4

      s_0 = k_1*k_2*E_T

      C_S0 = C_R2

      RATE_FACTOR_1 = C_S0*R_S

      sig_o = max(E_T * k_1 *
     $(c_12_p - c_13_p*max(eps_V, 0.0d0)/k_1),0.0d0)

      fsp = c_6_p * max(-snf + sig_o, 0.0d0) / (1.d0 + c_6_p / s_0 *
     $max(-snf + sig_o, 0.0d0)) * (1.d0 + RATE_FACTOR_1)

      epsf = sqrt(eps_L*eps_L + eps_M*eps_M)
      depsf = sqrt(del*del + dem*dem)
      do j=1, isize
      call c_shear_tens(sig_o, s_0, fsp_0(j),
     . !not modified
     . c_10
     . !modified
     . )
      fsp_0(j)= fsp_0(j) * (1.d0 + RATE_FACTOR_1)
      end do

      do j=1,isize
      if(sl0(j)*del(j) < 0.0d0 .or. sm0(j)*dem(j) < 0.0d0) then
      C_t(j) = C_d(j)
      else
      C_t(j) = E_T
      end if
      end do

      slf = sl0*C0 + C_t*del
      smf = sm0*C0 + C_t*dem
      stf = sqrt((slf*slf + smf*smf))
      do j = 1, isize
      if (stf(j) /= 0.0d0) then
      slf(j) = slf(j) / stf(j)
      smf(j) = smf(j) / stf(j)
      end if
      end do

      do j=1, isize
      if (snf(j) < 0.0d0 ) then
      if (stf(j) > fsp(j)) then
      stf(j) = fsp(j)
      end if
      else
      if (stf(j) > fsp_0(j)) then
      stf(j) = fsp_0(j)
      end if
      end if
      slf(j) = stf(j) * slf(j)
      smf(j) = stf(j) * smf(j)
      end do

      return
      end subroutine c_shear2

C +--------------------------------------------------------------------+
C |         SUBROUTINE C_VOL          |
C +--------------------------------------------------------------------+
      subroutine c_vol(dev, ev0, sv0, deps_N, eps_N, R_N, svneg,
     . !not modified
     . young, poisson, k_1,k_3, k_4, k_6, k_7, sv0_p, C_R2,
     . !modified
     . MaxPstrain,MinPstrain)
      use, intrinsic :: iso_fortran_env, only: RK => real64
      implicit none
!================================================input
      !not modified
      double precision :: young, poisson, k_1,k_3, k_4,
     . k_6, k_7, sv0_p, C_R2,
     . !modified
     . MaxPstrain,MinPstrain
!====================================================
      double precision fb
      integer, parameter ::  np=37
      double precision :: dev, ev0, sv0
      double precision, dimension(1:np) :: eps_N, deps_N
      double precision :: svneg,R_N,Cv0, evf, svb_neg,xk0, 
     . e0_V, prStrainDiff, xk4
      xk0 = k_3*k_1*young
      e0_V = k_4*k_1
      Cv0 = young / (1.0d0 - 2.0d0 * poisson)
      MaxPstrain = maxval(eps_N+deps_N,1)
      MinPstrain = minval(eps_N+deps_N,1)
      prStrainDiff = MaxPstrain - MinPstrain

      xk4 = (k_6*(prStrainDiff/k_1)**k_7)
     $/(1.0d0 + min(max(-sv0,0.d0),sv0_p)/Cv0) + k_4

      e0_V = xk4*k_1

      evf = ev0+dev

      svb_neg = fb(0,evf,xk0,e0_V)
      svneg=svb_neg*(1.0d0 + C_R2*R_N)

      end subroutine c_vol

C +--------------------------------------------------------------------+
C |          FUNCTION FB           |
C +--------------------------------------------------------------------+
      function fb(i,ef,xk0,e0_V)
      use, intrinsic :: iso_fortran_env, only: RK => real64
      implicit none
      integer ::  i
      double precision :: ef,xk0,e0_V
      double precision :: fb, aux

      aux=-xk0*exp(-ef/e0_V)
      if (i.eq.0) then
      fb = aux
      else if (i.eq.1) then
            fb = -aux/e0_V
      else
            print *, "ERROR at 610 while computing fb function"
            fb = 0.0
      endif
      return
      end function fb

C **********************************************************************
C *** SUBROUTINE SETSYSTEM *********************************************
C **********************************************************************

      subroutine setsystem(qn,ql,qm,w)
      use, intrinsic :: iso_fortran_env, only: RK => real64
      implicit none
      double precision, parameter :: PI=3.1415926535897932384626433832d0
      integer, parameter ::  np=37, nvhm=5, nvhf=2, nvhi=np*nvhm+nvhf
      double precision, dimension(1:6,1:np) :: qn, ql, qm
      double precision, dimension(1:np)     :: w
      integer :: jp, ij(1:2,1:6), i, j, k
      double precision, dimension(1:4,1:np) :: te
      double precision, dimension(1:3) :: xn, xm, xl, rand_vec
      double precision :: lengthn, lengthm, lengthl
      double precision, dimension(1:100,1:3) :: rand_list
      integer :: rand_counter
C
      te=0.0d0

      ij=reshape((/1,1,2,2,3,3,1,2,2,3,3,1/),(/2,6/))

      if(np.eq.37) then

      te = reshape(
     $(/0.000000000000D+00,0.000000000000D+00,1.000000000000D+00,
     $1.072388573030D-02,
     $0.000000000000D+00,1.000000000000D+00,0.000000000000D+00,
     $1.072388573030D-02,
     $1.000000000000D+00,0.000000000000D+00,0.000000000000D+00,
     $1.072388573030D-02,
     $0.000000000000D+00,7.071067811870D-01,7.071067811870D-01,
     $2.114160951980D-02,
     $0.000000000000D+00,-7.071067811870D-01,7.071067811870D-01,
     $2.114160951980D-02,
     $7.071067811870D-01,0.000000000000D+00,7.071067811870D-01,
     $2.114160951980D-02,
     $-7.071067811870D-01,0.000000000000D+00,7.071067811870D-01,
     $2.114160951980D-02,
     $7.071067811870D-01,7.071067811870D-01,0.000000000000D+00,
     $2.114160951980D-02,
     $-7.071067811870D-01,7.071067811870D-01,0.000000000000D+00,
     $2.114160951980D-02,
     $0.000000000000D+00,3.089512677750D-01,9.510778696510D-01,
     $5.355055908370D-03,
     $0.000000000000D+00,-3.089512677750D-01,9.510778696510D-01,
     $5.355055908370D-03,
     $0.000000000000D+00,9.510778696510D-01,3.089512677750D-01,
     $5.355055908370D-03,
     $0.000000000000D+00,-9.510778696510D-01,3.089512677750D-01,
     $5.355055908370D-03,
     $3.089512677750D-01,0.000000000000D+00,9.510778696510D-01,
     $5.355055908370D-03,
     $-3.089512677750D-01,0.000000000000D+00,9.510778696510D-01,
     $5.355055908370D-03,
     $9.510778696510D-01,0.000000000000D+00,3.089512677750D-01,
     $5.355055908370D-03,
     $-9.510778696510D-01,0.000000000000D+00,3.089512677750D-01,
     $5.355055908370D-03,
     $3.089512677750D-01,9.510778696510D-01,0.000000000000D+00,
     $5.355055908370D-03,
     $-3.089512677750D-01,9.510778696510D-01,0.000000000000D+00,
     $5.355055908370D-03,
     $9.510778696510D-01,3.089512677750D-01,0.000000000000D+00,
     $5.355055908370D-03,
     $-9.510778696510D-01,3.089512677750D-01,0.000000000000D+00,
     $5.355055908370D-03,
     $8.805355183100D-01,3.351545919390D-01,3.351545919390D-01,
     $1.677709091560D-02,
     $-8.805355183100D-01,3.351545919390D-01,3.351545919390D-01,
     $1.677709091560D-02,
     $8.805355183100D-01,-3.351545919390D-01,3.351545919390D-01,
     $1.677709091560D-02,
     $-8.805355183100D-01,-3.351545919390D-01,3.351545919390D-01,
     $1.677709091560D-02,
     $3.351545919390D-01,8.805355183100D-01,3.351545919390D-01,
     $1.677709091560D-02,
     $-3.351545919390D-01,8.805355183100D-01,3.351545919390D-01,
     $1.677709091560D-02,
     $3.351545919390D-01,-8.805355183100D-01,3.351545919390D-01,
     $1.677709091560D-02,
     $-3.351545919390D-01,-8.805355183100D-01,3.351545919390D-01,
     $1.677709091560D-02,
     $3.351545919390D-01,3.351545919390D-01,8.805355183100D-01,
     $1.677709091560D-02,
     $-3.351545919390D-01,3.351545919390D-01,8.805355183100D-01,
     $1.677709091560D-02,
     $3.351545919390D-01,-3.351545919390D-01,8.805355183100D-01,
     $1.677709091560D-02,
     $-3.351545919390D-01,-3.351545919390D-01,8.805355183100D-01,
     $1.677709091560D-02,
     $5.773502691900D-01,5.773502691900D-01,5.773502691900D-01,
     $1.884823095080D-02,
     $-5.773502691900D-01,5.773502691900D-01,5.773502691900D-01,
     $1.884823095080D-02,
     $5.773502691900D-01,-5.773502691900D-01,5.773502691900D-01,
     $1.884823095080D-02,
     $-5.773502691900D-01,-5.773502691900D-01,5.773502691900D-01,
     $1.884823095080D-02/)
     $,(/4,np/))

      end if

      !a list of random numbers (c++ compiler has issues finding the
      !intrinsic random generator function)
      rand_list = reshape(
     $(/0.9828567504882D-01, 0.75995391607285D+00, 0.3873291015620D+00,
     $0.24157404899597D-01, 0.81224930286407D+00, 0.11562412977219D+00,
     $0.13448733091354D+00, 0.17732083797455D-01, 0.28574603796005D+00,
     $0.48333787918091D+00, 0.24193179607391D+00, 0.72784507274628D+00,
     $0.37433540821075D+00, 0.20422506332397D+00, 0.86996459960938D+00,
     $0.11290234327316D+00, 0.74607300758362D+00, 0.27045130729675D+00,
     $0.30495953559875D+00, 0.46466767787933D-01, 0.18602609634399D+00,
     $0.94590300321579D+00, 0.93080466985703D+00, 0.56326377391815D+00,
     $0.93744999170303D+00, 0.17016625404358D+00, 0.83172667026520D+00,
     $0.63673776388168D+00, 0.90509170293808D+00, 0.80781996250153D-01,
     $0.68735098838806D+00, 0.98827022314072D+00, 0.12427037954330D+00,
     $0.41297256946564D-01, 0.58746260404587D+00, 0.22738790512085D+00,
     $0.36429065465927D+00, 0.76504719257355D+00, 0.89212298393250D-01,
     $0.29480618238449D+00, 0.15214323997498D-01, 0.46741801500320D+00,
     $0.80854523181915D+00, 0.29494619369507D+00, 0.52679270505905D+00,
     $0.82890880107880D+00, 0.82799196243286D+00, 0.40654093027115D+00,
     $0.22281312942505D+00, 0.55143457651138D+00, 0.87557244300842D+00,
     $0.53988575935364D+00, 0.79690814018250D-01, 0.45475929975510D+00,
     $0.41936284303665D+00, 0.90433490276337D+00, 0.38097238540649D+00,
     $0.21176517009735D+00, 0.45847445726395D+00, 0.37407201528549D+00,
     $0.58247053623199D+00, 0.72523069381714D+00, 0.38294112682343D+00,
     $0.79963022470474D+00, 0.64878582954407D+00, 0.91113871335983D+00,
     $0.22048246860504D+00, 0.73326307535172D+00, 0.33082985877991D+00,
     $0.54518377780914D+00, 0.68625313043594D+00, 0.65344959497452D+00,
     $0.48842376470566D+00, 0.93612748384476D+00, 0.73010480403900D+00,
     $0.14421451091766D+00, 0.28172445297241D+00, 0.38532310724258D+00,
     $0.71462720632553D+00, 0.28777182102203D-01, 0.28689461946487D+00,
     $0.33359837532043D+00, 0.58077365159988D+00, 0.86077439785004D+00,
     $0.19048988819122D+00, 0.73045420646667D+00, 0.77487111091614D-01,
     $0.92815947532654D+00, 0.96135336160660D+00, 0.81339657306671D+00,
     $0.77290034294128D+00, 0.39289391040802D+00, 0.83064794540405D-01,
     $0.81436491012573D+00, 0.27464348077774D+00, 0.69785475730896D+00,
     $0.72622752189636D+00, 0.57056963443756D+00, 0.53373098373413D+00,
     $0.89052766561508D+00, 0.10439497232437D+00, 0.44370949268341D+00,
     $0.83371174335480D+00, 0.76660126447678D+00, 0.43666344881058D+00,
     $0.52194917201996D+00, 0.74672120809555D+00, 0.47988361120224D+00,
     $0.56423664093018D-01, 0.48717916011810D-01, 0.19093745946884D+00,
     $0.55903953313828D+00, 0.89476704597473D-01, 0.60659587383270D-01,
     $0.12076753377914D+00, 0.43578207492828D-01, 0.26229143142700D+00,
     $0.54208844900131D+00, 0.67706847190857D+00, 0.34729516506195D+00,
     $0.28784900903702D+00, 0.70081567764282D+00, 0.52492463588715D+00,
     $0.65629631280899D+00, 0.28640425205231D+00, 0.98382699489594D+00,
     $0.10179662704468D+00, 0.34105694293976D+00, 0.93533742427826D+00,
     $0.44870615005493D+00, 0.90919679403305D+00, 0.36869311332703D+00,
     $0.39715743064880D+00, 0.11970198154449D+00, 0.96504551172256D+00,
     $0.49681144952774D+00, 0.22615516185760D+00, 0.46194493770599D-01,
     $0.28909325599670D-01, 0.96926856040955D+00, 0.97261154651642D+00,
     $0.17903053760529D+00, 0.96532094478607D+00, 0.29403424263000D+00,
     $0.11548119783401D+00, 0.66581213474274D+00, 0.51700073480606D+00,
     $0.98381501436234D+00, 0.37690371274948D+00, 0.33054304122925D+00,
     $0.63346850872040D+00, 0.78833806514740D+00, 0.97572487592697D+00,
     $0.31361341476440D-01, 0.35295128822327D-01, 0.77682971954346D-01,
     $0.51909381151199D+00, 0.55720746517181D+00, 0.15429109334946D+00,
     $0.11066728830338D+00, 0.56098663806915D+00, 0.69668126106262D+00,
     $0.96221327781677D-01, 0.36642336845398D+00, 0.71501332521439D+00,
     $0.97437757253647D+00, 0.35445624589920D+00, 0.83102440834045D+00,
     $0.47594881057739D+00, 0.92477875947952D+00, 0.38927114009857D+00,
     $0.37412768602371D+00, 0.17116397619247D+00, 0.81798398494720D+00,
     $0.37329262495041D+00, 0.71876001358032D+00, 0.96178060770035D+00,
     $0.71546298265457D+00, 0.14500439167023D-01, 0.86644506454468D+00,
     $0.99677217006683D+00, 0.23112386465073D+00, 0.49214643239975D+00,
     $0.68020576238632D+00, 0.40924060344696D+00, 0.26264148950577D+00,
     $0.31371957063675D+00, 0.99572634696960D+00, 0.97336316108704D+00,
     $0.35379034280777D+00, 0.41385596990585D+00, 0.11745524406433D+00,
     $0.98883545398712D+00, 0.47425985336304D+00, 0.37458121776581D+00,
     $0.30905687808990D+00, 0.97352373600006D+00, 0.76403617858887D-01,
     $0.57613718509674D+00, 0.66554796695709D+00, 0.46882772445679D+00,
     $0.39577203989029D+00, 0.86222189664841D+00, 0.75813245773315D+00,
     $0.42119151353836D+00, 0.52386933565140D+00, 0.75105327367783D+00,
     $0.74452406167984D+00, 0.29345220327377D+00, 0.59694403409958D+00,
     $0.13800948858261D+00, 0.22201073169708D+00, 0.72979748249054D-01,
     $0.57565605640411D+00, 0.34912282228470D+00, 0.24422127008438D+00,
     $0.29019010066986D+00, 0.79590600728989D+00, 0.58781945705414D+00,
     $0.42059338092804D+00, 0.63911634683609D+00, 0.44644081592560D+00,
     $0.58577197790146D+00, 0.25425022840500D+00, 0.50456845760345D+00,
     $0.95587545633316D+00, 0.54837560653687D+00, 0.99197554588318D+00,
     $0.65904629230499D+00, 0.72674357891083D+00, 0.56522768735886D+00,
     $0.36917006969452D+00, 0.71006077528000D+00, 0.66961377859116D+00,
     $0.16462564468384D+00, 0.40215849876404D-02, 0.22394019365311D+00,
     $0.12845385074615D+00, 0.58799272775650D+00, 0.21983397006989D+00,
     $0.54381072521210D-01, 0.78455865383148D-01, 0.35051292181015D+00,
     $0.91608244180679D+00, 0.30886912345886D+00, 0.17647904157639D+00,
     $0.56688421964645D+00, 0.26218175888062D+00, 0.37638366222382D+00,
     $0.73815274238586D+00, 0.44898778200150D+00, 0.37472015619278D+00,
     $0.24758487939835D+00, 0.14043205976486D+00, 0.33018738031387D+00,
     $0.33373850584030D+00, 0.16700303554535D+00, 0.89206743240356D+00,
     $0.40011823177338D+00, 0.59081101417542D+00, 0.65070927143097D-01,
     $0.85806685686111D+00, 0.82695364952087D+00, 0.24115985631943D+00,
     $0.70538848638535D+00, 0.31734997034073D+00, 0.63752508163452D+00,
     $0.62122488021851D+00, 0.21648150682449D+00, 0.67541533708572D+00,
     $0.74601423740387D+00, 0.10921978950500D+00, 0.67933183908463D+00,
     $0.60487753152847D+00, 0.88282966613770D+00, 0.26811838150024D-01,
     $0.39658391475677D+00, 0.67536926269531D+00, 0.47988468408585D+00,
     $0.18132925033569D+00, 0.58720421791077D+00, 0.85313922166824D+00,
     $0.27428126335144D+00, 0.40831637382507D+00, 0.56734585762024D+00,
     $0.88356167078018D+00, 0.66589576005936D+00, 0.97150665521622D+00,
     $0.93599849939346D+00, 0.58542084693909D+00, 0.16564524173737D+00,
     $0.68012356758118D+00, 0.19504541158676D+00, 0.55237507820129D+00,
     $0.82808983325958D+00, 0.49910986423492D+00, 0.58869194984436D+00,
     $0.38523459434509D+00, 0.79932337999344D+00, 0.28904652595520D+00/)
     $,(/100,3/))

C     ---------------------------------------
C ... Assemble tensors from direction cosines
C     ---------------------------------------
      qn = 0.0d0
      ql = 0.0d0
      qm = 0.0d0
      w = 0.0d0
      !call random_seed(size=rs_size)
      rand_counter = 0
      do jp=1,np
      rand_counter = rand_counter+1
      w(jp) = te(4,jp)*6.0D0
      xn(1) = te(3,jp)
      xn(2) = te(2,jp)
      xn(3) = te(1,jp)

      lengthm = 0.0d0
      do while (lengthm .lt. epsilon(lengthm))
      !call random_number(rand_vec)
      do k=1,3
      rand_vec(k) = rand_list(mod(rand_counter,100),k)
      end do
      xm = rand_vec - dot_product(xn,rand_vec)*xn
      lengthm = sqrt(dot_product(xm,xm))
      end do
      xm = xm/lengthm

      xl(1) = xn(2)*xm(3)-xn(3)*xm(2)
      xl(2) = xn(3)*xm(1)-xn(1)*xm(3)
      xl(3) = xn(1)*xm(2)-xn(2)*xm(1)
      lengthl = sqrt(dot_product(xl,xl))
      xl=xl/lengthl
      lengthn = sqrt(dot_product(xn,xn))
      lengthm = sqrt(dot_product(xm,xm))
      lengthl = sqrt(dot_product(xl,xl))
      do k=1,6
      i=ij(1,k)
      j=ij(2,k)
      qn(k,jp) = xn(i)*xn(j)
      qm(k,jp) = 0.5D0*(xn(i)*xm(j)+xn(j)*xm(i))
      ql(k,jp) = 0.5D0*(xn(i)*xl(j)+xn(j)*xl(i))
      end do
      end do
      return
      end subroutine setsystem

C *******************************************************************
C *** SUBROUTINE M7fMATERIAL *****************************************
C *******************************************************************

      subroutine m7fmaterial(
     $  dt, strainInc, epsOld,
     $  stressOld, stateOld, enerInternOld, enerInelasOld, damageOld,
     $  stressNew, stateNew, enerInternNew, enerInelasNew, damageNew,
     $  weightedDamage,MaxPstrain,MinPstrain,
     $  young, poisson, k_1,k_2,k_3, k_4, k_6,
     $  qn, ql, qm, w)
      
      use, intrinsic :: iso_fortran_env, only: RK => real64
      implicit none
INCLUDE 'params.inc' 
      double precision :: dt,MaxPstrain,MinPstrain
      integer, parameter ::  np=37, nvhm=5, nvhf=2, nvhi=np*nvhm+nvhf
      double precision, dimension(1:6) :: strainInc, epsOld
      double precision, dimension(1:6) :: stressOld, stressNew
      double precision, dimension(1:189) :: stateOld, stateNew
      double precision :: enerInternOld, enerInternNew
      double precision :: enerInelasOld, enerInelasNew
      double precision, dimension(1:37) :: damageOld, damageNew
      double precision :: weightedDamage
!================================================input
      !not modified
      double precision :: young, poisson, k_1,k_2,k_3, k_4,k_6
      double precision, dimension(1:6,1:np) :: qn, ql, qm
      double precision, dimension(1:np)     :: w
      !modified
      double precision :: c_1, c_5, c_6, c_7, c_8, c_9
!================================================
      double precision, dimension(1:np)     :: deps_N,
     $deps_L, deps_M, eps_N,  eps_L,  eps_M,  sig_N, sd0, 
     $sn0, ded, ed0, edf, enf, snf, snb, sdneg, sdpos, 
     $snb_fib, sl0, sm0, slf, smf, del, dem, eps_N0_neg, 
     $ eps_N0_pos,C_d, E_N
      double precision, dimension(1:nvhi+2)    :: vh_ini, vh_fin
      double precision RateFunc
      double precision :: equivStress, fractureWorkInc, smean,
     $stressPower
      double precision, parameter :: PI=3.1415926535897932384626433832d0
      double precision :: sum_bulk, bulk
      double precision :: zeta, zeta0
      double precision, dimension(1:6) :: eps, deps, sig_old
      double precision, dimension(1:6) :: sig, sig1, sig2, sig3
      double precision, dimension(1:6) :: deps_s
      double precision :: R_N
      double precision :: sv0, ev0, dev, evf, sum_snf, devv
      double precision :: svfB, svf, phi0
      integer ::  jp, NTENS,IC
      double precision :: svneg, sum1, strain_rate

      sd0 = 0.0d0
      NTENS = 6
      vh_ini = stateOld(:)
      
      DO IC = 1, NTENS
      sig(IC)=stressOld(IC)
      sig_old(IC)=stressOld(IC)
      sig(IC)=sig(IC)*unit_conv
      sig_old(IC)=sig_old(IC)*unit_conv
      IF (IC > 3) THEN
      deps(IC)=2.0d0*strainInc(IC)
      ELSE
      deps(IC)=strainInc(IC)
      END IF
      END DO

      DO IC = 1, NTENS
      IF (IC > 3) THEN ! SHEAR STRAINS ARE SHEAR ANGLES IN MICROPLANE MODEL
      eps(IC)=2.0d0*epsOld(IC)
      ELSE
      eps(IC)=epsOld(IC)
      END IF
      END DO

      deps_s(1)=deps(1)*deps(1)+deps(4)*deps(4)+deps(6)*deps(6)
      deps_s(2)=deps(4)*deps(4)+deps(2)*deps(2)+deps(5)*deps(5)
      deps_s(3)=deps(6)*deps(6)+deps(5)*deps(5)+deps(3)*deps(3)
      deps_s(4)=deps(1)*deps(4)+deps(2)*deps(4)+deps(5)*deps(6)
      deps_s(5)=deps(4)*deps(6)+deps(2)*deps(5)+deps(5)*deps(3)
      deps_s(6)=deps(1)*deps(6)+deps(4)*deps(5)+deps(6)*deps(3)
      sum1=0.d0
      do jp = 1, np
      sum1=sum1+(qn(1,jp)*deps_s(1)+qn(4,jp)*deps_s(4)+
     $qn(6,jp)*deps_s(6)+qn(4,jp)*deps_s(4)+
     $qn(2,jp)*deps_s(2)+qn(5,jp)*deps_s(5)+
     $qn(6,jp)*deps_s(6)+qn(5,jp)*deps_s(5)+
     $qn(3,jp)*deps_s(3))*w(jp)
      end do
      strain_rate=sqrt(sum1/2.d0)/dt

C     -------------------------------------------------------------------
C     Compute material response :
C     Bounding curve formulation for volumetric,
C ... deviatoric and shear variables.
C     Shear strains are shear angles.
C     -------------------------------------------------------------------

C     ---------------
C ... Initializations
C     ---------------

      jp = size(w)

      sig = 0.0d0

C     -------------------------------------------------------------------
C ... Calculate the rate coefficients using the rate function
C     -------------------------------------------------------------------

      if(on_rate_dep_frac.eq.0) then
      R_N=0.0d0
      else

      R_N = RateFunc(strain_rate/C_R1)

      end if

C     --------------------------------------------
C ... Evaluate local strains (eps_N, eps_L, eps_M)
C     --------------------------------------------

      eps_N = 0.0d0
      eps_L = 0.0d0
      eps_M = 0.0d0
      deps_N = 0.0d0
      deps_L = 0.0d0
      deps_M = 0.0d0

      call vecmatmull(eps,qn,eps_N,6,np)
      call vecmatmull(eps,ql,eps_L,6,np)
      call vecmatmull(eps,qm,eps_M,6,np)
      call vecmatmull(deps,qn,deps_N,6,np)
      call vecmatmull(deps,ql,deps_L,6,np)
      call vecmatmull(deps,qm,deps_M,6,np)

      if (maxval(eps_N)>th_del) then
      vh_fin(188)=0.0d0
      else
      vh_fin(188)=1.0d0
      end if

C     -----------------------------------------
C ... Volumetric law (same for all microplanes)
C     -----------------------------------------

      sv0 = vh_ini(1)
      phi0 = vh_ini(2)

      ev0 = (eps(1)+ eps(2)+ eps(3))/3.0d0
      dev = (deps(1)+deps(2)+deps(3))/3.0d0
      evf = ev0+dev
      zeta0 = vh_ini(189)
      zeta = zeta0

      call c_vol(dev, ev0, sv0, deps_N, eps_N, R_N, svneg,
     . !not modified
     . young, poisson, k_1,k_3, k_4, k_6, k_7, sv0_p, C_R2,
     . !modified
     . MaxPstrain,MinPstrain)

C     ---------------------
C ... Normal deviatoric law
C     ---------------------

      sn0 = vh_ini(3:nvhi:nvhm)

      ded = deps_N - dev
      ed0 = eps_N - ev0
      edf = ed0 + ded

      call c_dev(ded, ed0, dev, ev0, C_d, R_N,
     $sdneg, sdpos,
     . !not modified
     . young, poisson, k_1, c_20, C_R2, 
     . c_5_0, c_5_1, c_5_M4, c_6_0, c_6_1, c_6_2, c_7_0, 
     . c_7_1, c_7_M4, c_8_0, c_8_1, c_8_M4, c_9_0, c_9_1, 
     . c_9_2, c_9_M4, c_6_M4, p_5,p_6,V_f,cf_1,cf_2,cf_3, 
     . cf_1_1, cf_2_2,cf_3_3,f_c0,E_0,f_cp,
     . !modified
     . c_5,c_6,c_7, c_8, c_9)

      enf = evf + edf

      eps_N0_pos=vh_ini(6:nvhi:nvhm)
      eps_N0_neg=vh_ini(7:nvhi:nvhm)

      call c_norm_elastic(eps_N, deps_N, sn0, eps_N0_pos,
     $eps_N0_neg, sv0, snf, E_N, zeta, damageOld, damageNew,
     . !not modified
     . young, poisson, C0, c_18, c_19, c_20, c_21
     . !modified
     . )

      weightedDamage = 0.0d+0
      do jp = 1, np
      weightedDamage = weightedDamage + w(jp)/3.0d0*damageNew(jp)
      end do

      call c_norm(enf, sv0, R_N, snb,
     . !not modified
     . young, poisson, k_1, c_2, 
     . c_3, c_4, C_R2, V_f,d_1, d_2, d_3, d_4, d_5,
     . d_6,
     . !modified
     . c_1)

      call c_norm_fib(enf, R_N, snb_fib,
     . !not Modified
     . young, k_1, C_R2, sig_fib_0,
     . p_1, p_2, xp_1, xp_2, xp_3
     . !modified
     . )

      phi0 = vh_ini(2)

      sig_N = max(min(snf,snb+snb_fib),svneg+sdneg)

      do ic = 1, np
      if (snf(ic) > snb(ic)+snb_fib(ic)) then
      eps_N0_pos(ic) = enf(ic)
      else if (snf(ic) < svneg+sdneg(ic)) then
      eps_N0_neg(ic) = enf(ic)
      end if
      end do
      sum_snf = dot_product(sig_N,w)
      svfB = sum_snf/3.0d0
      svf = svfB
      sum_bulk = dot_product(E_N,w)
      bulk = sum_bulk/3.0d0
      if(sv0 > 0.d0) then
      if(svf > 0.d0) then
      devv = abs(dev - ((svf-sv0)/bulk))
      zeta = zeta0 + abs(devv)
      end if
      end if
C     -----------------
C ... History variables
C     -----------------
      sl0 = vh_ini(4:nvhi:nvhm)
      sm0 = vh_ini(5:nvhi:nvhm)
      snf = sig_N
      enf = eps_N + deps_N
      del = deps_L
      dem = deps_M
C     ---------
C ... Shear law
C     ---------
      call c_shear2(eps_L, eps_M, C_d, snf, R_N, del, dem, sl0, sm0,
     $slf, smf, evf,
     . !not modified
     . young, poisson, k_1, C0, k_2, c_10, c_11, c_12,C_R2
     . !modified
     . )
C     -------------------------------------------------
C     Final Stress Vector
C     -------------------------------------------------
      !sig = matmul(qn,snf*w) + matmul(qm,smf*w) + matmul(ql,slf*w)
      call matvecmull(qn,snf*w,sig1,6,np)
      call matvecmull(qm,smf*w,sig2,6,np)
      call matvecmull(ql,slf*w,sig3,6,np)
      sig = sig1 + sig2 + sig3
      sig = sig/unit_conv
C
C     -------------------------------------------
C ... Update microplane normal and shear stresses
C     -------------------------------------------
      vh_fin(1) = svf
      vh_fin(3:nvhi:nvhm) = snf
      vh_fin(4:nvhi:nvhm) = slf
      vh_fin(5:nvhi:nvhm) = smf
      vh_fin(6:nvhi:nvhm) = eps_N0_pos
      vh_fin(7:nvhi:nvhm) = eps_N0_neg
      vh_fin(189)=zeta

      stateNew(:)=vh_fin
      stressNew(:)=sig

      stressPower = (( stressOld(1)+stressNew(1) )*
     $strainInc(1) + ( stressOld(2)+stressNew(2) )*
     $strainInc(2) + ( stressOld(3)+stressNew(3) )*
     $strainInc(3) + 2.0d0*( stressOld(4)+stressNew(4) )*
     $strainInc(4) )/2.0d0

      enerInternNew = enerInternOld + stressPower

      smean = ( stressNew(1) + stressNew(2) +
     $stressNew(3) )/3.0d0
      equivStress = sqrt( 3.0d0/2.0d0 * ( (stressNew(1)-smean)**2.0d0 +
     $(stressNew(2)-smean)**2.0d0 + (stressNew(3)-smean)**2.0d0+
     $2.0d0*stressNew(4)**2.0d0 + 2.0d0*stressNew(5)**2.0d0 +
     $2.0d0*stressNew(6)**2.0d0 ) )

      fractureWorkInc = stressNew(1)*deps(1) +
     $stressNew(2)*deps(2) + stressNew(3)*deps(3) +
     $stressNew(4)*deps(4) + stressNew(5)*deps(5) +
     $stressNew(6)*deps(6)
      enerInelasNew = enerInelasOld + fractureWorkInc
      return
      end subroutine m7fmaterial

C***********************************************************************
C**** This function characterizes the theoretical rate effect      *****
C***********************************************************************
      function RateFunc(arg)
      use, intrinsic :: iso_fortran_env, only: RK => real64
      implicit none
      double precision :: RateFunc
      double precision :: arg

      RateFunc = log(arg + sqrt(arg * arg + 1.0d0))
      return
      end function RateFunc

C **********************************************************************
C *** SUBROUTINE vecmatmul *********************************************
C **********************************************************************
      subroutine vecmatmull(A,B,C,m,n)
C+---------------------------------------------------------------------+
C|  This function is added to help the c++ compiler as it is unable to
C|  find the matmul intrinsic function of Fortran.
C+---------------------------------------------------------------------+
      use, intrinsic :: iso_fortran_env, only: RK => real64
      implicit none
      integer :: i,j,m,n
      double precision, dimension(m) :: A
      double precision, dimension(m, n) :: B
      double precision, dimension(n) :: C

      do i=1,n
      C(i) = 0.000d+0
      do j=1,m
      C(i) = C(i) + A(j)*B(j,i)
      end do
      end do

      return
      end subroutine vecmatmull

C **********************************************************************
C *** SUBROUTINE matvecmul *********************************************
C **********************************************************************
      subroutine matvecmull(A,B,C,m,n)
C+---------------------------------------------------------------------+
C|  This function is added to help the c++ compiler as it is unable to
C|  find the matmul intrinsic function of Fortran.
C+---------------------------------------------------------------------+
      use, intrinsic :: iso_fortran_env, only: RK => real64
      implicit none
      integer :: i,j,m,n
      double precision, dimension(m, n) :: A
      double precision, dimension(n) :: B
      double precision, dimension(m) :: C

      do i=1,m
      C(i) = 0.000d+0
      do j=1,n
      C(i) = C(i) + A(i,j)*B(j)
      end do
      end do

      return
      end subroutine matvecmull

      !if (totalTime<=dt) then
      !call inputparams()
      !call setsystem()
      !stateOld(2)=1.0d0
      !end if
