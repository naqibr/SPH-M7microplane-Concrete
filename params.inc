      double precision :: sig_fib_0, p_1, p_2, p_5, p_6
      double precision :: xp_1, xp_2, xp_3, V_f
      double precision :: cf_1, cf_2, cf_3, cf_1_1, cf_2_2, cf_3_3
      double precision :: d_1, d_2, d_3, d_4, d_5, d_6
      double precision :: k_7, sv0_p
      double precision :: c_2, c_3, c_4, c_10
      double precision :: c_11,c_12,c_18, c_19,c_20,c_21,C_R2
      double precision :: c_5_0, c_5_1, c_5_M4, c_6_0, c_6_1, c_6_2
      double precision :: c_7_0, c_7_1, c_7_M4, c_8_0, c_8_1, c_8_M4
      double precision :: c_9_0, c_9_1, c_9_2, c_9_M4, c_6_M4
      double precision :: C0,th_del,unit_conv
      integer      :: on_rate_dep_frac
      double precision :: C_R1,f_c0,E_0,f_cp

      PARAMETER(f_c0=15.08d0)
      PARAMETER(E_0=20000.d0)
      PARAMETER(f_cp=90.3)
      PARAMETER(d_1 = 0.095)
      PARAMETER(c_2 = 1.76d-1)
      PARAMETER(c_3 = 4.0d0)
      PARAMETER(c_4 = 50.0d0)
      PARAMETER(c_7_1 = 3500.0d0)
      PARAMETER(c_8_1 = 20.0d0)
      PARAMETER(c_7_M4 = 10.d-1)
      PARAMETER(c_8_M4 = 8.0)
      PARAMETER(c_5_0 = 1.3d-2)
      PARAMETER(c_10 = 3.3d-1)
      PARAMETER(c_11 = 0.5d0)
      PARAMETER(c_12 = 7.00d3)
      PARAMETER(c_19 = 4500.0d0) 
      PARAMETER(c_20 = 300.0d0)
      PARAMETER(c_18 = 4000.d0)
      PARAMETER(c_21 = 60.d0)
      PARAMETER(d_4 = 1.7d0)
      PARAMETER(d_6 = 25.5d0)
      PARAMETER(d_5 = 1000.0d0)
      PARAMETER(k_7 = 1.8d0)
!
      PARAMETER(c_6_1 = 4.0d1)
      PARAMETER(c_6_2 = 13.0d0)
      PARAMETER(c_9_1 = 4.0d1)
      PARAMETER(c_9_2 = 13.0d0)
      PARAMETER(c_6_M4 = 1.30d0)
      PARAMETER(c_9_M4 = 0.d-1)
      
      PARAMETER(on_rate_dep_frac = 0)
      PARAMETER(C_R1 = 4.0d-6) 
      PARAMETER(C_R2 = 22.0d-3) 
      PARAMETER(C0 = 1.0d0)
      PARAMETER(c_5_1 = 4.0)
      PARAMETER(c_7_0 = 1.2d-2) 
      PARAMETER(c_8_0 = 1.2d-2) 
      PARAMETER(c_6_0 = 4.0d2) 
      PARAMETER(c_9_0 = 4.0d2) 
      PARAMETER(c_5_M4 = 3.0d0) 

      PARAMETER(V_f = 0.0d0)
      PARAMETER(d_2 = 35.0d0)
      PARAMETER(d_3 = 1.7d0) !1.7
      PARAMETER(sv0_p = 250.d0)
      PARAMETER(unit_conv = 1.0d0)
      PARAMETER(th_del = 0.005d10) 
      PARAMETER(cf_2_2 = 1.5d+1)
      PARAMETER(cf_1 = 104.0d-2)
      PARAMETER(cf_2 = 54.0d+1)
      PARAMETER(cf_1_1 = 3.5d+1)
      PARAMETER(cf_3 = 104.0d-2)
      PARAMETER(cf_3_3 = 3.5d+1)
      PARAMETER(sig_fib_0 = -4.5d0)
      PARAMETER(p_1 = 0.0d0)
      PARAMETER(p_2 = 300.d-3)
      PARAMETER(p_5 = 0.0d0)
      PARAMETER(p_6 = 0.0d0)
      PARAMETER(xp_1 = 3.34d0)
      PARAMETER(xp_2 = 0.0d0)
      PARAMETER(xp_3 = 1.0d0)