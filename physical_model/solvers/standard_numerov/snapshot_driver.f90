
Program CCFULL_Snapshot
  ! One exact Numerov solve, appended to ccfull.f90 in place of its main
  ! program.  ccfull.inp sets the reaction, the radial grid, and the IWBC
  ! radius; snapshot.dat holds beta2,beta4,beta6,V0,r0,a0,E,J.  The program
  ! writes the grid to snapshot_meta.txt and the wave function to
  ! snapshot_wave.bin.  The propagation is complex(8); only the exported
  ! wave function is complex(4).
  Use ccfull_initialization_mod
  Implicit None
  Integer :: input_unit, meta_unit, wave_unit
  Integer :: ir, i, j, it, jt, angular_momentum
  Real(8) :: beta2_value, beta4_value, beta6_value, v0_value, r0_value, a0_value
  Real(8) :: energy_value, p, r, rh

  Call Initialize_CCFull()
  Open(Newunit=meta_unit,File='snapshot_meta.txt',Status='Replace',Action='Write')
  Write(meta_unit,*) Nlevel,R_iterat+2,R_min,dr,ReduceMass,Zpro,Ztar
  Close(meta_unit)

  Open(Newunit=input_unit,File='snapshot.dat',Status='Old',Action='Read')
  Read(input_unit,*) beta2_value,beta4_value,beta6_value,v0_value,r0_value,a0_value, &
                     energy_value,angular_momentum
  Close(input_unit)
  Beta2T=beta2_value
  Beta4T=beta4_value
  Beta6T=beta6_value
  V0=v0_value
  R0=r0_value
  A0=a0_value
  ! PotShape sets the stabilization point of Numerov_core for these V0, r0,
  ! and a0; the radial grid and the IWBC radius stay those of ccfull.inp.
  Call PotShape(0,R_barrier,V_barrier,curv,R_bottom,V_bottom)
  Do it=0,Ntar
    Do jt=0,Ntar
      bett(it,jt)=0.d0
      If (it == jt-1 .Or. it == jt+1 .Or. it == jt) bett(it,jt)=Beta2T
    End Do
  End Do
  Call coupled_matrix0()
  Do ir=0,R_iterat+1
    r=R_min+dr*ir
    Call coupled_matrix(r,CPOT0)
    Do i=1,Nlevel
      Do j=1,Nlevel
        CPOT(i,j,ir)=CPOT0(i,j)
      End Do
    End Do
  End Do
  rh=R_min+dr/2.d0
  Call coupled_matrix(rh,CPOTH)
  E=Real(energy_value)
  L_i=angular_momentum
  Call PotShape(L_i,R_barrier_l,V_barrier_l,curv_l,R_bottom_l,V_bottom_l)
  Call NumerovSnapshot(p)

  Open(Newunit=wave_unit,File='snapshot_wave.bin',Status='Replace', &
       Access='Stream',Form='Unformatted',Action='Write')
  Write(wave_unit) Cmplx(wf(1:Nx2),Kind=4)
  Close(wave_unit)
End Program CCFULL_Snapshot
