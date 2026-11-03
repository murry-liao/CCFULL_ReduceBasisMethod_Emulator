
Subroutine CCBROM_Noninteractive_Setup()
  Use ccfull_initialization_mod
  Implicit None
  Integer :: it,jt,ip,jp

  ! The linked prototype supports the standard axial target rotor with an
  ! inert projectile.  Initialize exactly the harmonic/default values that the
  ! interactive Output_information -> Anharmonicity/grotation path would set.
  Ntrans=0
  iq=0
  NphononT=0
  NphononT2=0
  LambdaT=0
  LambdaT2=0
  OmegaT=0.d0
  OmegaT2=0.d0
  BetaT=0.d0
  BetaTn=0.d0
  BetaT2=0.d0
  BetaT2n=0.d0
  OmegaP=0.d0
  BetaP=0.d0
  BetaPn=0.d0
  E2P=0.d0
  Beta2P=0.d0
  Beta4P=0.d0
  imutual=0
  imutual(0,0)=1
  betnahv=0.d0
  betcahv=0.d0
  betnahv2=0.d0
  betcahv2=0.d0
  betnahvp=0.d0
  betcahvp=0.d0
  omeahv=0.d0
  omeahv2=0.d0
  omeahvp=0.d0
  bett=0.d0
  erott=0.d0
  betp=0.d0
  erotp=0.d0

  Do it=0,Ntar
    Do jt=0,Ntar
      If (it == jt-1 .Or. it == jt+1 .Or. it == jt) bett(it,jt)=Beta2T
    End Do
    erott(it)=2.d0*it*(2.d0*it+1.d0)/6.d0*E2T
  End Do
  Do ip=0,Npro
    Do jp=0,Npro
      If (ip == jp-1 .Or. ip == jp+1 .Or. ip == jp) betp(ip,jp)=Beta2P
    End Do
    erotp(ip)=2.d0*ip*(2.d0*ip+1.d0)/6.d0*E2P
  End Do
End Subroutine CCBROM_Noninteractive_Setup

Subroutine ccbrom_initialize(status,nlevel_out,nradial_out,rmin_out,dr_out) Bind(C)
  Use Iso_C_Binding
  Use ccfull_initialization_mod
  Implicit None
  Integer(C_Int), Intent(Out) :: status,nlevel_out,nradial_out
  Real(C_Double), Intent(Out) :: rmin_out,dr_out
  status=0
  Call Initialize_CCFull()
  If (IVIBROTT /= 1 .Or. IVIBROTP /= -1 .Or. Npro /= 0) status=1
  nlevel_out=Nlevel
  nradial_out=R_iterat+2
  rmin_out=R_min
  dr_out=dr
End Subroutine ccbrom_initialize

Subroutine ccbrom_evaluate_fields(beta2_value,beta4_value,beta6_value,v0_value,r0_value,a0_value,npoint,radii,values,status) Bind(C)
  Use Iso_C_Binding
  Use ccfull_initialization_mod
  Implicit None
  Real(C_Double), Value :: beta2_value,beta4_value,beta6_value,v0_value,r0_value,a0_value
  Integer(C_Int), Value :: npoint
  Real(C_Double), Intent(In) :: radii(npoint)
  Real(C_Double), Intent(Out) :: values(*)
  Integer(C_Int), Intent(Out) :: status
  Integer :: point,i,j,it,jt,offset
  status=0
  Beta2T=beta2_value
  Beta4T=beta4_value
  Beta6T=beta6_value
  V0=v0_value
  R0=r0_value
  A0=a0_value
  Do it=0,Ntar
    Do jt=0,Ntar
      bett(it,jt)=0.d0
      If (it == jt-1 .Or. it == jt+1 .Or. it == jt) bett(it,jt)=Beta2T
    End Do
  End Do
  Call coupled_matrix0()
  Do point=1,npoint
    Call coupled_matrix(radii(point),CPOT0)
    offset=(point-1)*Nlevel*Nlevel
    Do i=1,Nlevel
      Do j=1,Nlevel
        values(offset+(i-1)*Nlevel+j)=CPOT0(i,j)
      End Do
    End Do
  End Do
End Subroutine ccbrom_evaluate_fields

Subroutine ccbrom_evaluate_factors(beta2_value,beta4_value,beta6_value,v0_value,r0_value,a0_value, &
    npoint,radii,nmode,decomposition_inverse,include_nuclear,values,status) Bind(C)
  Use Iso_C_Binding
  Use ccfull_initialization_mod
  Implicit None
  Real(C_Double), Value :: beta2_value,beta4_value,beta6_value,v0_value,r0_value,a0_value
  Integer(C_Int), Value :: npoint,nmode,include_nuclear
  Real(C_Double), Intent(In) :: radii(npoint),decomposition_inverse(*)
  Real(C_Double), Intent(Out) :: values(*)
  Integer(C_Int), Intent(Out) :: status
  Integer :: point,mode,i,j,k,it,jt,matrix_index,weight_index,value_index
  Real(8) :: factor,bare,weight,quadrupole_geometry,hexadecapole_geometry,beta6_geometry
  Real(8) :: nuclear_weights(nmode,Nlevelmax),bare_weights(nmode)
  Real(8) :: quadrupole_weights(nmode),hexadecapole_weights(nmode),beta6_weights(nmode)
  Real(8) :: Vn,VnCC,Fct2,Fct4,Fct6,CG
  External Vn,VnCC,Fct2,Fct4,Fct6,CG

  status=0
  If (npoint < 1 .Or. nmode < 1) Then
    status=2
    Return
  End If
  Beta2T=beta2_value
  Beta4T=beta4_value
  Beta6T=beta6_value
  V0=v0_value
  R0=r0_value
  A0=a0_value
  Do it=0,Ntar
    Do jt=0,Ntar
      bett(it,jt)=0.d0
      If (it == jt-1 .Or. it == jt+1 .Or. it == jt) bett(it,jt)=Beta2T
    End Do
  End Do
  Call coupled_matrix0()

  ! Contract the fixed channel-space decomposition with CCFULL's eigenchannel
  ! representation once per parameter point.  For the standard axial target
  ! rotor, coupled_matrix constructs the nuclear field as
  !
  !   sum_k Evec(i,k) VnCC(r,Ev(k)) Evec(j,k),
  !
  ! so its projection can be evaluated at any radius without first forming the
  ! complete Nlevel-by-Nlevel matrix.  The two Coulomb contractions below use
  ! exactly the same geometries and beta-dependent factors as coupled_matrix.
  nuclear_weights=0.d0
  bare_weights=0.d0
  quadrupole_weights=0.d0
  hexadecapole_weights=0.d0
  beta6_weights=0.d0
  If ((Nlevel-Ntrans) > 1) Then
    Do mode=1,nmode
      Do i=1,Nlevel
        Do j=1,Nlevel
          matrix_index=(i-1)*Nlevel+j
          weight_index=(mode-1)*Nlevel*Nlevel+matrix_index
          weight=decomposition_inverse(weight_index)
          If (i == j) bare_weights(mode)=bare_weights(mode)+weight

          ! coupled_matrix fills the upper triangle and copies it to the lower
          ! triangle.  Use the same ordered channel pair for identical angular
          ! momentum algebra in both halves of the contraction.
          it=Min(i,j)-1
          jt=Max(i,j)-1
          quadrupole_geometry= &
              Sqrt((2*2*it+1)*5.d0*(2*2*jt+1)/(4.d0*Pi)) &
              *CG(2*it,0,2,0,2*jt,0)**2/(2.d0*2*jt+1) &
              *(bett(it,jt)+2.d0*Sqrt(5.d0/Pi)*bett(it,jt)**2/7.d0)
          hexadecapole_geometry= &
              Sqrt((2*2*it+1)*9.d0*(2*2*jt+1)/(4.d0*Pi)) &
              *CG(2*it,0,4,0,2*jt,0)**2/(2.d0*2*jt+1)
          beta6_geometry= &
              Sqrt((2*2*it+1)*13.d0*(2*2*jt+1)/(4.d0*Pi)) &
              *CG(2*it,0,6,0,2*jt,0)**2/(2.d0*2*jt+1)
          quadrupole_weights(mode)=quadrupole_weights(mode) &
              +weight*quadrupole_geometry
          hexadecapole_weights(mode)=hexadecapole_weights(mode) &
              +weight*hexadecapole_geometry
          beta6_weights(mode)=beta6_weights(mode)+weight*beta6_geometry

          Do k=1,Nlevel-Ntrans
            nuclear_weights(mode,k)=nuclear_weights(mode,k) &
                +weight*Evec(i,k)*Evec(j,k)
          End Do
        End Do
      End Do
    End Do
  End If

  Do point=1,npoint
    If (radii(point) <= 0.d0) Then
      status=3
      Return
    End If
    bare=Vn(radii(point))
    Do mode=1,nmode
      If ((Nlevel-Ntrans) == 1) Then
        factor=0.d0
        If (include_nuclear /= 0) factor= &
            decomposition_inverse((mode-1)*Nlevel*Nlevel+1)*bare
      Else
        factor=0.d0
        Do k=1,Nlevel-Ntrans
          factor=factor+nuclear_weights(mode,k)*VnCC(radii(point),Ev(k))
        End Do
        factor=factor+quadrupole_weights(mode)*Fct2(radii(point)) &
            +hexadecapole_weights(mode)*Fct4(radii(point)) &
            +beta6_weights(mode)*Fct6(radii(point))
        If (include_nuclear == 0) factor=factor-bare_weights(mode)*bare
      End If
      value_index=(point-1)*nmode+mode
      values(value_index)=factor
    End Do
  End Do
End Subroutine ccbrom_evaluate_factors

Subroutine ccbrom_coulomb_values(ell,npoint,energies,radii,up_real,up_imag,down_real,down_imag,status) Bind(C)
  Use Iso_C_Binding
  Use ccfull_initialization_mod
  Implicit None
  Integer(C_Int), Value :: ell,npoint
  Real(C_Double), Intent(In) :: energies(npoint),radii(npoint)
  Real(C_Double), Intent(Out) :: up_real(npoint),up_imag(npoint),down_real(npoint),down_imag(npoint)
  Integer(C_Int), Intent(Out) :: status
  Integer :: point,iv
  Real(8) :: rho,eta,f,fp,g,gp
  status=0
  Do point=1,npoint
    If (energies(point) <= 0.d0) Then
      status=3
      up_real(point)=0.d0
      up_imag(point)=0.d0
      down_real(point)=0.d0
      down_imag(point)=0.d0
    Else
      rho=Sqrt(2.d0*ReduceMass*energies(point))/Hbar*radii(point)
      eta=(Zpro*Ztar/137.d0)*Sqrt(ReduceMass/(2.d0*energies(point)))
      Call myCOULFG(Real(ell,8),eta,rho,f,fp,g,gp,iv)
      up_real(point)=g
      up_imag(point)=f
      down_real(point)=g
      down_imag(point)=-f
    End If
  End Do
End Subroutine ccbrom_coulomb_values

Subroutine ccbrom_exact_penetrabilities(beta2_value,beta4_value,beta6_value,v0_value,r0_value,a0_value, &
    energy_value,nwave,angular_momenta,values,status) Bind(C)
  Use Iso_C_Binding
  Use ccfull_initialization_mod
  Implicit None
  Real(C_Double), Value :: beta2_value,beta4_value,beta6_value,v0_value,r0_value,a0_value
  Real(C_Double), Value :: energy_value
  Integer(C_Int), Value :: nwave
  Integer(C_Int), Intent(In) :: angular_momenta(nwave)
  Real(C_Double), Intent(Out) :: values(nwave)
  Integer(C_Int), Intent(Out) :: status
  Integer :: wave,ir,i,j,it,jt
  Real(8) :: p,r,rh

  status=0
  If (energy_value <= 0.d0 .Or. nwave < 1) Then
    status=3
    Return
  End If
  Beta2T=beta2_value
  Beta4T=beta4_value
  Beta6T=beta6_value
  V0=v0_value
  R0=r0_value
  A0=a0_value
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
  Do wave=1,nwave
    If (angular_momenta(wave) < 0) Then
      status=4
      Return
    End If
    L_i=angular_momenta(wave)
    Call PotShape(L_i,R_barrier_l,V_barrier_l,curv_l,R_bottom_l,V_bottom_l)
    Call NumerovSnapshot(p)
    values(wave)=p
  End Do
End Subroutine ccbrom_exact_penetrabilities

Subroutine ccbrom_active_partial_waves(energy_value,v0_value,r0_value,a0_value,lmax,nactive) Bind(C)
  Use Iso_C_Binding
  Use ccfull_initialization_mod
  Implicit None
  Real(C_Double), Value :: energy_value,v0_value,r0_value,a0_value
  Integer(C_Int), Value :: lmax
  Integer(C_Int), Intent(Out) :: nactive
  Integer :: ell
  Real(8) :: V
  External V
  V0=v0_value
  R0=r0_value
  A0=a0_value
  nactive=0
  Do ell=0,lmax
    L_i=ell
    Call PotShape(L_i,R_barrier_l,V_barrier_l,curv_l,R_bottom_l,V_bottom_l)
    If (R_bottom_l < 0.d0 .Or. V(R_bottom_l,L_i) > energy_value) Exit
    nactive=nactive+1
  End Do
End Subroutine ccbrom_active_partial_waves
