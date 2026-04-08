!'''
!@File    :   CC_python.py Ec_version solved by modified Numerove or Discrete Basis method 
!@Time    :   2025/06/22 
!@Author  :   Liao ZeHong, Kouichi Hagino
!@Version :   1.0
!@Contact :   liaozh26@mail2.sysu.edu.cn
!'''


Module ccfull_initialization_mod
    Implicit None   
  !  """Principal global variables
  !  P                   - penetrability 
  !  Sigma               - fusion cross section, unit mb
  !  Spin                - mean angular momentum
  !    
  !  Apro, Zpro, Rpro    - atomic #, proton # and radius of the projectile
  !  Atar, Ztar, Rtar    - those of the target
  !  ReduceMass          - reduced mass, unit MeV/c**2
  !  Amu                 - nucleon mass, unit MeV/c**2
  !  E                   - bombarding energy in the center of mass frame, unit MeV
  !    
  !  V0,R0,A0            - depth, range, and dissuseness parameters of uncoupled 
  !                        nuclear potential, which is assumed to be a Woods-Saxon form
  !    
  !  IVIBROTT (IVIBROTP) - option for intrinsic degree of freedom 
  !                        = -1; no excitation (inert)
  !                        =  0; vibrational coupling 
  !                        =  1; rotational coupling
  !
  !  Ntar (Npro)         - the number of levels to be included    
  !
  !  BetaTn (BetaPn)     - i am not very clear about it (Liao)
  !  BetaT (BetaP)       - defeormation parameter
  !  OmegaT (OmegaP)     - excitation energy of the oscillator
  !  LambdaT (LambdaP)   - multipolarity
  !  NphononT (NphononP) - the number of phonon to be included
  !  Beta2T (Beta2P)     - static quadrupole deformation parameter
  !  Beta4T (Beta4P)     - static hexadecapole deformation parameter
  !  BetaT2n             - is the same as BetaTn but for the second mode
  !  BetaT2              - is the same as BETAT but for the second mode
  !  OmegaT2             - is the same as OmegaT but for the second mode
  !  LambdaT2            - is the same as LambdaT but for the second mode
  !  NphononT2           - is the same as NphononT but for the second mode
  !
  !  E2T (E2P)           - excitation energy of 2+ state in a rotational band
  !  NrotT (NrotP)       - the number of levels in the rotational band to be included 
  !                        (up to I^pi=2*NROT+ states are included)
  !
  !  L                   - angular momentum of the relative motion
  !  CPOT                - coupling matrix
  !  """
  Real(8), parameter :: Hbar = 197.329d0
  Real(8), parameter :: Pi   = 3.141592653d0
  Real(8), parameter :: Amu  = 938.0
  Real(8), parameter :: e2   = 1.44

  integer, parameter :: Nlevelmax = 30
  integer, Public ::  Apro, Zpro, Atar, Ztar, L_i
  Real(8), Public ::  R0P, R0T, Rpro, Rtar, ReduceMass, h2m
  Real(8), Public ::  OmegaT, BetaT, OmegaT2, BetaT2, E2T, Beta2T, Beta4T
  Real(8), Public ::  OmegaP, BetaP, E2P, Beta2P, Beta4P
  Real(8), Public ::  V0, R0, A0, Emin, Emax, dE, R_max, R_min, dr
  Real(8), Public ::  BetaTn, BetaPn, BetaT2n, BetaP2n
  Real(8), Public ::  R_barrier, V_barrier, curv, R_bottom, V_bottom
  Real(8), Public ::  R_barrier_l, V_barrier_l, curv_l, R_bottom_l, V_bottom_l
  Real(8), Public ::  t
  Real,    Public ::  Ftr, Qtrans
  Real,    Public ::  E
  Integer, Public ::  IVIBROTP, IVIBROTT, LambdaT, NphononT, LambdaT2, NphononT2
  Integer, Public ::  NrotT, Ntar, LambdaP, NphononP, NrotP, Npro
  Integer, Public ::  Ntrans, iq
  Integer, Public ::  Nlevel, R_iterat
  Real(8), Allocatable, Public :: Pot_para1(:), Pot_para2(:), Pot_para3(:)
  Real(8), Allocatable, Public :: Ev(:), eps(:)
  Real(8), Allocatable, Public :: omeahv(:), omeahv2(:), omeahvp(:), erott(:), erotp(:)
  Real(8), Allocatable, Public :: betnahv(:,:), betcahv(:,:), betnahv2(:,:), betcahv2(:,:)
  Real(8), Allocatable, Public :: betnahvp(:,:), betcahvp(:,:), bett(:,:), betp(:,:)
  Real(8), Allocatable, Public :: Evec(:,:), CPOT0(:,:), CPOTH(:,:), CPOT(:,:,:)
  
  Integer, Allocatable, Public :: imutual(:,:)


  Integer, Public ::  NumEmulator, Nx, Nx_cc, Nx1, Nx2
  Real(8),    Allocatable, Public ::  Echannel(:)
  complex*16, Allocatable, Public :: wf(:)
  complex*16, Allocatable, Public :: wf_base(:,:), wf_base2(:,:), wf_base3_dagger(:,:)
  complex*16, Allocatable, Public :: Hamilton(:,:), unit(:,:)
  complex*16, Allocatable, Public :: Ec_psi(:,:), Ec_psi0(:), Ec_psi0_2(:)
  complex*16, Allocatable, Public :: cc_ec(:,:), cc_ec_inv(:,:),dd_ec(:), dd2_ec(:), a_tra(:)

  Real(8), allocatable :: training_set1(:), training_set2(:)
  Real(8), allocatable :: testing_set1(:), testing_set2(:)

  Real :: start_time1, end_time1, time_inteval_tranning
  Real :: start_time2, end_time2, time_inteval_ture
  Real :: start_time3, end_time3, time_inteval_Ec
  Integer ::  num_traning_set1, num_traning_set2
  Integer ::  Num_testing, num_testing_set1, num_testing_set2
  
  character(len=2), dimension(0:109) :: Element

  data Element / &
       'O ','H ','He','Li','Be','B ','C ','N ','O ','F ','Ne', &
       'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca','Sc','Ti','V ', &
       'Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr', &
       'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In', &
       'Sn','Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm', &
       'Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W ','Re', &
       'Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra', &
       'Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md', &
       'No','Lr','XX','X1','X2','X3','X4','04' /

Contains

  Subroutine Read_Input()

    Open(10, File='ccfull.inp', Status = 'Unknown')
    !Open(20, File='sigma.dat',  Status = 'Unknown')
    !Open(30, File='Output.dat', Status = 'Unknown')
    !Open(40, File='Wavefunction_Db.dat', Status = 'Unknown')
    !Open(50, File='Wavefunction_Nu.dat', Status = 'Unknown')
    Read(10,*)Apro, Zpro, Atar, Ztar
    Read(10,*)R0P, IVIBROTP, R0T, IVIBROTT

    Rpro = R0P * Apro ** (1.d0 / 3.d0)
    Rtar = R0T * Atar ** (1.d0 / 3.d0)
    ReduceMass = Apro * Atar * 1.d0 / (Apro + Atar) * Amu

    h2m = Hbar**2 / (2 * ReduceMass)
    
    If (IVIBROTT == 0) Then

      Read(10,*)OmegaT, BetaT, LambdaT, NphononT
      Read(10,*)OmegaT2, BetaT2, LambdaT2, NphononT2
      Ntar = NphononT

    Else If (IVIBROTT == 1) Then

      Read(10,*)E2T, Beta2T, Beta4T, NrotT 
      Read(10,*)
      Ntar = NrotT
    Else
      Read(10,*)
      Read(10,*)
      Ntar = 0
    End If

    If (IVIBROTP == 0) Then

      Read(10,*)OmegaP, BetaP, LambdaP, NphononP
      Npro = NphononP

    Else If (IVIBROTP == 1) Then

      Read(10,*)E2P, Beta2P, Beta4P, NrotP
      Npro = NrotP
    Else
      Read(10,*)
      Npro = 0
    End If
    Nlevel=(Npro + 1)*(Ntar + 1)

    Read(10,*)Qtrans, Ftr
    Read(10,*)V0, R0, A0
    Read(10,*)Emin, Emax, dE
    Read(10,*)R_max, dr
  End Subroutine Read_Input

  Subroutine Output_information()
      Character(len=1) :: ans
      !    print '(A, I0 ,A, A,I0 ,A)', 'System: ', Int(Apro), Element(Zpro),'+', Int(Atar), Element(Ztar)
    !
      !    print '(A, F8.3, A, F8.3, A, F8.3, A, F8.3, A)', &
      !          "Simulation range E from ", Emin, " MeV To ", Emax, " MeV, dE = ", dE, " MeV"
    !
      !    print '(A, F8.4, A, F8.4, A, F8.3, A, F8.3, A)', &
      !          "Simulation range R from ", R_min, " fm To ", R_max, " fm, dR = ", dr, " fm"
    !
      !    print '(A, F8.3, A, F8.3, A, F8.3, A, F8.3, A)', &
      !          "Potential parameters: V0= ", V0, "(MeV), a= ", A0, "(fm), r0= ", R0, "(fm), R= ", R0*(Apro**(1.0/3.0) + Atar**(1.0/3.0)), "(fm)"
    !
      !    print '(A, F8.4, A)', "Coulomb barrier position :", R_barrier, "fm"
    !
      !    print '(A, F8.4, A)', "Coulomb barrier energy   :", V_barrier, "MeV"
    !
      !    print '(A, F8.4, A)', "Coulomb barrier curv     :", curv, "MeV"
    !
      !    print '(A, F8.4, A)', "Coulomb bottom position  :", R_bottom, "fm"
    !
      !    print '(A, F10.4, A)', "Coulomb bottom energy    :", V_bottom, "MeV"

      if (Ntar /= 0) then
          if (IVIBROTT == 0) then
              !write(*,'(A, F6.3, A, F6.3, A, I0, A, I0)') "Phonon Excitation in the targ.: beta=", BetaT, ", omega=", OmegaT, " (MeV), Lambda=", LambdaT, ", Nph=", NphononT
              BetaTn = BetaT
              !Write(*,*)' Different beta_N from beta_C for this mode(n/y)?'
              !Read(*,*)ans
              ans = 'n'
              IF(ans .Eq. 'Y' .Or. ans .Eq. 'y') Then
                Write(*,*)'beta_N=?'
                Read(*,*)BetaTn
              End If
              !write(*,'(A, F6.3, A, F6.3, A, F6.3, A)') "Phonon Excitation in the targ.: beta_N=", BetaTn, ", beta_C=", BetaT, ", r0=", R0T, " (fm),"
              !write(*,'(A, F6.3, A, I0, A, I0)') "                              omega=", OmegaT, " (MeV), Lambda=", LambdaT, ", Nph=", NphononT
              !write(30,'(A, F6.3, A, F6.3, A, F6.3, A)') "Phonon Excitation in the targ.: beta_N=", BetaTn, ", beta_C=", BetaT, ", r0=", R0T, " (fm),"
              !write(30,'(A, F6.3, A, I0, A, I0)') "                              omega=", OmegaT, " (MeV), Lambda=", LambdaT, ", Nph=", NphononT
          else if (IVIBROTT == 1) then
              !write(*,'(A, F6.3, A, F6.3, A, F6.3, A)') "Rotational Excitation in the targ.: beta2=", Beta2T, ", beta4=", Beta4T, ", r0=", R0T, " (fm),"
              !write(*,'(A, F6.3, A, I0)') "                                   E2=", E2T, " (MeV), Nrot=", NrotT
              !write(30,'(A, F6.3, A, F6.3, A, F6.3, A)') "Rotational Excitation in the targ.: beta2=", Beta2T, ", beta4=", Beta4T, ", r0=", R0T, " (fm),"
              !write(30,'(A, F6.3, A, I0)') "                                   E2=", E2T, " (MeV), Nrot=", NrotT
          end if
      end if

      if (NphononT2 /= 0) then
          !write(*,'(A, F6.3, A, F6.3, A, I0, A, I0)') "Phonon Excitation in the targ.: beta=", BetaT2, ", omega=", OmegaT2, " (MeV), Lambda=", LambdaT2, ", Nph=", NphononT2
          BetaT2n = BetaT2
          !Write(*,*)' Different beta_N from beta_C for this mode(n/y)?'
          Read(*,*)ans
          IF(ans .Eq. 'Y' .Or. ans .Eq. 'y') Then
            Write(*,*)'beta_N=?'
            Read(*,*)BetaT2n
          End If
          !write(*,'(A, F6.3, A, F6.3, A, F6.3, A)') "Phonon Excitation in the targ.: beta_N=", BetaT2n, ", beta_C=", BetaT2, ", r0=", R0T, " (fm),"
          !write(*,'(A, F6.3, A, I0, A, I0)') "                              omega=", OmegaT2, " (MeV), Lambda=", LambdaT2, ", Nph=", NphononT2
          !write(30,'(A, F6.3, A, F6.3, A, F6.3, A)') "Phonon Excitation in the targ.: beta_N=", BetaT2n, ", beta_C=", BetaT2, ", r0=", R0T, " (fm),"
          !write(30,'(A, F6.3, A, I0, A, I0)') "                              omega=", OmegaT2, " (MeV), Lambda=", LambdaT2, ", Nph=", NphononT2
          Call Mutual
      End if

      !print *, "------------------------------"
      !print *, "Mode of excitation for projectile"

      if (Npro /= 0) then
          if (IVIBROTP == 0) then
              write(*,'(A, F6.3, A, F6.3, A, I0, A, I0)') "Phonon Excitation in the proj.: beta=",&
                           BetaP, ", omega=", OmegaP, " (MeV), Lambda=", LambdaP, ", Nph=", NphononP
              BetaPn = BetaP
              Write(*,*)' Different beta_N from beta_C for this mode(n/y)?'
              Read(*,*)ans
              IF(ans .Eq. 'Y' .Or. ans .Eq. 'y') Then
                Write(*,*)'beta_N=?'
                Read(*,*)BetaPn
              End If
              write(*,'(A, F6.3, A, F6.3, A, F6.3, A)') "Phonon Excitation in the proj.: beta_N=",&
                                                 BetaPn, ", beta_C=", BetaP, ", r0=", R0P, " (fm),"
              write(*,'(A, F6.3, A, F6.3, A, I0)') "omega=", OmegaP, " (MeV), Lambda=", LambdaP, ", Nph=", NphononP
              write(30,'(A, F6.3, A, F6.3, A, F6.3, A)') "Phonon Excitation in the proj.: beta_N=",&
                                                     BetaPn, ", beta_C=", BetaP, ", r0=", R0P, " (fm),"
              write(30,'(A, F6.3, A, F6.3, A, I0)') "omega=", OmegaP, " (MeV), Lambda=", LambdaP, ", Nph=", NphononP
          else if (IVIBROTP == 1) then
              write(*,'(A, F6.3, A, F6.3, A, F6.3, A)') "Rotational Excitation in the proj.: beta2=",&
                                                   Beta2P, ", beta4=", Beta4P, ", r0=", R0P, " (fm),"
              write(*,'(A, F6.3, A, I0)') " E2=", E2P, " (MeV), Nrot=", NrotP
              write(30,'(A, F6.3, A, F6.3, A, F6.3, A)') "Rotational Excitation in the proj.: beta2=",&
                                                     Beta2P, ", beta4=", Beta4P, ", r0=", R0P, " (fm),"
              write(30,'(A, F6.3, A, I0)') "  E2=", E2P, " (MeV), Nrot=", NrotP
          end if
      end if
      Call Anharmonicity
      Call grotation

      !print *, "------------------------------"
      !print *, "Transfer coupled mode"

      if (Ntrans /= 0) then
          write(*,'(A, F6.3, A, F6.3, A)') "Transfer channel: Strength= ",&
                                               Ftr, ", Q = ", Qtrans, " MeV"
          write(30,'(A, F6.3, A, F6.3, A)') "Transfer channel: Strength= ",&
                                               Ftr, ", Q = ", Qtrans, " MeV"
      end if
      !print *, "------------------------------"
  End Subroutine Output_information


  Subroutine Output_information_old()
      Character(len=1) :: ans


      print '(A, I0 ,A, A,I0 ,A)', 'System: ', Int(Apro), Element(Zpro),'+', Int(Atar), Element(Ztar)

      print '(A, F8.3, A, F8.3, A, F8.3, A, F8.3, A)', &
            "Simulation range E from ", Emin, " MeV To ", Emax, " MeV, dE = ", dE, " MeV"

      print '(A, F8.4, A, F8.4, A, F8.3, A, F8.3, A)', &
            "Simulation range R from ", R_min, " fm To ", R_max, " fm, dR = ", dr, " fm"

      print '(A, F8.3, A, F8.3, A, F8.3, A, F8.3, A)', &
            "Potential parameters: V0= ", V0, "(MeV), a= ", A0, "(fm), r0= ", R0, "(fm), R= ",&
                                         R0*(Apro**(1.0/3.0) + Atar**(1.0/3.0)), "(fm)"

      print '(A, F8.4, A)', "Coulomb barrier position :", R_barrier, "fm"

      print '(A, F8.4, A)', "Coulomb barrier energy   :", V_barrier, "MeV"

      print '(A, F8.4, A)', "Coulomb barrier curv     :", curv, "MeV"

      print '(A, F8.4, A)', "Coulomb bottom position  :", R_bottom, "fm"

      print '(A, F10.4, A)', "Coulomb bottom energy    :", V_bottom, "MeV"

      write(30,'(A, I0 ,A, A,I0 ,A)') 'System: ',&
               Int(Apro), Element(Zpro),'+', Int(Atar), Element(Ztar)

      write(30,'(A, F8.3, A, F8.3, A, F8.3, A, F8.3, A)') &
            "Simulation range E from ", Emin, " MeV To ", Emax, " MeV, dE = ", dE, " MeV"

      write(30,'(A, F8.4, A, F8.4, A, F8.3, A, F8.3, A)') &
            "Simulation range R from ", R_min, " fm To ", R_max, " fm, dR = ", dr, " fm"

      write(30,'(A, F8.3, A, F8.3, A, F8.3, A, F8.3, A)') &
            "Potential parameters: V0= ", V0, "(MeV), a= ", A0, "(fm), r0= ", R0, "(fm), R= ",&
                                       R0*(Apro**(1.0/3.0) + Atar**(1.0/3.0)), "(fm)"

      write(30,'(A, F8.4, A)') "Coulomb barrier position :", R_barrier, "fm"

      write(30,'(A, F8.4, A)') "Coulomb barrier energy   :", V_barrier, "MeV"

      write(30,'(A, F8.4, A)') "Coulomb barrier curv     :", curv, "MeV"

      write(30,'(A, F8.4, A)') "Coulomb bottom position  :", R_bottom, "fm"

      write(30,'(A, F10.4, A)') "Coulomb bottom energy    :", V_bottom, "MeV"

      if (Ntar /= 0) then
          if (IVIBROTT == 0) then
              write(*,'(A, F6.3, A, F6.3, A, I0, A, I0)') "Phonon Excitation in the targ.: beta=",&
                           BetaT, ", omega=", OmegaT, " (MeV), Lambda=", LambdaT, ", Nph=", NphononT
              BetaTn = BetaT
              !Write(*,*)' Different beta_N from beta_C for this mode(n/y)?'
              !Read(*,*)ans
              ans = 'n'
              IF(ans .Eq. 'Y' .Or. ans .Eq. 'y') Then
                Write(*,*)'beta_N=?'
                Read(*,*)BetaTn
              End If
              write(*,'(A, F6.3, A, F6.3, A, F6.3, A)') "Phonon Excitation in the targ.: beta_N=",&
                                             BetaTn, ", beta_C=", BetaT, ", r0=", R0T, " (fm),"
              write(*,'(A, F6.3, A, I0, A, I0)') "omega=", OmegaT, " (MeV), Lambda=", LambdaT, ", Nph=", NphononT
              write(30,'(A, F6.3, A, F6.3, A, F6.3, A)') "Phonon Excitation in the targ.: beta_N=",&
                                                   BetaTn, ", beta_C=", BetaT, ", r0=", R0T, " (fm),"
              write(30,'(A, F6.3, A, I0, A, I0)') "  omega=", OmegaT, " (MeV), Lambda=", LambdaT, ", Nph=", NphononT
          else if (IVIBROTT == 1) then
              write(*,'(A, F6.3, A, F6.3, A, F6.3, A)') "Rotational Excitation in the targ.: beta2=",&
                                       Beta2T, ", beta4=", Beta4T, ", r0=", R0T, " (fm),"
              write(*,'(A, F6.3, A, I0)') " E2=", E2T, " (MeV), Nrot=", NrotT
              write(30,'(A, F6.3, A, F6.3, A, F6.3, A)') "Rotational Excitation in the targ.: beta2=",&
                                                 Beta2T, ", beta4=", Beta4T, ", r0=", R0T, " (fm),"
              write(30,'(A, F6.3, A, I0)') "  E2=", E2T, " (MeV), Nrot=", NrotT
          end if
      end if

      if (NphononT2 /= 0) then
          write(*,'(A, F6.3, A, F6.3, A, I0, A, I0)') "Phonon Excitation in the targ.: beta=",&
           BetaT2, ", omega=", OmegaT2, " (MeV), Lambda=", LambdaT2, ", Nph=", NphononT2
          BetaT2n = BetaT2
          Write(*,*)' Different beta_N from beta_C for this mode(n/y)?'
          Read(*,*)ans
          IF(ans .Eq. 'Y' .Or. ans .Eq. 'y') Then
            Write(*,*)'beta_N=?'
            Read(*,*)BetaT2n
          End If
          write(*,'(A, F6.3, A, F6.3, A, F6.3, A)') "Phonon Excitation in the targ.: beta_N=",&
                                   BetaT2n, ", beta_C=", BetaT2, ", r0=", R0T, " (fm),"
          write(*,'(A, F6.3, A, I0, A, I0)') "omega=", OmegaT2, " (MeV), Lambda=", LambdaT2, ", Nph=", NphononT2
          write(30,'(A, F6.3, A, F6.3, A, F6.3, A)') "Phonon Excitation in the targ.: beta_N=",&
                                     BetaT2n, ", beta_C=", BetaT2, ", r0=", R0T, " (fm),"
          write(30,'(A, F6.3, A, I0, A, I0)') "  omega=", OmegaT2, " (MeV), Lambda=", LambdaT2, ", Nph=", NphononT2
          Call Mutual
      End if

      print *, "------------------------------"
      print *, "Mode of excitation for projectile"

      if (Npro /= 0) then
          if (IVIBROTP == 0) then
              write(*,'(A, F6.3, A, F6.3, A, I0, A, I0)') "Phonon Excitation in the proj.: beta=",&
                         BetaP, ", omega=", OmegaP, " (MeV), Lambda=", LambdaP, ", Nph=", NphononP
              BetaPn = BetaP
              Write(*,*)' Different beta_N from beta_C for this mode(n/y)?'
              Read(*,*)ans
              IF(ans .Eq. 'Y' .Or. ans .Eq. 'y') Then
                Write(*,*)'beta_N=?'
                Read(*,*)BetaPn
              End If
              write(*,'(A, F6.3, A, F6.3, A, F6.3, A)') "Phonon Excitation in the proj.: beta_N=",&
                       BetaPn, ", beta_C=", BetaP, ", r0=", R0P, " (fm),"
              write(*,'(A, F6.3, A, F6.3, A, I0)') "omega=", OmegaP, " (MeV), Lambda=", LambdaP, ", Nph=", NphononP
              write(30,'(A, F6.3, A, F6.3, A, F6.3, A)') "Phonon Excitation in the proj.: beta_N=",&
                                         BetaPn, ", beta_C=", BetaP, ", r0=", R0P, " (fm),"
              write(30,'(A, F6.3, A, F6.3, A, I0)') "  omega=", OmegaP, " (MeV), Lambda=", LambdaP, ", Nph=", NphononP
          else if (IVIBROTP == 1) then
              write(*,'(A, F6.3, A, F6.3, A, F6.3, A)') "Rotational Excitation in the proj.: beta2=",&
                                                 Beta2P, ", beta4=", Beta4P, ", r0=", R0P, " (fm),"
              write(*,'(A, F6.3, A, I0)') "   E2=", E2P, " (MeV), Nrot=", NrotP
              write(30,'(A, F6.3, A, F6.3, A, F6.3, A)') "Rotational Excitation in the proj.: beta2=",&
                                                 Beta2P, ", beta4=", Beta4P, ", r0=", R0P, " (fm),"
              write(30,'(A, F6.3, A, I0)') "  E2=", E2P, " (MeV), Nrot=", NrotP
          end if
      end if
      Call Anharmonicity
      Call grotation

      print *, "------------------------------"
      print *, "Transfer coupled mode"

      if (Ntrans /= 0) then
          write(*,'(A, F6.3, A, F6.3, A)') "Transfer channel: Strength= ",&
                                       Ftr, ", Q = ", Qtrans, " MeV"
          write(30,'(A, F6.3, A, F6.3, A)') "Transfer channel: Strength= ",&
                                       Ftr, ", Q = ", Qtrans, " MeV"
      end if
      print *, "------------------------------"
  End Subroutine Output_information_old



  Subroutine Initialize_CCFull()

    Implicit None
    Real(8) ::  step
    Integer ::  i
    Allocate(imutual(0:Nlevelmax,0:Nlevelmax))
    Allocate(Ev(Nlevelmax))
    Allocate(Evec(Nlevelmax,Nlevelmax))
    Allocate(betnahv(0:Nlevelmax,0:Nlevelmax))
    Allocate(betcahv(0:Nlevelmax,0:Nlevelmax))
    Allocate(omeahv(0:Nlevelmax))
    Allocate(betnahv2(0:Nlevelmax,0:Nlevelmax))
    Allocate(betcahv2(0:Nlevelmax,0:Nlevelmax))
    Allocate(omeahv2(0:Nlevelmax+1))
    Allocate(betnahvp(0:Nlevelmax,0:Nlevelmax))
    Allocate(betcahvp(0:Nlevelmax,0:Nlevelmax))
    Allocate(omeahvp(0:Nlevelmax))
    Allocate(bett(0:Nlevelmax,0:Nlevelmax))
    Allocate(erott(0:Nlevelmax+1))
    Allocate(betp(0:Nlevelmax,0:Nlevelmax))
    Allocate(erotp(0:Nlevelmax))
    Allocate(eps(0:Nlevelmax))

    Call Read_Input()

    Call PotShape(0, R_barrier, V_barrier, curv, R_bottom, V_bottom)
    R_iterat = int((R_max - R_bottom) / dr)
    R_min = R_bottom
    R_max = R_bottom + dr * R_iterat
    t = Hbar**2 / (2 * ReduceMass * dr**2)
    
    Nx = R_iterat
    Nx_cc = Nlevel * Nx
    Nx1 = Nx_cc + Nlevel
    Nx2 = Nx_cc + 2 * Nlevel
    Write(*,*)'Nx_cc = ', Nx_cc
    num_traning_set1 = 3
    num_traning_set2 = 3
    num_testing_set1 = 11
    num_testing_set2 = 11
    NumEmulator = num_traning_set1 * num_traning_set2
    Num_testing = num_testing_set1 * num_testing_set2
    Allocate(training_set1(num_traning_set1))
    Allocate(training_set2(num_traning_set2))  
    Allocate(testing_set1(num_testing_set1))
    Allocate(testing_set2(num_testing_set2))

    Do i = 1, num_traning_set1
      step = 0.05
      training_set1(i) = 0.25d0 + (i - 1) * step
    End Do
    training_set1(1) = 0.25
    training_set1(2) = 0.30
    training_set1(3) = 0.35
    !training_set1(4) = 0.20
    !training_set1(4) = 0.25
 
   

    Write(*, '(10F8.2)') (training_set1(i), i=1, size(training_set1))

    Do i = 1, num_traning_set2
      step = 0.025
      training_set2(i) = 0.00d0  + (i - 1) * step
    End Do
    training_set2(1) = -0.025
    training_set2(2) = 0.00
    training_set2(3) = 0.025

    Write(*, '(10F8.3)') (training_set2(i), i=1, size(training_set2))
    Write(*,*)'training number = ', num_traning_set1*num_traning_set2
    Do i = 1, num_testing_set1
      step = 0.01
      testing_set1(i) = 0.25d0 + (i - 1) * step
    End Do
    Write(*,'(A, F8.3, A, F8.3, A, F8.3)')'testing_set_1 from', testing_set1(1), 'to',&
                                           testing_set1(num_testing_set1),', step=',step

    Do i = 1, num_testing_set2
      step = 0.01
      testing_set2 (i) = -0.05d0 + (i - 1) * step
    End Do
    Write(*,'(A, F8.3, A, F8.3, A, F8.3)')'testing_set_2 from', testing_set2(1), 'to',&
                                           testing_set2(num_testing_set2),', step=', step
    Write(*,*)'Testnumber = ', num_testing_set1*num_testing_set2



    Allocate(CPOT(Nlevelmax,Nlevelmax,0:R_iterat+2))
    Allocate(CPOT0(Nlevelmax,Nlevelmax))
    Allocate(CPOTH(Nlevelmax,Nlevelmax))
    Allocate(wf(Nx2))
    Allocate(wf_base(Nx2, 1:NumEmulator+Nlevel))
    Allocate(wf_base2(Nx1, 1:NumEmulator+Nlevel))
    Allocate(wf_base3_dagger(1:NumEmulator+Nlevel, Nx1))
    Allocate(Hamilton(Nx1, Nx2))
    Allocate(Ec_psi(Nx2, Nx1))
    Allocate(Ec_psi0(Nx2))
    Allocate(Ec_psi0_2(Nx1))
    Allocate(Echannel(1:Nlevelmax))
    Allocate(cc_ec(NumEmulator+Nlevel, NumEmulator+Nlevel))
    Allocate(cc_ec_inv(NumEmulator+Nlevel, NumEmulator+Nlevel))
    Allocate(dd_ec(NumEmulator+Nlevel))
    Allocate(dd2_ec(NumEmulator+Nlevel))
    Allocate(a_tra(Nlevel))

    



  End Subroutine Initialize_CCFull


End Module ccfull_initialization_mod


Subroutine PotShape(l, rb, vb, curv_out, rmin, vmin)
  Use ccfull_initialization_mod
  
  Implicit None
  Integer :: n, i
  Integer, Intent(In) :: l
  Real(8), Intent(Out) :: rb, vb, curv_out, rmin, vmin
  Real(8) :: r, u0, u1, ra, rb_tmp, tolk, u, ddv00
  Real(8) :: dV_dr, V
  external dV_dr, V
  r = 50.5D0
  u0 = dV_dr(r, l)

  Do
    r = r - 1.0D0
    If (r < 0.0D0) Then
      rb = -5.0D0
      vb = -5.0D0
      curv_out = -1.0D0
      rmin = -1.0D0
      vmin = -1.0D0
      Return
    End If
    u1 = dV_dr(r, l)
    If (u0 * u1 <= 0.0D0) Exit
    u0 = u1
  End Do

  ra = r + 1.0D0
  rb_tmp = r
  tolk = 1.0D-6
  n = Int(Log10(Abs(rb_tmp - ra) / tolk) / Log10(2.0D0) + 0.5D0)

  Do i = 1, n
    r = (ra + rb_tmp) / 2.0D0
    u = dV_dr(r, l)
    If (u0 * u < 0.0D0) Then
      rb_tmp = r
    Else
      ra = r
    End If
  End Do

  rb = r
  ddv00 = (dV_dr(rb + 1.0D-5, l) - dV_dr(rb - 1.0D-5, l)) / (2.0D-5)

  If (ddv00 > 0.0D0) Then
    Print *, "Error: Second derivative positive at barrier top."
    Stop
  End If

  curv_out = Hbar * Sqrt(Abs(ddv00) / ReduceMass)
  vb = V(rb, l)

  ra = rb - 0.5D0
  u0 = dV_dr(ra, l)
  rb_tmp = 0.5D0
  u1 = dV_dr(rb_tmp, l)
  n = Int(Log10(Abs(rb - ra) / tolk) / Log10(2.0D0) + 0.5D0)

  Do i = 1, n
    r = (ra + rb_tmp) / 2.0D0
    u = dV_dr(r, l)
    If (u0 * u < 0.0D0) Then
      rb_tmp = r
    Else
      ra = r
    End If
  End Do

  rmin = r
  vmin = V(rmin, l)

End Subroutine PotShape

Real(8) Function Vn(r)
  Use ccfull_initialization_mod
  Real(8), Intent(In) :: r
  Real(8) :: R12
  R12 = R0 * (Apro**(1.0D0/3.0D0) + Atar**(1.0D0/3.0D0))
  Vn = -V0 / (1.0D0 + Exp((r - R12) / A0))
End Function Vn

Real(8) Function W(r)
  Use ccfull_initialization_mod
  Real(8), Intent(In) :: r
  Real(8) :: rw
  rw = 1.0D0 * (Apro**(1.0D0/3.0D0) + Atar**(1.0D0/3.0D0))
  W = 50.0D0 / (1.0D0 + Exp((r - rw) / 0.3D0))
End Function W

Real(8) Function dVn_dr(r)
  Use ccfull_initialization_mod
  Real(8), Intent(In) :: r
  Real(8) :: R12, exp_term
  R12 = R0 * (Apro**(1.0D0/3.0D0) + Atar**(1.0D0/3.0D0))
  exp_term = Exp((r - R12) / A0)
  dVn_dr = V0 / A0 * exp_term / (1.0D0 + exp_term)**2
End Function dVn_dr

Real(8) Function dVcent_dr(r, l)
  Use ccfull_initialization_mod
  Real(8), Intent(In) :: r
  Integer, Intent(In) :: l
  dVcent_dr = -2.0D0 * l * (l + 1.0D0) * Hbar**2 / (2.0D0 * ReduceMass * r**3)
End Function dVcent_dr

Real(8) Function dVc_dr(r)
  Use ccfull_initialization_mod
  Real(8), Intent(In) :: r
  dVc_dr = -Zpro * Ztar * Hbar / 137.0D0 / r**2
End Function dVc_dr

Real(8) Function dV_dr(r, l)
  Use ccfull_initialization_mod
  Real(8), Intent(In) :: r
  Integer, Intent(In) :: l
  Real(8)             :: dVc_dr, dVn_dr, dVcent_dr
  external dVc_dr, dVn_dr, dVcent_dr
  dV_dr = dVn_dr(r) + dVc_dr(r) + dVcent_dr(r, l)
End Function dV_dr

Real(8) Function Vc(r)
  Use ccfull_initialization_mod
  Real(8), Intent(In) :: r
  Vc = Zpro * Ztar * Hbar / 137.0D0 / r
End Function Vc

Real(8) Function Vr(r, l)
  Use ccfull_initialization_mod
  Real(8), Intent(In) :: r
  Integer, Intent(In) :: l
  Vr = h2m * l * (l + 1.0D0) / r**2
End Function Vr

Real(8) Function V(r, l)
  Use ccfull_initialization_mod
  Real(8), Intent(In) :: r
  Integer, Intent(In) :: l
  Real(8)             :: Vn, Vc, Vr
  external Vn, Vc, Vr
  V = Vn(r) + Vc(r) + Vr(r, l)
End Function V

Real(8) Function VnCC(r, Xt)
  Use ccfull_initialization_mod
  Real(8), Intent(In) :: r, Xt
  Real(8) :: R12
  R12 = R0 * (Apro**(1.0D0/3.0D0) + Atar**(1.0D0/3.0D0))
  VnCC = -V0 / (1.0D0 + Exp((r - R12 - Xt) / A0))
End Function VnCC

Real(8) Function Fct(r)
  Use ccfull_initialization_mod
  IMPLICIT NONE
  Real(8), Intent(In) :: r
  Real(8) ::  result, Lambda

  Lambda = LambdaT
  If (r > Rtar) Then
     result = 3.D0 / (2.D0*Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / r * (Rtar / r)**Lambda
  Else
     result = 3.D0 / (2.D0*Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / Rtar * (r / Rtar)**Lambda
  End If

  result = result * BetaT / Sqrt(4.D0 * Pi)

  Fct = result
  Return
End Function Fct

Real(8) Function Fct2(r)
  Use ccfull_initialization_mod
  Implicit None
  Real(8), Intent(In) :: r
  Real(8) :: result
  Integer, Parameter :: Lambda = 2

  If (r > Rtar) Then
     result = 3.D0 / (2.D0*Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / r * (Rtar / r) ** Lambda
  Else
     result = 3.D0 / (2.D0*Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / Rtar * (r / Rtar) ** Lambda
  End If

  Fct2 = result
  Return
End Function Fct2

Real(8) Function Fct3(r)
  Use ccfull_initialization_mod
  Implicit None
  Real(8), Intent(In) :: r
  Real(8) :: result
  Integer, Parameter :: Lambda = 3

  If (r > Rtar) Then
     result = 3.D0 / (2.D0*Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / r * (Rtar / r) ** Lambda
  Else
     result = 3.D0 / (2.D0*Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / Rtar * (r / Rtar) ** Lambda
  End If

  result = result * (Beta4T + 7.D0 * Beta2T**2 / 7.D0 / Sqrt(Pi))

  Fct3 = result
  Return
End Function Fct3

Real(8) Function Fct4(r)
  Use ccfull_initialization_mod
  Implicit None
  Real(8), Intent(In) :: r
  Real(8) :: result
  Integer, Parameter :: Lambda = 4

  If (r > Rtar) Then
     result = 3.D0 / (2.D0*Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / r * (Rtar / r) ** Lambda
  Else
     result = 3.D0 / (2.D0*Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / Rtar * (r / Rtar) ** Lambda
  End If

  result = result * (Beta4T + 9.D0 * Beta2T**2 / 7.D0 / Sqrt(Pi))

  Fct4 = result
  Return
End Function Fct4

Real(8) Function Fcp(r)
  Use ccfull_initialization_mod
  Implicit None
  Real(8), Intent(In) :: r
  Real(8) :: result, Lambda

  Lambda = LambdaP

  If (r > Rpro) Then
     result = 3.D0 / (2.D0*Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / r * (Rpro / r) ** Lambda
  Else
     result = 3.D0 / (2.D0*Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / Rpro * (r / Rpro) ** Lambda
  End If

  result = result * BetaP / SQRT(4.D0 * Pi)

  Fcp = result
  Return
End Function Fcp

Real(8) Function Fcp2(r)
  Use ccfull_initialization_mod
  Implicit None
  Real(8), Intent(In) :: r
  Real(8) :: result
  Integer :: Lambda

  Lambda = 2

  If (r > Rpro) Then
     result = 3.D0 / (2.D0 * Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / r * (Rpro / r)**Lambda
  Else
     result = 3.D0 / (2.D0 * Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / Rpro * (r / Rpro)**Lambda
  End If

  Fcp2 = result
  Return
End Function Fcp2

Real(8) Function Fcp3(r)
  Use ccfull_initialization_mod
  Implicit None
  Real(8), Intent(In) :: r
  Real(8) :: result
  Integer :: Lambda

  Lambda = 3
  result = 0.D0

  If (r > Rpro) Then
     result = 3.D0 / (2.D0 * Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / r * (Rpro / r)**Lambda
  Else
     result = 3.D0 / (2.D0 * Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / Rpro * (r / Rpro)**Lambda
  End If

  result = result * (Beta4P + 7.D0 * Beta2P**2 / 7.D0 / Sqrt(Pi))

  Fcp3 = result
  Return
End Function Fcp3

Real(8) Function Fcp4(r)
  Use ccfull_initialization_mod
  Implicit None
  Real(8), Intent(In) :: r
  Real(8) :: result
  Integer :: Lambda

  Lambda = 4
  result = 0.D0

  If (r > Rpro) Then
     result = 3.D0 / (2.D0 * Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / r * (Rpro / r)**Lambda
  Else
     result = 3.D0 / (2.D0 * Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / Rpro * (r / Rpro)**Lambda
  End If

  result = result * (Beta4P + 9.D0 * Beta2P**2 / 7.D0 / Sqrt(pi))

  Fcp4 = result
  Return
End Function Fcp4

Real(8) Function Fctt(r)
  Use ccfull_initialization_mod
  Implicit None
  Real(8), Intent(In) :: r
  Real(8) :: result
  Integer :: Lambda

  Lambda = LambdaT2
  result = 0.D0

  If (r > Rtar) Then
     result = 3.D0 / (2.D0 * Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / r * (Rtar / r)**Lambda
  Else
     result = 3.D0 / (2.D0 * Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / Rtar * (r / Rtar)**Lambda
  End If

  result = result * BetaT2 / SQRT(4.D0 * Pi)

  Fctt = result
  Return
End Function Fctt

Real(8) Function Fct2v(r)
  Use ccfull_initialization_mod
  Implicit None
  Real(8), Intent(In) :: r
  Real(8) :: result
  Integer :: Lambda

  Lambda = 2
  result = 0.D0

  If (r > Rtar) Then
     result = 3.D0 / (2.D0 * Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / r * (Rtar / r)**Lambda
  Else
     result = 3.D0 / (2.D0 * Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / Rtar * (r / Rtar)**Lambda
  End If

  Fct2v = result
  Return
End Function Fct2v

Real(8) Function Fcp2v(r)
  Use ccfull_initialization_mod
  Implicit None
  Real(8), Intent(In) :: r
  Real(8) :: result
  Integer :: Lambda

  Lambda = 2
  result = 0.D0

  If (r > Rpro) Then
     result = 3.D0 / (2.D0 * Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / r * (Rpro / r)**Lambda
  Else
     result = 3.D0 / (2.D0 * Lambda + 1.D0) * Zpro * Ztar / 137.D0 * Hbar / Rpro * (r / Rpro)**Lambda
  End If

  Fcp2v = result
  Return
End Function Fcp2v

Real(8) Function Ftrans(r)
  Use ccfull_initialization_mod
  Implicit None
  Real(8), Intent(In) :: r
  Real(8) :: result, dVn_dr
  external dVn_dr

  result = 0.D0
  result = Ftr * dVn_dr(r)

  Ftrans = result
  Return
End Function Ftrans

Real(8) Function CG(j1, m1, j2, m2, j3, m3)
    ! Clebsch-Gordan coefficient <j1 m1 j2 m2 | j3 m3>
    implicit real*8 (a-h, o-z)
    external fact

    if (m1 + m2 .ne. m3) then
        cg = 0.d0
        return
    end if

    if (j3 .lt. abs(j1 - j2)) then
        cg = 0.d0
        return
    end if

    if (j3 .gt. j1 + j2) then
        cg = 0.d0
        return
    end if

    ka = j1 + j2 - j3
    kb = j3 + j1 - j2
    kc = j2 + j3 - j1
    kd = j1 + j2 + j3 + 1

    del = sqrt(fact(ka) * fact(kb) * fact(kc) / fact(kd))

    cg = 0.d0
    do n = 0, max(j1 + j2 - j3, j1 - m1, j2 + m2)
        ka1 = j1 + j2 - j3 - n
        if (ka1 .lt. 0.d0) cycle

        ka2 = j3 - j2 + m1 + n
        if (ka2 .lt. 0.d0) cycle

        ka3 = j3 - j1 - m2 + n
        if (ka3 .lt. 0.d0) cycle

        ka4 = j1 - m1 - n
        if (ka4 .lt. 0.d0) cycle

        ka5 = j2 + m2 - n
        if (ka5 .lt. 0.d0) cycle

        cg = cg + (-1.d0)**n / &
                 (fact(n) * fact(ka1) * fact(ka2) * fact(ka3) * fact(ka4) * fact(ka5))
    end do

    cg = cg * sqrt(fact(j1 + m1) * fact(j1 - m1))
    cg = cg * sqrt(fact(j2 + m2) * fact(j2 - m2))
    cg = cg * sqrt(fact(j3 + m3) * fact(j3 - m3))

    cg = cg * sqrt(2.d0 * j3 + 1.d0) * del

    return
End Function CG

Subroutine mdiag(a,n,d,v)

	implicit real*8(a-h,o-z)
	parameter(nmax=100)
  parameter (nlevelmax=30)
  dimension a(nlevelmax,nlevelmax)
	dimension b(nmax),z(nmax),d(nlevelmax),v(nlevelmax,nlevelmax)

	do 12 ip=1,n
	do 11 iq=1,n
	v(ip,iq)=0.d0
 11	continue
	v(ip,ip)=1.d0
 12	continue

	do 13 ip=1,n
	b(ip)=a(ip,ip)
	d(ip)=b(ip)
	z(ip)=0.d0
 13	continue

	nrot=0

	do 24 i=1,50
	sm=0.d0

	do 15 ip=1,n-1
	do 14 iq=ip+1,n
	sm=sm+abs(a(ip,iq))
 14	continue
 15     continue

	if(sm.eq.0.d0) return

	if(i.lt.4) then
	tresh=0.2d0*sm/n**2
	else
	tresh=0.d0
	endif

	do 22 ip=1,n-1
	do 21 iq=ip+1,n
	g=100.d0*abs(a(ip,iq))
	
	if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq)))) then
	   a(ip,iq)=0.d0
	elseif(abs(a(ip,iq)).gt.tresh) then
	   h=d(iq)-d(ip)

	   if(abs(h)+g.eq.abs(h)) then
	     t=a(ip,iq)/h
	   else
	     theta=0.5d0*h/a(ip,iq)
	     t=1.d0/(abs(theta)+sqrt(1.d0+theta**2))
	     if(theta.lt.0.d0) t=-t
	   endif

	c=1.d0/sqrt(1.d0+t**2)
	s=t*c
	tau=s/(1.d0+c)
	h=t*a(ip,iq)
	z(ip)=z(ip)-h
	z(iq)=z(iq)+h
	d(ip)=d(ip)-h
	d(iq)=d(iq)+h
	a(ip,iq)=0.d0

	do 16 j=1,ip-1
	g=a(j,ip)
	h=a(j,iq)
	a(j,ip)=g-s*(h+g*tau)
	a(j,iq)=h+s*(g-h*tau)
 16	continue

	do 17 j=ip+1,iq-1
	g=a(ip,j)
	h=a(j,iq)
	a(ip,j)=g-s*(h+g*tau)
	a(j,iq)=h+s*(g-h*tau)
 17	continue

	do 18 j=iq+1,n
	g=a(ip,j)
	h=a(iq,j)
	a(ip,j)=g-s*(h+g*tau)
	a(iq,j)=h+s*(g-h*tau)
 18	continue

	do 19 j=1,n
	g=v(j,ip)
	h=v(j,iq)
	v(j,ip)=g-s*(h+g*tau)
	v(j,iq)=h+s*(g-h*tau)
 19	continue
	
	nrot=nrot+1
	endif

 21	continue
 22	continue

	do 23 ip=1,n
	b(ip)=b(ip)+z(ip)
	d(ip)=b(ip)
	z(ip)=0.d0

 23	continue
 24	continue

	return

End Subroutine  mdiag

SUBROUTINE matinv(nmax, c, d)
  IMPLICIT NONE
  
  ! Argument declarations
  INTEGER, INTENT(IN) :: nmax
  COMPLEX(8), DIMENSION(30,30), INTENT(INOUT) :: c, d
  
  ! Parameter and local variable declarations
  COMPLEX(8) :: u(30), v(30), t, a, b
  REAL(8) :: deter
  INTEGER :: m, n, j, k, l
  
  ! Initialize
  deter = 1.0_8
  d = (0.0_8, 0.0_8)
  DO m = 1, nmax
      d(m, m) = (1.0_8, 0.0_8)
  END DO
  
  ! Main inversion loop
  DO n = 1, nmax
      t = c(n, n)
      
      ! Check for zero pivot element
      IF (ABS(t) < 1.0E-10_8) THEN
          j = n
          DO
              j = j + 1
              IF (j > nmax) THEN
                  PRINT *, 'matrix not invertible'
                  RETURN
              END IF
              
              t = c(n, j)
              deter = -deter
              IF (ABS(t) > 1.0E-10_8) EXIT
          END DO
          
          ! Swap rows
          DO k = 1, nmax
              u(k) = c(n, k)
              v(k) = d(n, k)
              c(n, k) = c(j, k)
              d(n, k) = d(j, k)
              c(j, k) = u(k)
              d(j, k) = v(k)
          END DO
      END IF
      
      ! Eliminate column
      DO k = 1, nmax
          IF (k == n) CYCLE
          
          a = c(k, n) / c(n, n)
          DO l = 1, nmax
              c(k, l) = c(k, l) - a * c(n, l)
              d(k, l) = d(k, l) - a * d(n, l)
          END DO
      END DO
      
      ! Normalize row
      b = c(n, n)
      deter = b * deter
      DO m = 1, nmax
          c(n, m) = c(n, m) / b
          d(n, m) = d(n, m) / b
      END DO
  END DO
  
END SUBROUTINE matinv

subroutine matinv_lapack(n, A, Ainv)
  implicit none
  integer, intent(in) :: n
  complex*16, intent(in) :: A(n,n)
  complex*16, intent(out) :: Ainv(n,n)

  ! 局部变量
  integer :: info, lwork
  integer, allocatable :: ipiv(:)
  complex*16, allocatable :: work(:)
  complex*16 :: tmp(n,n)

  ! 拷贝 A 到临时变量
  tmp = A
  Ainv = (0.0d0, 0.0d0)

  allocate(ipiv(n))

  ! LU 分解：tmp 会被覆盖为 LU
  call zgetrf(n, n, tmp, n, ipiv, info)
  if (info /= 0) then
    print *, "Error in zgetrf: info = ", info
    stop
  end if

  ! 工作空间大小（可调优）
  lwork = n * n
  allocate(work(lwork))

  ! 计算逆矩阵
  call zgetri(n, tmp, n, ipiv, work, lwork, info)
  if (info /= 0) then
    print *, "Error in zgetri: info = ", info
    stop
  end if

  Ainv = tmp

  deallocate(ipiv)
  deallocate(work)
end subroutine matinv_lapack

Real(8) Function fact(n)
    Implicit none
    Integer n, i

    If (n .lt. 0) Then
        fact = 0.d0
        Return
    End If

    If (n .eq. 0) Then
        fact = 1.d0
        Return
    End If

    fact = 1.d0
    Do i = 1, n
        fact = fact * i * 1.d0
    End Do

    Return
End Function fact

Subroutine grotation()
    Use ccfull_initialization_mod
    Implicit None
    Integer :: it, jt, ip, jp, i
    Character(len=1) :: ans
    Real(8) :: erott0, erotp0

    Do it = 0, Ntar
      Do jt = 0, Ntar
        bett(it,jt) = 0.0D0
        If (it == jt - 1 .Or. it == jt + 1 .Or. it == jt) bett(it,jt) = Beta2T
      End Do
      erott(it) = 2.0D0 * it * (2.0D0 * it + 1.0D0) / 6.0D0 * E2T
    End Do

    If (IVIBROTT == 1 .And. Ntar > 1) Then
      !Print *, ''
      !Print *, 'Generalised E2 couplings for the target rotor (n/y)?'
      !Read(*,'(A)') ans
      ans = 'n'
      If (ans == 'y' .Or. ans == 'Y') Then
        Print *, 'Beta2_targ=', Beta2T
        Do it = 1, Ntar - 1
          jt = it + 1
          Print *, 'Modify beta2 for transition from', 2*it, '+ to', 2*jt, '+ ? (n/y)'
          Read(*,'(A)') ans
          If (ans == 'y' .Or. ans == 'Y') Then
            Print *, 'BETA2=?'
            Read(*,*) bett(it,jt)
            bett(jt,it) = bett(it,jt)
            Print *, 'Beta2 set to: ', bett(it,jt)
          End If
        End Do

        Print *, 'Reorientation terms:'
        Do it = 1, Ntar
          jt = it
          Print *, 'Modify beta2 for transition from', 2*it, '+ to', 2*jt, '+ ? (n/y)'
          Read(*,'(A)') ans
          If (ans == 'y' .Or. ans == 'Y') Then
            Print *, 'BETA2=?'
            Read(*,*) bett(it,jt)
            bett(jt,it) = bett(it,jt)
          End If
        End Do

        Print *, 'Excitation energies:'
        Do i = 2, Ntar
          Print *, 'Energy of the', 2*i, '+ state for pure rotor =', erott(i)
          Print *, 'Modify this energy? (n/y)'
          Read(*,'(A)') ans
          If (ans == 'y' .Or. ans == 'Y') Then
            erott0 = erott(i)
            Print *, 'Energy=?'
            Read(*,*) erott(i)
            Print *, 'Modified from', erott0, ' to ', erott(i)
          End If
        End Do
      End If
    End If

    Do ip = 0, Npro
      Do jp = 0, Npro
        betp(ip,jp) = 0.0D0
        If (ip == jp - 1 .Or. ip == jp + 1 .Or. ip == jp) betp(ip,jp) = Beta2P
      End Do
      erotp(ip) = 2.0D0 * ip * (2.0D0 * ip + 1.0D0) / 6.0D0 * E2P
    End Do

    If (IVIBROTP == 1 .And. Npro > 1) Then
      !Print *, 'Generalised E2 couplings for the projectile rotor (n/y)?'
      !Read(*,'(A)') ans
      ans = 'n'
      If (ans == 'y' .Or. ans == 'Y') Then
        Print *, 'Beta2_proj=', Beta2P
        Do ip = 1, Npro - 1
          jp = ip + 1
          Print *, 'Modify beta2 for transition from', 2*ip, '+ to', 2*jp, '+ ? (n/y)'
          Read(*,'(A)') ans
          If (ans == 'y' .Or. ans == 'Y') Then
            Print *, 'BETA2=?'
            Read(*,*) betp(ip,jp)
            betp(jp,ip) = betp(ip,jp)
          End If
        End Do

        Print *, 'Reorientation terms:'
        Do ip = 1, Npro
          jp = ip
          Print *, 'Modify beta2 for transition from', 2*ip, '+ to', 2*jp, '+ ? (n/y)'
          Read(*,'(A)') ans
          If (ans == 'y' .Or. ans == 'Y') Then
            Print *, 'BETA2=?'
            Read(*,*) betp(ip,jp)
            betp(jp,ip) = betp(ip,jp)
          End If
        End Do

        Print *, 'Excitation energies:'
        Do i = 2, Npro
          Print *, 'Energy of the', 2*i, '+ state for pure rotor =', erotp(i)
          Print *, 'Modify this energy? (n/y)'
          Read(*,'(A)') ans
          If (ans == 'y' .Or. ans == 'Y') Then
            erotp0 = erotp(i)
            Print *, 'Energy=?'
            Read(*,*) erotp(i)
            Print *, 'Modified from', erotp0, ' to ', erotp(i)
          End If
        End Do
      End If
    End If

End Subroutine grotation

Subroutine Mutual()
    Use ccfull_initialization_mod
    Implicit None
    Integer :: i, j
    Character(len=1) :: ans, sign1, sign2
    Integer :: imut

    Do i = 0, Nlevelmax
      Do j = 0, Nlevelmax
        imutual(i,j) = 0
      End Do
    End Do

    Print *, ''
    Print *, 'Mutual excitations in the *target* nucleus'
    Print *, 'Include the mutual excitations (y/n)?'
    Read(*,'(A)') ans

    If (ans == 'n' .Or. ans == 'N') Then
      imut = 0
      Nlevel = (Ntar + NphononT2 + 1) * (Npro + 1)
      Do i = 0, Ntar
        Do j = 0, NphononT2
          If (i == 0 .Or. j == 0) imutual(i,j) = 1
        End Do
      End Do
      Return
    End If

    Print *, 'All the possible mutual excitation channels (n/y)?'
    Read(*,'(A)') ans

    If (ans == 'y' .Or. ans == 'Y') Then
      imut = 1
      Nlevel = (Ntar + 1) * (NphononT2 + 1) * (Npro + 1)
      Do i = 0, Ntar
        Do j = 0, NphononT2
          imutual(i,j) = 1
        End Do
      End Do
      Return
    End If

    imut = 2
    If (Mod(LambdaT, 2) == 0) Then
      sign1 = '+'
    Else
      sign1 = '-'
    End If
    If (Mod(LambdaT2, 2) == 0) Then
      sign2 = '+'
    Else
      sign2 = '-'
    End If

    Do i = 0, Ntar
      Do j = 0, NphononT2
        If (i == 0 .Or. j == 0) Then
          imutual(i,j) = 1
        Else
          Print *, 'Include (', LambdaT, sign1, '^', i, ',', LambdaT2, sign2, '^', j, ') state ? (y/n)'
          Read(*,'(A)') ans
          If (ans == 'n' .Or. ans == 'N') Then
            imutual(i,j) = 0
          Else
            imutual(i,j) = 1
          End If
        End If
      End Do
    End Do

    Print *, 'Excited states in the target to be included:'
    Nlevel = 0

    Do i = 0, Ntar
      Do j = 0, NphononT2
        If (imutual(i,j) == 0) Cycle
        Print *, '(', LambdaT, sign1, '^', i, ',', LambdaT2, sign2, '^', j, ') state'
        Nlevel = Nlevel + 1
      End Do
    End Do

    Nlevel = Nlevel * (Npro + 1)

End Subroutine Mutual

Subroutine Anharmonicity()
  Use ccfull_initialization_mod
  Implicit None

  Character(1) :: ans, sign
  Integer :: it, jt, i

  ! Target: AHV Couplings (1st mode)
  Do it = 0, Ntar
    Do jt = 0, Ntar
      betnahv(it, jt) = 0.0
      betcahv(it, jt) = 0.0

      If (it == jt - 1) Then
        betnahv(it, jt) = BetaTn * Sqrt(Real(jt))
        betcahv(it, jt) = BetaT  * Sqrt(Real(jt))
      Else If (jt == it - 1) Then
        betnahv(it, jt) = BetaTn * Sqrt(Real(it))
        betcahv(it, jt) = BetaT  * Sqrt(Real(it))
      End If
    End Do
    omeahv(it) = it * OmegaT
  End Do

  If (IVIBROTT == 0 .And. Ntar > 1) Then
    Write(*, '(A)') 'AHV couplings for the first mode in the target phonon (n/y)?'
    Read(*, '(A1)') ans
    If (ans == 'Y' .Or. ans == 'y') Then
      Write(*, '(A)') '**** AHV Couplings in the target (the 1st mode)'
      If ((-1.)**LambdaT ==  1) sign = '+'
      If ((-1.)**LambdaT == -1) sign = '-'

      Do it = 0, Ntar
        Do jt = 0, Ntar
          If (it > jt)Cycle
          If ((it == 0 .And. jt <= 1)) Cycle
          Write(*, *) ' '
          Write(*, '(A,I2,2A,I2,A,I2,2A,I2,A)') 'Transition from the', LambdaT, sign, '^', it, &
                                                      ' to the', LambdaT, sign, '^', jt, 'state:'
          Write(*, '(A,2G15.5)') '   beta_N and beta_C in the HO limit=', betnahv(it, jt), betcahv(it, jt)
          Write(*, '(2A)') '    Modify these beta_N a/o beta_C (n/y)?'
          Read(*, '(A1)') ans
          If (ans == 'Y' .Or. ans == 'y') Then
            Write(*, *) '   beta_N and beta_C =?'
            Read(*, *) betnahv(it, jt), betcahv(it, jt)
            Write(*, '(G15.5, G15.5)') betnahv(it, jt), betcahv(it, jt)
          End If
        End Do
      End Do

      Write(*, *) 'Excitation energy for the first mode:'
      Do i = 2, Ntar
        Write(*, '(A,I2,A,G15.5)') '    Energy of the', i, '-phonon state in the HO=', i * OmegaT
        Write(*, *) '   Modify this energy(n/y)?'
        Read(*, '(A1)') ans
        If (ans == 'Y' .Or. ans == 'y') Then
          Write(*, *) '   Energy=?'
          Read(*, *) omeahv(i)
          Write(*, '(A,I2,A,G15.5,A,G12.5,A)') 'Energy of the', i, '-phonon state=', omeahv(i), '(HO: ', i * OmegaT, ')'
        End If
      End Do
    End If
  End If


  ! Target: AHV Couplings (2st mode)
  Do it = 0, NphononT2
    Do jt = 0, NphononT2
      betnahv2(it, jt) = 0.0
      betcahv2(it, jt) = 0.0

      If (it == jt - 1) Then
        betnahv2(it, jt) = BetaT2n * Sqrt(Real(jt))
        betcahv2(it, jt) = BetaT2  * Sqrt(Real(jt))
      Else If (jt == it - 1) Then
        betnahv(it, jt) = BetaT2n * Sqrt(Real(it))
        betcahv(it, jt) = BetaT2  * Sqrt(Real(it))
      End If
    End Do
    omeahv2(it) = it * OmegaT2
  End Do

  If (NphononT2 > 1) Then
    Write(*, '(A)') 'AHV couplings for the second mode in the target phonon (n/y)?'
    Read(*, '(A1)') ans
    If (ans == 'Y' .Or. ans == 'y') Then
      Write(*, '(A)') '**** AHV Couplings in the target (the 1st mode)'
      If ((-1.)**LambdaT2 ==  1) sign = '+'
      If ((-1.)**LambdaT2 == -1) sign = '-'

      Do it = 0, NphononT2
        Do jt = 0, NphononT2
          If (it > jt)Cycle
          If ((it == 0 .And. jt <= 1)) Cycle
          Write(*, *) ' '
          Write(*, '(A,I2,2A,I2,A,I2,2A,I2,A)') 'Transition from the', LambdaT2, sign, '^', it, &
                                                    ' to the', LambdaT2, sign, '^', jt, 'state:'
          Write(*, '(A,2G15.5)') '   beta_N and beta_C in the HO limit=', betnahv2(it, jt), betcahv2(it, jt)
          Write(*, '(2A)') '    Modify these beta_N a/o beta_C (n/y)?'
          Read(*, '(A1)') ans
          If (ans == 'Y' .Or. ans == 'y') Then
            Write(*, *) '   beta_N and beta_C =?'
            Read(*, *) betnahv2(it, jt), betcahv2(it, jt)
            Write(*, '(G15.5, G15.5)') betnahv2(it, jt), betcahv2(it, jt)
          End If
        End Do
      End Do

      Write(*, *) 'Excitation energy for the second mode:'
      Do i = 2, Ntar
        Write(*, '(A,I2,A,G15.5)') '    Energy of the', i, '-phonon state in the HO=', i * OmegaT2
        Write(*, *) '   Modify this energy(n/y)?'
        Read(*, '(A1)') ans
        If (ans == 'Y' .Or. ans == 'y') Then
          Write(*, *) '   Energy=?'
          Read(*, *) omeahv2(i)
          Write(*, '(A,I2,A,G15.5,A,G12.5,A)') 'Energy of the', i, '-phonon state=', omeahv2(i), '(HO: ', i * OmegaT2, ')'
        End If
      End Do
    End If
  End If


  Do it = 0, Npro
    Do jt = 0, Npro
      betnahvp(it, jt) = 0.0
      betcahvp(it, jt) = 0.0

      If (it == jt - 1) Then
        betnahvp(it, jt) = BetaPn * Sqrt(Real(jt))
        betcahvp(it, jt) = BetaP  * Sqrt(Real(jt))
      Else If (jt == it - 1) Then
        betnahvp(it, jt) = BetaPn * Sqrt(Real(it))
        betcahvp(it, jt) = BetaP  * Sqrt(Real(it))
      End If
    End Do
    omeahvp(it) = it * OmegaP
  End Do

  If (IVIBROTP == 0 .And. Npro > 1) Then
      Write(*, '(A)') 'AHV couplings for the first mode in the projectile phonon (n/y)?'
      Read(*, '(A1)') ans
      If (ans == 'Y' .Or. ans == 'y') Then
        Write(*, '(A)') '**** AHV Couplings in the projectile'
        If ((-1.)**LambdaP ==  1) sign = '+'
        If ((-1.)**LambdaP == -1) sign = '-'

        Do it = 0, Npro
          Do jt = 0, Npro
            If (it > jt)Cycle
            If ((it == 0 .And. jt <= 1)) Cycle
            Write(*, *) ' '
            Write(*, '(A,I2,2A,I2,A,I2,2A,I2,A)') 'Transition from the', LambdaP, sign, '^', it, &
                                                        ' to the', LambdaP, sign, '^', jt, 'state:'
            Write(*, '(A,2G15.5)') '   beta_N and beta_C in the HO limit=', betnahvp(it, jt), betcahvp(it, jt)
            Write(*, '(2A)') '    Modify these beta_N a/o beta_C (n/y)?'
            Read(*, '(A1)') ans
            If (ans == 'Y' .Or. ans == 'y') Then
              Write(*, *) '   beta_N and beta_C =?'
              Read(*, *) betnahvp(it, jt), betcahvp(it, jt)
              Write(*, '(G15.5, G15.5)') betnahvp(it, jt), betcahvp(it, jt)
            End If
          End Do
        End Do

        Write(*, *) 'Excitation energy for the first mode:'
        Do i = 2, Npro
          Write(*, '(A,I2,A,G15.5)') '    Energy of the', i, '-phonon state in the HO=', i * OmegaP
          Write(*, *) '   Modify this energy(n/y)?'
          Read(*, '(A1)') ans
          If (ans == 'Y' .Or. ans == 'y') Then
            Write(*, *) '   Energy=?'
            Read(*, *) omeahvp(i)
            Write(*, '(A,I2,A,G15.5,A,G12.5,A)') 'Energy of the', i, '-phonon state=', omeahvp(i), '(HO: ', i * OmegaP, ')'
          End If
        End Do
      End If
    End If

End Subroutine Anharmonicity

Subroutine coupled_matrix0()
  Use ccfull_initialization_mod
  implicit none
  real(8) :: A(Nlevelmax, Nlevelmax)
  Real(8) :: C, CG
  integer :: i, j, ip, it, it2, jp, jt, jt2
  external CG
  Call Output_information()
  !initialization of A matrix
  Do i = 1, Nlevelmax
    Do j = 1, Nlevelmax
     A(i,j) = 0.0d0
    End Do
  End Do

  If((Nlevel - Ntrans) == 1)Return
  i = 0
  Do ip = 0, Npro
    Do it = 0, Ntar
      Do it2 = 0, NphononT2

        if (NphononT2 .Ne. 0)Then
          if(imutual(it, it2) == 0)Cycle
        End if

        i = i + 1 
        j = 0
        Do jp = 0, Npro
          Do jt = 0, Ntar
            Do jt2 = 0, NphononT2

              if (NphononT2 .Ne. 0)Then
                if(imutual(jt, jt2) == 0)Cycle
              End if

              j = j + 1 

              If (i > j)Then
                A(i,j) = A(j,i)
                Cycle
              End If

              C = 0.0d0

              if (ip == jp .and. it2 == jt2) then
                  if (IVIBROTT == 0) then
                      C = Rtar * betnahv(it, jt) / sqrt(4.0d0 * pi)

                  else
                      C = Rtar * bett(it, jt) * sqrt((2 * 2 * it + 1) * 5.0d0 * (2 * 2 * jt + 1) / (4.0d0 * pi)) &
                          * CG(2 * it, 0, 2, 0, 2 * jt, 0)**2 / (2.0d0 * 2 * jt + 1)

                      C = C + Rtar * Beta4T * sqrt((2 * 2 * it + 1) * 9.0d0 * (2 * 2 * jt + 1) / (4.0d0 * pi)) &
                          * CG(2 * it, 0, 4, 0, 2 * jt, 0)**2 / (2.0d0 * 2 * jt + 1)
                  end if
              end if

              if (ip == jp .and. it == jt) then
                  C = C + Rtar * betnahv2(it2, jt2) / sqrt(4.0d0 * pi)
              end if

              if (it == jt .and. it2 == jt2) then
                  if (IVIBROTP == 0) then
                      C = C + Rpro * betnahvp(ip, jp) / sqrt(4.0d0 * pi)
                  else
                      C = C + Rpro * betp(ip, jp) * sqrt((2 * 2 * ip + 1) * 5.0d0 * (2 * 2 * jp + 1) / (4.0d0 * pi)) &
                          * CG(2 * ip, 0, 2, 0, 2 * jp, 0)**2 / (2.0d0 * 2 * jp + 1)

                      C = C + Rpro * Beta4P * sqrt((2 * 2 * ip + 1) * 9.0d0 * (2 * 2 * jp + 1) / (4.0d0 * pi)) &
                          * CG(2 * ip, 0, 4, 0, 2 * jp, 0)**2 / (2.0d0 * 2 * jp + 1)
                  end if
              end if

              A(i,j) = C


            End Do
          End Do
        End Do

      End Do
    End Do
  End Do

  Call mdiag(A, Nlevel - Ntrans, Ev, Evec)

  Return

End Subroutine coupled_matrix0

Subroutine coupled_matrix(r, cpot_matrix)
  Use ccfull_initialization_mod
  Implicit None
  Real(8), Intent(In) :: r
  real(8), intent(out) :: cpot_matrix(Nlevelmax, Nlevelmax)
  Real(8) :: C, A, C0
  Real(8) :: Fct, Fct2, Fct3, Fct4, Fct2v, Fctt
  Real(8) :: Ftrans, Fcp, Fcp2, Fcp3, Fcp4, Fcp2v
  Real(8) :: CG, Vn, VnCC
  integer :: i, j, k, ip, it, it2, jp, jt, jt2
  external CG, Fct, Fct2, Fct3, Fct4, Fct2v, Fctt
  external Ftrans, Fcp, Fcp2, Fcp3, Fcp4, Fcp2v, Vn, VnCC
  
  Do i = 1, nlevelmax
    eps(i) = 0.
  End Do

  Do i = 1, Nlevelmax
    Do j = 1, Nlevelmax
     cpot_matrix(i,j) = 0.0d0
    End Do
  End Do

  If((Nlevel - Ntrans) == 1)Return

  i = 0
  Do ip = 0, Npro
    Do it = 0, Ntar
      Do it2 = 0, NphononT2

        if (NphononT2 .Ne. 0)Then
          if(imutual(it, it2) == 0.)Cycle
        End if

        i = i + 1 
        j = 0
        Do jp = 0, Npro
          Do jt = 0, Ntar
            Do jt2 = 0, NphononT2

            if (NphononT2 .Ne. 0)Then
              if(imutual(jt, jt2) == 0.)Cycle
            End if

            j = j + 1 

            If (i > j)Then

              cpot_matrix(i,j) = cpot_matrix(j,i)
              Cycle
            End If

            C = 0.0d0
            ! nuclear coupling
            Do k = 1, Nlevel-Ntrans
              C = C + VnCC(r, Ev(k))*Evec(i, k)*Evec(j, k)

            End Do

            ! ---- Target contribution ----
            if (ip == jp .and. it2 == jt2) then
                if (IVIBROTT == 0) then
                    if (it /= jt) then
                        A = betcahv(it, jt) * Fct(r)
                        if (BetaT /= 0.d0) A = A / BetaT
                        
                    else
                        A = betcahv(it, jt) * Fct2v(r) / sqrt(4.d0 * pi)

                    end if
                    C = C + A
                else
                    C = C + sqrt((2 * 2 * it + 1) * 5.d0 * (2 * 2 * jt + 1) / (4.d0 * pi)) &
                          * CG(2 * it, 0, 2, 0, 2 * jt, 0)**2 / (2.d0 * 2 * jt + 1) * Fct2(r) &
                          * (bett(it, jt) + 2.d0 * sqrt(5.d0 / pi) * bett(it, jt)**2 / 7.d0)

                    C = C + sqrt((2 * 2 * it + 1) * 9.d0 * (2 * 2 * jt + 1) / (4.d0 * pi)) &
                          * CG(2 * it, 0, 4, 0, 2 * jt, 0)**2 / (2.d0 * 2 * jt + 1) * Fct4(r)
                end if
            end if

            ! ---- Projectile contribution ----
            if (it == jt .and. it2 == jt2) then
                if (IVIBROTP == 0) then
                    if (ip /= jp) then
                        A = betcahvp(ip, jp) * Fcp(r)
                        if (BetaP /= 0.d0) A = A / BetaP
                    else
                        A = betcahvp(ip, jp) * Fcp2v(r) / sqrt(4.d0 * pi)
                    end if
                    C = C + A
                else
                    C = C + sqrt((2 * 2 * ip + 1) * 5.d0 * (2 * 2 * jp + 1) / (4.d0 * pi)) &
                          * CG(2 * ip, 0, 2, 0, 2 * jp, 0)**2 / (2.d0 * 2 * jp + 1) * Fcp2(r) &
                          * (betp(ip, jp) + 2.d0 * sqrt(5.d0 / pi) * betp(ip, jp)**2 / 7.d0)

                    C = C + sqrt((2 * 2 * ip + 1) * 9.d0 * (2 * 2 * jp + 1) / (4.d0 * pi)) &
                          * CG(2 * ip, 0, 4, 0, 2 * jp, 0)**2 / (2.d0 * 2 * jp + 1) * Fcp4(r)
                end if
            end if

            ! ---- BetaT2 part (vibrational excitation of target) ----
            if (ip == jp .and. it == jt) then
                if (it2 /= jt2) then
                    A = betcahv2(it2, jt2) * Fctt(r)
                    if (BetaT2 /= 0.d0) A = A / BetaT2
                else
                    A = betcahv2(it2, jt2) * Fct2v(r) / sqrt(4.d0 * pi)
                end if
                C = C + A
            end if

            ! ---- Excitation energy shift ----
            If (it == jt .and. ip == jp .and. it2 == jt2) then
                if (IVIBROTT == 0) then
                    C = C + omeahv(it)
                    eps(i) = eps(i) + omeahv(it)
                else
                    C = C + erott(it)
                    eps(i) = eps(i) + erott(it)
                end if

                if (IVIBROTP == 0) then
                    C = C + omeahvp(ip)
                    eps(i) = eps(i) + omeahvp(ip)
                else
                    C = C + erotp(ip)
                    eps(i) = eps(i) + erotp(ip)
                end if

                C = C + omeahv2(it2)
                eps(i) = eps(i) + omeahv2(it2)
            End if

            cpot_matrix(i, j) = C


            End Do
          End Do
        End Do

      End Do
    End Do
  End Do

  C0 = cpot_matrix(1,1)
  Do i = 1, Nlevel-Ntrans
    cpot_matrix(i,i) = cpot_matrix(i,i) - Vn(r)
    !cpot_matrix(i,i) = cpot_matrix(i,i) - C0
  End Do

  ! transfer coupling
  If (Ntrans == 1) Then
      cpot_matrix(1, Nlevel) = Ftrans(r)
      cpot_matrix(Nlevel, 1) = Ftrans(r)
      cpot_matrix(Nlevel, Nlevel) = cpot_matrix(1, 1) - Qtrans + (Zpro + iq) * (Ztar - iq) &
                                        / r * Hbar / 137.0 - Zpro * Ztar / r * Hbar / 137.0
  End If

  Return

End Subroutine coupled_matrix

Subroutine rkutta00(psi0, phi0, psi1)
  Use ccfull_initialization_mod
  Implicit None
  complex*16, Intent(In)  ::  psi0(Nlevelmax), phi0(Nlevelmax)
  complex*16, intent(out) ::  psi1(Nlevelmax)
  complex*16 :: ai
  complex*16, dimension(Nlevelmax) :: ak1, ak2, ak3, ak4
  complex*16, dimension(Nlevelmax) :: bk1, bk2, bk3, bk4

  Real(8) ::  fac, r, rh, rpp, V
  Integer ::  j1, i0, ic, is
  external V

  ai = (0.d0, 1.d0)
  fac = dr * (2.d0 * ReduceMass / Hbar**2)

  r  = R_min
  rh = R_min + dr/2.d0
  rpp= R_min + dr
  Do j1 = 1, Nlevel
      ak1(j1) = 0.0D0
      ak2(j1) = 0.0D0
      ak3(j1) = 0.0D0
      ak4(j1) = 0.0D0
      bk1(j1) = 0.0D0
      bk2(j1) = 0.0D0
      bk3(j1) = 0.0D0
      bk4(j1) = 0.0D0
  End Do

  Do i0 = 1, Nlevel
      Do ic = 1, Nlevel
        ak1(i0) = ak1(i0) + fac * CPOT(i0, ic, 0) * psi0(ic)
      End Do
      ak1(i0) = ak1(i0) - fac * (E - V(r, L_i)) * psi0(i0)
      bk1(i0) = dr * phi0(i0)
  End Do

  Do i0 = 1, Nlevel
      Do ic = 1, Nlevel
        ak2(i0) = ak2(i0) + fac * CPOTH(i0, ic) * (psi0(ic) + 1.0D0 / 2.0D0 * bk1(ic))
      End Do
      ak2(i0) = ak2(i0) - fac * (E - V(rh, L_i)) * (psi0(i0) + 1.0D0 / 2.0D0 * bk1(i0))
      bk2(i0) = dr * (phi0(i0) + 1.0D0 / 2.0D0 * ak1(i0))
  End Do

  Do i0 = 1, Nlevel
      Do ic = 1, Nlevel
        ak3(i0) = ak3(i0) + fac * CPOTH(i0, ic) * (psi0(ic) + 1.0D0 / 2.0D0 * bk2(ic))
      End Do
      ak3(i0) = ak3(i0) - fac * (E - V(rh, L_i)) * (psi0(i0) + 1.0D0 / 2.0D0 * bk2(i0))
      bk3(i0) = dr * (phi0(i0) + 1.0D0 / 2.0D0 * ak2(i0))
  End Do

  Do i0 = 1, Nlevel
      Do ic = 1, Nlevel
        ak4(i0) = ak4(i0) + fac * CPOT(i0, ic, 1) * (psi0(ic) + bk3(ic))
      End Do
      ak4(i0) = ak4(i0) - fac * (E - V(rpp, L_i)) * (psi0(i0) + bk3(i0))
      bk4(i0) = dr * (phi0(i0) + ak3(i0))
  End Do

  Do is = 1, nlevel
      psi1(is) = psi0(is) + (1.0D0 / 6.0D0) * (bk1(is) + 2.0D0 * bk2(is) + 2.0D0 * bk3(is) + bk4(is))
  End Do

  Return

End Subroutine rkutta00

Subroutine stabilize(xi1, xi, aa, ir1)
  Use ccfull_initialization_mod
  Implicit None
  Integer, Intent(In)  ::ir1
  complex*16, Intent(Inout)  ::  aa(Nlevelmax, Nlevelmax)
  complex*16, Intent(Inout)  ::  xi(Nlevelmax, Nlevelmax)
  complex*16, Intent(Inout)  ::  xi1(Nlevelmax, Nlevelmax)
  complex*16 :: ai
  complex*16, dimension(Nlevelmax, Nlevelmax)  :: psi
  complex*16, dimension(Nlevelmax, Nlevelmax)  :: cc, cin
  complex*16, dimension(Nlevelmax, Nlevelmax)  :: aa0
  complex*16, dimension(Nlevelmax, Nlevelmax)  :: xid, xid1
  Real(8) ::  r, V, fac
  Integer ::  i, j, k, i0, ic, ich
  external  V

  Do i = 1, nlevel
      Do j = 1, nlevel
        aa0(i, j)  = aa(i, j)
        xid(i, j)  = xi(i, j)
        xid1(i, j) = xi1(i, j)
      End Do
  End Do

  ai = (0.d0, 1.d0)
  fac = dr**2 * (2.d0 * ReduceMass / Hbar**2)
  r = R_min + ir1*dr

  Do i0 = 1, Nlevel
    Do ic = 1, Nlevel
        cc(i0, ic) = -fac / 12.0D0 * CPOT(i0, ic, ir1)
        If (i0 == ic) Then
          cc(i0, ic) = cc(i0, ic) - fac / 12.0D0 * (V(r, L_i) - E) + 1.0D0
        End If
    End Do
  End Do

  Call matinv(Nlevel,cc,cin)

  Do ich = 1, nlevel
      Do i0 = 1, nlevel
        psi(i0, ich) = 0.0D0
        Do ic = 1, nlevel
            psi(i0, ich) = psi(i0, ich) + cin(i0, ic) * xi(ic, ich)
        End Do
      End Do
  End Do

  Do i = 1, nlevel
      Do j = 1, nlevel
        aa(i, j) = 0.0D0
        Do k = 1, nlevel
            aa(i, j) = aa(i, j) + psi(i, k) * aa0(k, j)
        End Do
      End Do
  End Do

  call matinv(Nlevel,psi,cin)

  Do i = 1, Nlevel
      Do j = 1, Nlevel
        xi(i, j)  = 0.0D0
        xi1(i, j) = 0.0D0
        Do k = 1, Nlevel
            xi(i, j)  = xi(i, j) + xid(i, k)  * cin(k, j)
            xi1(i, j) = xi1(i, j) + xid1(i, k) * cin(k, j)
        End Do
      End Do
  End Do

  Return
End Subroutine stabilize

Subroutine Numerov(P)
  Use ccfull_initialization_mod
  Implicit None
  real(8), intent(out) :: P
  complex*16, dimension(Nlevelmax) :: psi, psi0, psi1
  complex*16, dimension(Nlevelmax, Nlevelmax) :: xi, xi0, xi1
  complex*16, dimension(Nlevelmax) :: phi0
  complex*16, dimension(Nlevelmax, Nlevelmax) :: bb, bin
  complex*16, dimension(Nlevelmax, Nlevelmax) :: bb2
  complex*16, dimension(Nlevelmax, Nlevelmax) :: cc, cin
  complex*16, dimension(Nlevelmax, Nlevelmax) :: dd0, dd1
  complex*16, dimension(Nlevelmax) :: dd
  real(8), dimension(0:200) :: fcw, gcw, fpcw, gpcw
  real(8), dimension(0:200) :: sigmad
  integer, dimension(0:200) :: iexp
  real(8), dimension(Nlevelmax) :: ech, ech2
  complex*16 :: k, kk, k2
  complex*16 :: ai
  complex*16 :: cwup0, cwdown0, cwup1, cwdown1
  complex*16 :: dummy, dummy2
  complex*16, dimension(Nlevelmax, Nlevelmax) :: xi1d, xi0d
  complex*16, dimension(Nlevelmax, Nlevelmax) :: aa, bb0, bb20
  Integer ::  i, j, lc, io, io2, i0, ir, ii, ik, ic, ich, ibarrier
  Integer ::  j1, j2
  Real(8) ::  V, fac, r, r1, r2
  Real(8) ::  ak, ec, eta, rho
  external V

  ai = (0.d0, 1.d0)
  fac = dr**2 * (2.d0 * ReduceMass / Hbar**2)
  ibarrier = (R_barrier - R_min)/dr
  Do lc = 0, 200
    fcw(lc)   = 0.0d0
    gcw(lc)   = 0.0d0
    fpcw(lc)  = 0.0d0
    gpcw(lc)  = 0.0d0
    sigmad(lc) = 0.0d0
    iexp(lc)   = 0
  End Do

  Do i = 1, Nlevelmax
      Do j = 1, Nlevelmax
          bb(i,j)  = 0.0d0
          bin(i,j) = 0.0d0
          cc(i,j)  = 0.0d0
          cin(i,j) = 0.0d0
          aa(i,j)  = 0.0d0
      End Do
      dd(i) = 0.0d0
  End Do

  Do i = 1, Nlevel
      aa(i,i) = 1.0d0
  End Do


  Do io = 1, Nlevel

    Do j1 = 1, Nlevel
        psi(j1)  = 0.0d0
        psi0(j1) = 0.0d0
        psi1(j1) = 0.0d0
        phi0(j1) = 0.0d0
    End Do

    If (io == 1) Then
        Do io2 = 1, Nlevel
            ech(io2) = E - V(R_min, L_i) - CPOT(io2, io2, 0)
            ech2(io2) = E - CPOT(io2, io2, R_iterat)

        End Do
    End If

    If (ech(io) > 0.d0) Then
        k = sqrt(2.d0 * ReduceMass / Hbar**2 * ech(io))
        psi0(io) = exp(-ai * k * R_min)
        phi0(io) = -ai * k * psi0(io)

    Else
        k = sqrt(2.d0 * ReduceMass / hbar**2 * abs(ech(io)))
        psi0(io) = exp(k * R_min)
        phi0(io) = k * psi0(io)
    End If


    call rkutta00(psi0,phi0,psi1)

    do i0 = 1, Nlevel
        xi0(i0, io) = (1.d0 - fac/12.d0 * (V(R_min, L_i) - e)) * psi0(i0)
        xi1(i0, io) = (1.d0 - fac/12.d0 * (V(R_min+dr, L_i) - e)) * psi1(i0)
        do ic = 1, Nlevel
            xi0(i0, io) = xi0(i0, io) - fac/12.d0 * CPOT(i0, ic, 0) * psi0(ic)
            xi1(i0, io) = xi1(i0, io) - fac/12.d0 * CPOT(i0, ic, 1) * psi1(ic)
        end do
    end do
  End Do

  !----------------------------------------------  iterations start
  Do ir = 2, R_iterat + 1
    r = R_min + dr * ir
    !r0 = R_min + dr * (ir - 2)
    r1 = R_min + dr * (ir - 1)

    Do i0 = 1, Nlevel
        Do ic = 1, Nlevel
            dd0(i0, ic) = fac / sqrt(12.d0) * CPOT(i0, ic, ir - 1)
            If (i0 == ic) Then
                dd0(i0, ic) = dd0(i0, ic) + fac / sqrt(12.d0) * (V(r1, L_i) - E) + sqrt(3.d0)
            End If
        End Do
    End Do

    Do i0 = 1, Nlevel
        Do ic = 1, Nlevel
            dd1(i0, ic) = 0.d0
            If (i0 == ic) Then
                dd1(i0, ic) = dd1(i0, ic) - 1.d0
            End If
            Do ik = 1, Nlevel
                dd1(i0, ic) = dd1(i0, ic) + dd0(i0, ik) * dd0(ik, ic)
            End Do
        End Do
    End Do

    Do ich = 1, Nlevel
        Do i0 = 1, Nlevel
            xi(i0, ich) = -xi0(i0, ich)
            Do ic = 1, nlevel
                xi(i0, ich) = xi(i0, ich) + dd1(i0, ic) * xi1(ic, ich)
            End Do
        End Do
    End Do

    If(ir == R_iterat+1)Exit
    If(ir == ibarrier) Call stabilize(xi1,xi,aa,ir)
    

    Do ich = 1, nlevel
        Do i0 = 1, nlevel
            xi0(i0, ich) = xi1(i0, ich)
            xi1(i0, ich) = xi(i0, ich)
        End Do
    End Do


  End Do


  !--------------------------------------------------------------
  !  matching to the coulomb wave function at rmax

  Do io = 1, Nlevel

    Do i0 = 1, Nlevel
        Do ic = 1, Nlevel
            cc(i0, ic) = -fac / 12.d0 * CPOT(i0, ic, R_iterat - 1)
            If (i0 == ic) Then
                cc(i0, ic) = cc(i0, ic) - fac / 12.d0 * (V(R_max - dr, L_i) - E) + 1.d0
            End If
        End Do
    End Do


    Call matinv(Nlevel,cc,cin)
    Do i0 = 1, Nlevel
        psi0(i0) = 0.d0
        Do ic = 1, Nlevel
            psi0(i0) = psi0(i0) + cin(i0, ic) * xi0(ic, io)
        End Do
    End Do


    Do i0 = 1, Nlevel
        Do ic = 1, Nlevel
            cc(i0, ic) = -fac / 12.d0 * CPOT(i0, ic, R_iterat + 1)
            If (i0 == ic) Then
                cc(i0, ic) = cc(i0, ic) - fac / 12.d0 * (V(R_max + dr, L_i) - E) + 1.d0
            End If
        End Do
    End Do

    Call matinv(Nlevel,cc,cin)

    Do i0 = 1, Nlevel
        psi(i0) = 0.d0
        Do ic = 1, Nlevel
            psi(i0) = psi(i0) + cin(i0, ic) * xi(ic, io)
        End Do
    End Do


    Do ii = 1, Nlevel
        ! coulomb wave function
        ec = E - CPOT(ii, ii, R_iterat - 1)

        If (ec < 0.d0) Then
            r1 = R_max - dr
            r2 = R_max + dr
            ak = sqrt(2.d0 * ReduceMass * abs(ec) / Hbar**2)
            
            bb0(ii, io) = (exp(-ak * r2) * psi0(ii) - exp(-ak * r1) * psi(ii)) &
                          / (exp(ak * (r1 - r2)) - exp(-ak * (r1 - r2)))
            bb20(ii, io) = -(exp(ak * r2) * psi0(ii) - exp(ak * r1) * psi(ii)) &
                          / (exp(ak * (r1 - r2)) - exp(-ak * (r1 - r2)))
        Else
            rho = sqrt(2.d0 * ReduceMass * ec) / hbar * (R_max - dr)
            eta = (Zpro * Ztar / 137.d0) * Sqrt(ReduceMass / (2.d0 * ec))


            !Call dfcoul(eta, rho, fcw, fpcw, gcw, gpcw, sigmad, L_i, iexp)
            Call myCOULFG(L_i*1.d0, eta, rho, fcw(L_i), fpcw(L_i), gcw(L_i), gpcw(L_i), iexp(L_i))
  
 
            cwup0 = (gcw(L_i) + ai * fcw(L_i))
            cwdown0 = gcw(L_i) - ai * fcw(L_i)

            ec = E - CPOT(ii, ii, R_iterat + 1)
            rho = Sqrt(2.d0 * ReduceMass * ec) / Hbar * (R_max + dr)
            eta = (Zpro * Ztar / 137.d0) * Sqrt(ReduceMass / (2.d0 * ec))
            !Call dfcoul(eta, rho, fcw, fpcw, gcw, gpcw, sigmad, L_i, iexp)
            Call myCOULFG(L_i*1.d0, eta, rho, fcw(L_i), fpcw(L_i), gcw(L_i), gpcw(L_i), iexp(L_i))
            cwup1 = (gcw(L_i) + ai * fcw(L_i))
            cwdown1 = gcw(L_i) - ai * fcw(L_i)

            bb0(ii, io) = (cwup0 * psi(ii) - cwup1 * psi0(ii)) &
                          / (cwup0 * cwdown1 - cwup1 * cwdown0)
            bb20(ii, io) = (cwdown1 * psi0(ii) - cwdown0 * psi(ii)) &
                          / (cwup0 * cwdown1 - cwup1 * cwdown0)
        End If
    End Do

  End Do
  !===============================================================
  !                                        penetration probability

  Do i = 1, Nlevel
      Do j = 1, Nlevel
        bb(i,j) = 0.d0
        bb2(i,j) = 0.d0
        Do j2 = 1, nlevel
            bb(i,j) = bb(i,j) + bb0(i,j2) * aa(j2,j)
            bb2(i,j) = bb2(i,j) + bb20(i,j2) * aa(j2,j)
        End Do
      End Do
  End Do

  Call matinv(Nlevel,bb,bin)

  P = 0.d0

  Do io = 1, Nlevel
      If (ech(io) .lt. 0.d0) Cycle
      k = sqrt((2.d0 * ReduceMass / Hbar**2 * ech(io)))
      kk = sqrt(2.d0 * ReduceMass / Hbar**2 * e)
      p = p + (abs(bin(io,1)))**2 * k / kk
  End Do

  Return

End Subroutine Numerov

Subroutine Numerov_Ec(P)
  Use ccfull_initialization_mod
  Implicit None
  real(8), intent(out) :: P
  complex*16, dimension(Nlevelmax) :: psi, psi0, psi1
  complex*16, dimension(Nlevelmax, Nlevelmax) :: xi, xi0, xi1
  complex*16, dimension(Nlevelmax) :: phi0
  complex*16, dimension(Nlevelmax, Nlevelmax) :: bb, bin
  complex*16, dimension(Nlevelmax, Nlevelmax) :: bb2
  complex*16, dimension(Nlevelmax, Nlevelmax) :: cc, cin
  complex*16, dimension(Nlevelmax, Nlevelmax) :: dd0, dd1
  complex*16, dimension(Nlevelmax) :: dd
  real(8), dimension(0:200) :: fcw, gcw, fpcw, gpcw
  real(8), dimension(0:200) :: sigmad
  integer, dimension(0:200) :: iexp
  real(8), dimension(Nlevelmax) :: ech, ech2
  complex*16 :: ai
  complex*16 :: cwup0, cwdown0, cwup1, cwdown1
  complex*16, dimension(Nlevelmax, Nlevelmax) :: xi1d, xi0d
  complex*16, dimension(Nlevelmax, Nlevelmax) :: aa, bb0, bb20
  complex*16, dimension(Nlevelmax, Nlevelmax, 0:R_iterat+1) :: wf_mn
  Integer ::  i, j, lc, io, io2, i0, ir, ii, ik, ic, ich, ibarrier
  Integer ::  j1, j2, istab
  Real(8) ::  V, fac, r, r1, r2, k, kk
  Real(8) ::  ak, ec, eta, rho, dummy
  external V

  istab=20
  wf_mn = (0.0d0,0.0d0)
  wf    = (0.0d0,0.0d0)

  ai = (0.d0, 1.d0)
  fac = dr**2 * (2.d0 * ReduceMass / Hbar**2)
  ibarrier = (R_barrier - R_min)/dr
  Do lc = 0, 200
    fcw(lc)   = 0.0d0
    gcw(lc)   = 0.0d0
    fpcw(lc)  = 0.0d0
    gpcw(lc)  = 0.0d0
    sigmad(lc) = 0.0d0
    iexp(lc)   = 0
  End Do

  Do i = 1, Nlevelmax
      Do j = 1, Nlevelmax
          bb(i,j)  = 0.0d0
          bin(i,j) = 0.0d0
          cc(i,j)  = 0.0d0
          cin(i,j) = 0.0d0
          aa(i,j)  = 0.0d0
      End Do
      dd(i) = 0.0d0
  End Do

  Do i = 1, Nlevel
      aa(i,i) = 1.0d0
  End Do


  Do io = 1, Nlevel

    Do j1 = 1, Nlevel
        psi(j1)  = 0.0d0
        psi0(j1) = 0.0d0
        psi1(j1) = 0.0d0
        phi0(j1) = 0.0d0
    End Do

    If (io == 1) Then
        Do io2 = 1, Nlevel
            ech(io2) = E - V(R_min, L_i) - CPOT(io2, io2, 0)
            ech2(io2) = E - CPOT(io2, io2, R_iterat)
        End Do
    End If

    If (ech(io) > 0.d0) Then
        k = sqrt(2.d0 * ReduceMass / Hbar**2 * ech(io))
        psi0(io) = exp(-ai * k * R_min)
        phi0(io) = -ai * k * psi0(io)
    Else
        k = sqrt(2.d0 * ReduceMass / hbar**2 * abs(ech(io)))
        psi0(io) = exp(k * R_min)
        phi0(io) = k * psi0(io)
    End If

    call rkutta00(psi0,phi0,psi1)

    do i0 = 1, Nlevel
        xi0(i0, io) = (1.d0 - fac/12.d0 * (V(R_min, L_i) - e)) * psi0(i0)
        xi1(i0, io) = (1.d0 - fac/12.d0 * (V(R_min+dr, L_i) - e)) * psi1(i0)
        do ic = 1, Nlevel
            xi0(i0, io) = xi0(i0, io) - fac/12.d0 * CPOT(i0, ic, 0) * psi0(ic)
            xi1(i0, io) = xi1(i0, io) - fac/12.d0 * CPOT(i0, ic, 1) * psi1(ic)
        end do
    end do
  End Do


  ! output the wave function
  Do io = 1, Nlevel
      Do i0 = 1, Nlevel
          Do ic = 1, Nlevel
              cc(i0, ic) = -fac / 12.d0 * CPOT(i0, ic, 1)
              If (i0 == ic) Then
                  cc(i0, ic) = cc(i0, ic) - fac / 12.d0 * (V(R_min+dr, L_i) - E) + 1.d0
              End If
          End Do
      End Do
      Call matinv(Nlevel,cc,cin)
      Do i0 = 1, nlevel
        Do ic = 1, Nlevel
          wf_mn(i0, io, 1) = wf_mn(i0, io, 1) + cin(i0, ic) * xi1(ic, io)
        End Do
      End Do
  End Do

  Do io = 1, Nlevel
      Do i0 = 1, Nlevel
          Do ic = 1, Nlevel
              cc(i0, ic) = -fac / 12.d0 * CPOT(i0, ic, 0)
              If (i0 == ic) Then
                  cc(i0, ic) = cc(i0, ic) - fac / 12.d0 * (V(R_min, L_i) - E) + 1.d0
              End If
          End Do
      End Do
      Call matinv(Nlevel,cc,cin)
      Do i0 = 1, nlevel
        Do ic = 1, Nlevel
          wf_mn(i0, io, 0) = wf_mn(i0, io, 0) + cin(i0, ic) * xi0(ic, io)
        End Do
      End Do
  End Do



 ! Do ich = 1, nlevel
 !   Do i0 = 1, Nlevel
 !     wf_mn(ich, i0, 0) = xi0(ich, i0)
 !     wf_mn(ich, i0, 1) = xi1(ich, i0)
 !   End Do
 ! End Do

  !wf_mn(1, 1, 0) = psi0(1)
  !wf_mn(1, 1, 1) = psi1(1)
  
  

  !----------------------------------------------  iterations start
  Do ir = 2, R_iterat + 1
    r = R_min + dr * ir
    !r0 = R_min + dr * (ir - 2)
    r1 = R_min + dr * (ir - 1)

    Do i0 = 1, Nlevel
        Do ic = 1, Nlevel
            dd0(i0, ic) = fac / sqrt(12.d0) * CPOT(i0, ic, ir - 1)
            If (i0 == ic) Then
                dd0(i0, ic) = dd0(i0, ic) + fac / sqrt(12.d0) * (V(r1, L_i) - E) + sqrt(3.d0)
            End If
        End Do
    End Do
 
    Do i0 = 1, Nlevel
        Do ic = 1, Nlevel
            dd1(i0, ic) = 0.d0
            If (i0 == ic) Then
                dd1(i0, ic) = dd1(i0, ic) - 1.d0
            End If
            Do ik = 1, Nlevel
                dd1(i0, ic) = dd1(i0, ic) + dd0(i0, ik) * dd0(ik, ic)
            End Do
        End Do
    End Do
 
    Do ich = 1, Nlevel
        Do i0 = 1, Nlevel
            xi(i0, ich) = -xi0(i0, ich)
            Do ic = 1, nlevel
                xi(i0, ich) = xi(i0, ich) + dd1(i0, ic) * xi1(ic, ich)
            End Do
        End Do
    End Do


    ! output the wave function
    Do io = 1, Nlevel
        Do i0 = 1, Nlevel
            Do ic = 1, Nlevel
                cc(i0, ic) = -fac / 12.d0 * CPOT(i0, ic, ir)
                If (i0 == ic) Then
                    cc(i0, ic) = cc(i0, ic) - fac / 12.d0 * (V(r, L_i) - E) + 1.d0
                End If
            End Do
        End Do
        Call matinv(Nlevel,cc,cin)
        Do i0 = 1, nlevel
          wf_mn(i0, io, ir) = (0.d0, 0.d0)
          Do ic = 1, Nlevel
            wf_mn(i0, io, ir) = wf_mn(i0, io, ir) + cin(i0, ic) * xi(ic, io)
          End Do
        End Do
    End Do

    !Do ich = 1, nlevel
    !  Do i0 = 1, Nlevel
    !    wf_mn(ich, i0, ir) = xi(ich, i0)
    !  End Do
    !End Do
    ! output the wave function

    If(ir == R_iterat+1)Exit
    !If(ir == ibarrier) Call stabilize(xi1,xi,aa,ir)

    Do ich = 1, nlevel
        Do i0 = 1, nlevel
            xi0(i0, ich) = xi1(i0, ich)
            xi1(i0, ich) = xi(i0, ich)
        End Do
    End Do



  End Do


  !--------------------------------------------------------------
  !  matching to the coulomb wave function at rmax

  Do io = 1, Nlevel

    Do i0 = 1, Nlevel
        Do ic = 1, Nlevel
            cc(i0, ic) = -fac / 12.d0 * CPOT(i0, ic, R_iterat - 1)
            If (i0 == ic) Then
                cc(i0, ic) = cc(i0, ic) - fac / 12.d0 * (V(R_max - dr, L_i) - E) + 1.d0
            End If
        End Do
    End Do


    Call matinv(Nlevel,cc,cin)
    Do i0 = 1, Nlevel
        psi0(i0) = 0.d0
        Do ic = 1, Nlevel
            psi0(i0) = psi0(i0) + cin(i0, ic) * xi0(ic, io)
        End Do
    End Do


    Do i0 = 1, Nlevel
        Do ic = 1, Nlevel
            cc(i0, ic) = -fac / 12.d0 * CPOT(i0, ic, R_iterat + 1)
            If (i0 == ic) Then
                cc(i0, ic) = cc(i0, ic) - fac / 12.d0 * (V(R_max + dr, L_i) - E) + 1.d0
            End If
        End Do
    End Do

    Call matinv(Nlevel,cc,cin)

    Do i0 = 1, Nlevel
        psi(i0) = 0.d0
        Do ic = 1, Nlevel
            psi(i0) = psi(i0) + cin(i0, ic) * xi(ic, io)
        End Do
    End Do


    Do ii = 1, Nlevel
        ! coulomb wave function
        ec = E - CPOT(ii, ii, R_iterat - 1) 

        If (ec < 0.d0) Then
            r1 = R_max - dr
            r2 = R_max + dr
            ak = sqrt(2.d0 * ReduceMass * abs(ec) / Hbar**2)
            
            bb0(ii, io) = (exp(-ak * r2) * psi0(ii) - exp(-ak * r1) * psi(ii)) &
                          / (exp(ak * (r1 - r2)) - exp(-ak * (r1 - r2)))
            bb20(ii, io) = -(exp(ak * r2) * psi0(ii) - exp(ak * r1) * psi(ii)) &
                          / (exp(ak * (r1 - r2)) - exp(-ak * (r1 - r2)))
        Else
            rho = sqrt(2.d0 * ReduceMass * ec) / hbar * (R_max - dr)
            eta = (Zpro * Ztar / 137.d0) * Sqrt(ReduceMass / (2.d0 * ec))
            

            !Call dfcoul(eta, rho, fcw, fpcw, gcw, gpcw, sigmad, L_i, iexp)
            Call myCOULFG(L_i*1.d0, eta, rho, fcw(L_i), fpcw(L_i), gcw(L_i), gpcw(L_i), iexp(L_i))
  
 
            cwup0 = (gcw(L_i) + ai * fcw(L_i))
            cwdown0 = gcw(L_i) - ai * fcw(L_i)

            ec = E - CPOT(ii, ii, R_iterat + 1) 
            rho = Sqrt(2.d0 * ReduceMass * ec) / Hbar * (R_max + dr)
            
            eta = (Zpro * Ztar / 137.d0) * Sqrt(ReduceMass / (2.d0 * ec))
            !Call dfcoul(eta, rho, fcw, fpcw, gcw, gpcw, sigmad, L_i, iexp)
            Call myCOULFG(L_i*1.d0, eta, rho, fcw(L_i), fpcw(L_i), gcw(L_i), gpcw(L_i), iexp(L_i))
            cwup1 = (gcw(L_i) + ai * fcw(L_i))
            cwdown1 = gcw(L_i) - ai * fcw(L_i)

            bb0(ii, io) = (cwup0 * psi(ii) - cwup1 * psi0(ii)) &
                          / (cwup0 * cwdown1 - cwup1 * cwdown0)
            bb20(ii, io) = (cwdown1 * psi0(ii) - cwdown0 * psi(ii)) &
                          / (cwup0 * cwdown1 - cwup1 * cwdown0)
        End If
    End Do
    
  End Do

 
  !===============================================================
  !                                        penetration probability

  Do i = 1, Nlevel
      Do j = 1, Nlevel
        bb(i,j) = 0.d0
        bb2(i,j) = 0.d0
        Do j2 = 1, nlevel
            bb(i,j) = bb(i,j) + bb0(i,j2) * aa(j2,j)        !C_nm
            bb2(i,j) = bb2(i,j) + bb20(i,j2) * aa(j2,j)     !D_nm
        End Do
      End Do
  End Do

  Call matinv(Nlevel,bb,bin)

  P = 0.d0
  ! Penetration probability

  Do io = 1, Nlevel
      If (ech(io) .lt. 0.d0) Cycle
      k = sqrt((2.d0 * ReduceMass / Hbar**2 * ech(io)))
      kk = sqrt(2.d0 * ReduceMass / Hbar**2 * e)
      p = p + (abs(bin(io,1)))**2 * k / kk
      !Write(*,*)'Numerove:', (abs(bin(io,1)))
  End Do

  ! Wavefunction 
  wf = (0.d0, 0.d0)
  Do ir = 0, R_iterat + 1
      Do i0 = 1, Nlevel
        Do ic = 1, Nlevel
          wf(ir*Nlevel+i0) =  wf(ir*Nlevel+i0) + bin(ic,1) * wf_mn(i0, ic, ir)!wf_mn(i0, ic, ir) !
        End Do
      End Do

  End Do

  Return

End Subroutine Numerov_Ec

Subroutine Prameterize_H0(para1, para2, H2, dbpsi,dbpsi0)
  Use ccfull_initialization_mod
  implicit none
  !H2, dbpsi, dbpsi0
  complex*16, dimension(Nlevelmax, Nlevelmax) :: dd0, dd1
  complex*16, dimension(Nlevelmax, Nlevelmax) :: B, B2, C, cc, tt, ee, Hc
  complex*16 H1(Nx_cc,Nx_cc)
  complex*16 H2(Nx1, Nx2)
  complex*16 dbpsi(Nx2, Nx1)
  complex*16 dbpsi0(Nx2)
  complex*16 ai
  real(8), dimension(Nlevelmax, Nlevelmax, -1:R_iterat + 1) :: A

  Integer ::  i, j, ir, ie, il, ik, ic, i0, idx, jdx
  Real(8) ::  r,rh, V, eta, rho, k, ec, fac,temp
  real(8), intent(in) :: para1, para2
  Real(8), dimension(0:200) :: fcw, gcw, fpcw, gpcw
  Real(8), dimension(0:200) :: sigmad
  Integer, dimension(0:200) :: iexp
  external V
  H1 = (0.d0, 0.d0)
  H2 = (0.d0, 0.d0)
  ai = (0.d0, 1.d0)
  fac = dr**2*(2.d0*ReduceMass/Hbar**2) 

  !BetaT = para1
  Beta2T = para1
  Beta4T = para2
  ! contruct the cpot
  Call coupled_matrix0

  !rh = R_min - dr 
  !Call coupled_matrix(rh, CPOTH)
  !Do  i = 1, Nlevel
  !  Do  j = 1, Nlevel
  !    A(i, j, -1) = 2*ReduceMass/Hbar**2 * CPOTH(i, j)
  !    if(i == j)Then
  !      A(i, j, -1) = A(i, j, -1) + 2*ReduceMass/Hbar**2 * (V(rh, L_i) - E)
  !    End If
  !  End Do
  !End Do

  rh = R_min + dr / 2.0
  Call coupled_matrix(rh, CPOTH)


  Do ir = 0, R_iterat + 1
    r = R_min + dr * ir
    Call coupled_matrix(r, CPOT0)
    Do  i = 1, Nlevel
      Do  j = 1, Nlevel
        CPOT(i, j, ir) = CPOT0(i, j)
        A(i, j, ir) = 2*ReduceMass/Hbar**2 * CPOT0(i, j)
        if(i == j)Then
          A(i, j, ir) = A(i, j, ir) + 2*ReduceMass/Hbar**2 * (V(r, L_i) - E)
        End If
      End Do
    End Do
  End Do
  
  !Do ir = 0, R_iterat + 1
  !  r = R_min + dr * ir
  !  
  !  Do  i = 1, Nlevel
  !    Do  j = 1, Nlevel
  !      B(i,j) = 0.0D0
  !      B(i,j) = dr**2 /SQRT(12.) * A(i, j, ir)
  !      if(i == j)B(i,j) = B(i,j) + sqrt(3.d0)
  !    End Do
  !  End Do
  !
  !  Do i = 1, Nlevel
  !    Do j = 1, Nlevel
  !      B2(i,j) = 0.0D0
  !      Do ik = 1, Nlevel
  !        B2(i,j) = B2(i,j) + B(i,ik) * B(ik,j)
  !      End Do
  !      if(i == j)B2(i,j) = B2(i,j) - 1
  !    End Do
  !  End Do
  !
  !  Do i = 1, Nlevel
  !    Do j = 1, Nlevel
  !        cc(i, j) = 0.0D0
  !        cc(i, j) = - dr**2 / 12.d0 * A(i, j, ir)
  !        If (i == j)cc(i, j) = cc(i, j) + 1.d0
  !    End Do
  !  End Do
  !
  !  Do i = 1, Nlevel
  !    Do j = 1, Nlevel
  !      C(i,j) = 0.0D0
  !      Do ik = 1, Nlevel
  !        C(i,j) = C(i,j) + B2(i,ik) * cc(ik,j)
  !      End Do
  !    End Do
  !  End Do
  !
  !  If(ir .lt. R_iterat)Then
  !    Do  i = 1, Nlevel
  !      Do  j = 1, Nlevel
  !        idx = Nlevel * ir + i
  !        jdx = Nlevel * ir + j
  !        H1(idx, jdx) = H1(idx, jdx) - C(i, j)
  !      End Do
  !    End Do
  !  End If
  !
  !  If(ir == R_iterat)Then
  !    Do i = 1, Nlevel
  !      Do j = 1, Nlevel
  !        Hc(i,j) = -C(i,j)
  !      End Do
  !    End Do
  !  End If
  !
  !End Do
  !
  !
  !Do ir = 0, (R_iterat-2)
  !  Do i = 1, Nlevel
  !    Do j = 1, Nlevel
  !        H1(ir * nlevel + i, ir * nlevel + nlevel + j) =  - dr**2 /12. * A(i, j, ir+1)
  !        H1(ir * nlevel + nlevel + i, ir * nlevel + j) =  - dr**2 /12. * A(i, j, ir)
  !        if(i == j)Then
  !          H1(ir * nlevel + i, ir * nlevel + nlevel + j) = H1(ir * nlevel + i, ir * nlevel + nlevel + j) + 1
  !          H1(ir * nlevel + nlevel + i, ir * nlevel + j) = H1(ir * nlevel + nlevel + i, ir * nlevel + j) + 1
  !        End If
  !    End Do
  !  End Do
  !End Do
  !
  !
  !Do i = 1, Nlevel
  !  Do j = 1, Nlevel
  !  Echannel(i) = E - V(R_min, L_i) - CPOT(i,i,0)
  !  ec = E - V(R_min, L_i) - CPOT(i,j,0)
  !
  !  If(ec > 0)Then
  !    !k = acos(1. - ec/2. / t)/dr
  !    k = sqrt(2. * ReduceMass / Hbar**2 * ec)
  !    if(i == j)Then
  !      H1(i, j) = H1(i, j) + (1 - dr**2 /12. * A(i, j, -1)) * exp(ai * k * dr)
  !    else
  !      H1(i, j) = H1(i, j) + (  - dr**2 /12. * A(i, j, -1)) * exp(ai * k * dr)
  !    End If
  !
  !  Else
  !    k = acos(1. - abs(ec)/2. / t)/dr
  !    H1(i, j) = H1(i, j) + ((- dr**2 /12. * A(i, j, -1)) * exp(k * dr))
  !    if(i == j)H1(i, j) = H1(i, j) + 1 * exp(k * dr)
  !
  !  End If 
  !
  !  End do
  !End do
  !
  !Do i = 1, Nlevel*R_iterat
  !    Do j = 1, Nlevel*R_iterat
  !      H2(i,j) = H1(i,j)
  !    End Do
  !End Do
  !
  !Do i = 1, Nlevel
  !  Do j = 1, Nlevel
  !
  !    H2(Nlevel*R_iterat + i, Nlevel*R_iterat + j) = Hc(i, j)
  !    H2(Nlevel*R_iterat - Nlevel + i, Nlevel*R_iterat + j) = (- dr**2 /12. * A(i, j, R_iterat))
  !    H2(Nlevel*R_iterat + i, Nlevel*R_iterat - Nlevel + j) = (- dr**2 /12. * A(i, j, R_iterat-1))
  !    H2(Nlevel*R_iterat + i, Nlevel*R_iterat + Nlevel + j) = (- dr**2 /12. * A(i, j, R_iterat+1))
  !
  !    If(i == j)Then
  !    H2(Nlevel*R_iterat - Nlevel + i, Nlevel*R_iterat + j) = H2(Nlevel*R_iterat - Nlevel + i, Nlevel*R_iterat + j) + 1
  !    H2(Nlevel*R_iterat + i, Nlevel*R_iterat - Nlevel + j) = H2(Nlevel*R_iterat + i, Nlevel*R_iterat - Nlevel + j) + 1
  !    H2(Nlevel*R_iterat + i, Nlevel*R_iterat + Nlevel + j) = H2(Nlevel*R_iterat + i, Nlevel*R_iterat + Nlevel + j) + 1
  !    End If
  ! 
  !  End Do
  !End Do
  !
  !Do i = 1, Nlevel*R_iterat
  !    dbpsi(i, i) = (1.0d0, 0.0d0)
  !End Do
  !
  !
  !Do i = 1, Nlevel
  !  ec = E - CPOT(i,i,R_iterat)
  !  !rho = acos(1. - ec/2. / t)/dr * (R_max)
  !  rho = sqrt(2.d0 * ReduceMass * ec) / Hbar * (R_max)
  !  eta = (Zpro * Ztar / 137.d0) * Sqrt(ReduceMass / (2.d0 * ec))
  !  Call myCOULFG(L_i*1.d0, eta, rho, fcw(L_i), fpcw(L_i), gcw(L_i), gpcw(L_i), iexp(L_i))
  !  dbpsi(Nlevel*R_iterat + i, Nlevel*R_iterat + i) = gcw(L_i) + ai*fcw(L_i)
  !
  !  ec = E - CPOT(i,i,R_iterat+1)
  !  !rho = acos(1. - ec/2. / t)/dr  * (R_max + dr)
  !  rho = sqrt(2.d0 * ReduceMass * ec) / Hbar * (R_max + dr)
  !  eta = (Zpro * Ztar / 137.d0) * Sqrt(ReduceMass / (2.d0 * ec))
  !  Call myCOULFG(L_i*1.d0, eta, rho, fcw(L_i), fpcw(L_i), gcw(L_i), gpcw(L_i), iexp(L_i))
  !  dbpsi(Nlevel*R_iterat + i + Nlevel, Nlevel*R_iterat + i) = gcw(L_i) + ai*fcw(L_i)
  !End Do
  !
  !do i = 1, Nx2
  !    dbpsi0(i) = dconjg(dbpsi(i, Nx1 - Nlevel + 1))
  !end do

End Subroutine Prameterize_H0

Subroutine Prameterize_H(para1, para2, H2, dbpsi,dbpsi0)
  Use ccfull_initialization_mod
  implicit none
  !H2, dbpsi, dbpsi0
  complex*16, dimension(Nlevelmax, Nlevelmax) :: dd0, dd1
  complex*16, dimension(Nlevelmax, Nlevelmax) :: B, B2, C, cc, tt, ee, Hc
  complex*16 H1(Nx_cc,Nx_cc)
  complex*16 H2(Nx1, Nx2)
  complex*16 dbpsi(Nx2, Nx1)
  complex*16 dbpsi0(Nx2)
  complex*16 ai
  real(8), dimension(Nlevelmax, Nlevelmax, -1:R_iterat + 1) :: A

  Integer ::  i, j, ir, ie, il, ik, ic, i0, idx, jdx
  Real(8) ::  r,rh, V, eta, rho, k, ec, fac,temp
  real(8), intent(in) :: para1, para2
  Real(8), dimension(0:200) :: fcw, gcw, fpcw, gpcw
  Real(8), dimension(0:200) :: sigmad
  Integer, dimension(0:200) :: iexp
  external V
  H1 = (0.d0, 0.d0)
  H2 = (0.d0, 0.d0)
  ai = (0.d0, 1.d0)
  fac = dr**2*(2.d0*ReduceMass/Hbar**2) 

  !BetaT = para1
  Beta2T = para1
  Beta4T = para2
  ! contruct the cpot
  Call coupled_matrix0

  rh = R_min - dr 
  Call coupled_matrix(rh, CPOTH)
  Do  i = 1, Nlevel
    Do  j = 1, Nlevel
      A(i, j, -1) = 2*ReduceMass/Hbar**2 * CPOTH(i, j)
      if(i == j)Then
        A(i, j, -1) = A(i, j, -1) + 2*ReduceMass/Hbar**2 * (V(rh, L_i) - E)
      End If
    End Do
  End Do

  rh = R_min + dr / 2.0
  Call coupled_matrix(rh, CPOTH)


  Do ir = 0, R_iterat + 1
    r = R_min + dr * ir
    Call coupled_matrix(r, CPOT0)
    Do  i = 1, Nlevel
      Do  j = 1, Nlevel
        CPOT(i, j, ir) = CPOT0(i, j)
        A(i, j, ir) = 2*ReduceMass/Hbar**2 * CPOT0(i, j)
        if(i == j)Then
          A(i, j, ir) = A(i, j, ir) + 2*ReduceMass/Hbar**2 * (V(r, L_i) - E)
        End If
      End Do
    End Do
  End Do
  
  Do ir = 0, R_iterat + 1
    r = R_min + dr * ir
    
    Do  i = 1, Nlevel
      Do  j = 1, Nlevel
        B(i,j) = 0.0D0
        B(i,j) = dr**2 /SQRT(12.) * A(i, j, ir)
        if(i == j)B(i,j) = B(i,j) + sqrt(3.d0)
      End Do
    End Do

    Do i = 1, Nlevel
      Do j = 1, Nlevel
        B2(i,j) = 0.0D0
        Do ik = 1, Nlevel
          B2(i,j) = B2(i,j) + B(i,ik) * B(ik,j)
        End Do
        if(i == j)B2(i,j) = B2(i,j) - 1
      End Do
    End Do

    Do i = 1, Nlevel
      Do j = 1, Nlevel
          cc(i, j) = 0.0D0
          cc(i, j) = - dr**2 / 12.d0 * A(i, j, ir)
          If (i == j)cc(i, j) = cc(i, j) + 1.d0
      End Do
    End Do

    Do i = 1, Nlevel
      Do j = 1, Nlevel
        C(i,j) = 0.0D0
        Do ik = 1, Nlevel
          C(i,j) = C(i,j) + B2(i,ik) * cc(ik,j)
        End Do
      End Do
    End Do

    If(ir .lt. R_iterat)Then
      Do  i = 1, Nlevel
        Do  j = 1, Nlevel
          idx = Nlevel * ir + i
          jdx = Nlevel * ir + j
          H1(idx, jdx) = H1(idx, jdx) - C(i, j)
        End Do
      End Do
    End If

    If(ir == R_iterat)Then
      Do i = 1, Nlevel
        Do j = 1, Nlevel
          Hc(i,j) = -C(i,j)
        End Do
      End Do
    End If

  End Do


  Do ir = 0, (R_iterat-2)
    Do i = 1, Nlevel
      Do j = 1, Nlevel
          H1(ir * nlevel + i, ir * nlevel + nlevel + j) =  - dr**2 /12. * A(i, j, ir+1)
          H1(ir * nlevel + nlevel + i, ir * nlevel + j) =  - dr**2 /12. * A(i, j, ir)
          if(i == j)Then
            H1(ir * nlevel + i, ir * nlevel + nlevel + j) = H1(ir * nlevel + i, ir * nlevel + nlevel + j) + 1
            H1(ir * nlevel + nlevel + i, ir * nlevel + j) = H1(ir * nlevel + nlevel + i, ir * nlevel + j) + 1
          End If
      End Do
    End Do
  End Do
  
  
  Do i = 1, Nlevel
    Do j = 1, Nlevel
    Echannel(i) = E - V(R_min, L_i) - CPOT(i,i,0)
    ec = E - V(R_min, L_i) - CPOT(i,j,0)

    If(ec > 0)Then
      !k = acos(1. - ec/2. / t)/dr
      k = sqrt(2. * ReduceMass / Hbar**2 * ec)
      if(i == j)Then
        H1(i, j) = H1(i, j) + (1 - dr**2 /12. * A(i, j, -1)) * exp(ai * k * dr)
      else
        H1(i, j) = H1(i, j) + (  - dr**2 /12. * A(i, j, -1)) * exp(ai * k * dr)
      End If

    Else
      k = acos(1. - abs(ec)/2. / t)/dr
      H1(i, j) = H1(i, j) + ((- dr**2 /12. * A(i, j, -1)) * exp(k * dr))
      if(i == j)H1(i, j) = H1(i, j) + 1 * exp(k * dr)

    End If 

    End do
  End do

  Do i = 1, Nlevel*R_iterat
      Do j = 1, Nlevel*R_iterat
        H2(i,j) = H1(i,j)
      End Do
  End Do

  Do i = 1, Nlevel
    Do j = 1, Nlevel

      H2(Nlevel*R_iterat + i, Nlevel*R_iterat + j) = Hc(i, j)
      H2(Nlevel*R_iterat - Nlevel + i, Nlevel*R_iterat + j) = (- dr**2 /12. * A(i, j, R_iterat))
      H2(Nlevel*R_iterat + i, Nlevel*R_iterat - Nlevel + j) = (- dr**2 /12. * A(i, j, R_iterat-1))
      H2(Nlevel*R_iterat + i, Nlevel*R_iterat + Nlevel + j) = (- dr**2 /12. * A(i, j, R_iterat+1))

      If(i == j)Then
      H2(Nlevel*R_iterat - Nlevel + i, Nlevel*R_iterat + j) = H2(Nlevel*R_iterat - Nlevel + i, Nlevel*R_iterat + j) + 1
      H2(Nlevel*R_iterat + i, Nlevel*R_iterat - Nlevel + j) = H2(Nlevel*R_iterat + i, Nlevel*R_iterat - Nlevel + j) + 1
      H2(Nlevel*R_iterat + i, Nlevel*R_iterat + Nlevel + j) = H2(Nlevel*R_iterat + i, Nlevel*R_iterat + Nlevel + j) + 1
      End If

    End Do
  End Do

  Do i = 1, Nlevel*R_iterat
      dbpsi(i, i) = (1.0d0, 0.0d0)
  End Do


  Do i = 1, Nlevel
    ec = E - CPOT(i,i,R_iterat)
    !rho = acos(1. - ec/2. / t)/dr * (R_max)
    rho = sqrt(2.d0 * ReduceMass * ec) / Hbar * (R_max)
    eta = (Zpro * Ztar / 137.d0) * Sqrt(ReduceMass / (2.d0 * ec))
    Call myCOULFG(L_i*1.d0, eta, rho, fcw(L_i), fpcw(L_i), gcw(L_i), gpcw(L_i), iexp(L_i))
    dbpsi(Nlevel*R_iterat + i, Nlevel*R_iterat + i) = gcw(L_i) + ai*fcw(L_i)

    ec = E - CPOT(i,i,R_iterat+1)
    !rho = acos(1. - ec/2. / t)/dr  * (R_max + dr)
    rho = sqrt(2.d0 * ReduceMass * ec) / Hbar * (R_max + dr)
    eta = (Zpro * Ztar / 137.d0) * Sqrt(ReduceMass / (2.d0 * ec))
    Call myCOULFG(L_i*1.d0, eta, rho, fcw(L_i), fpcw(L_i), gcw(L_i), gpcw(L_i), iexp(L_i))
    dbpsi(Nlevel*R_iterat + i + Nlevel, Nlevel*R_iterat + i) = gcw(L_i) + ai*fcw(L_i)
  End Do

  do i = 1, Nx2
      dbpsi0(i) = dconjg(dbpsi(i, Nx1 - Nlevel + 1))
  end do

End Subroutine Prameterize_H

Subroutine DiscreteBasis(P)
  Use ccfull_initialization_mod

  Implicit None
  real(8), intent(out) :: P
  Real(8) ::  k, kk
  Integer ::  i, j
  complex*16 Hpsi(Nx1, Nx1)
  complex*16 Hpsi0(Nx1, 1)
  complex*16 cc(Nx1)
  Call MATMUL_CPLX(Hamilton, Ec_psi,  Hpsi, Nx1,  Nx2, Nx1)
  Call MATMUL_CPLX(Hamilton, Ec_psi0, Hpsi0, Nx1, Nx2, 1)
  call CSOLVE_EQUATION(Nx1,Hpsi,Hpsi0, cc)

  P = 0.

  Do i = 1, Nlevel
    if(Echannel(i) .gt. 0)then
      !k  = acos(1. - Echannel(i)/2. / t)/dr
      !kk = acos(1. - E/2. / t)/dr

      k=sqrt((2.d0*ReduceMass/Hbar**2*Echannel(i)))
      kk=sqrt(2.d0*ReduceMass/hbar**2*E)

      P = P+(abs(cc(i)))**2*k/kk
      !Write(*,*)'Discrete basis:', (abs(cc(i)))
    end if 
  end do

  Do i = 1, Nx2
    wf(i) = Ec_psi0(i)
  End Do

  Do i = 1, Nx1
    Do j = 1, Nx2
      wf(j) = wf(j) + cc(i) * Ec_psi(j, i)
    End Do
  End Do
End Subroutine

SUBROUTINE matmul_cplx(H2, PSI, HPSI, M, N, P)
    IMPLICIT NONE
    INTEGER M, N, P
    COMPLEX*16 H2(M, N), PSI(N, P), HPSI(M, P)

    COMPLEX*16 ALPHA, BETA
    CHARACTER*1 TRANSA, TRANSB

    ALPHA = (1.0D0, 0.0D0)
    BETA  = (0.0D0, 0.0D0)
    TRANSA = 'N'
    TRANSB = 'N'

    !CALL ZGEMM(TRANSA, TRANSB, M, P, N, ALPHA,H2,M,PSI, N,BETA, HPSI, M)
    !Return
    IF (P == 1) THEN
        ! 调用 ZGEMV，矩阵-向量乘法
        CALL ZGEMV(TRANSA, M, N, ALPHA, H2, M, PSI(:,1), 1, BETA, HPSI(:,1), 1)
    ELSE
        ! 调用 ZGEMM，矩阵-矩阵乘法
        CALL ZGEMM(TRANSA, TRANSB, M, P, N, ALPHA, H2, M, PSI, N, BETA, HPSI, M)
    END IF

    RETURN
End Subroutine

Subroutine matmul_cplx_sparse(H2, PSI, HPSI, M, N, P, Span)
    IMPLICIT NONE
    INTEGER M, N, P, Span
    COMPLEX*16 H2(M, N), PSI(N, P), HPSI(M, P)
    Real(8) ::  t00,t11,t22,t33,t44
    ! 稀疏矩阵相关变量
    INTEGER, ALLOCATABLE :: row_ptr(:), col_ind(:)
    COMPLEX*16, ALLOCATABLE :: values(:)
    INTEGER nnz, i, j, k, po, row_start, row_end
    INTEGER leftrange, rightrange
    COMPLEX*16 temp

    Call CPU_TIME(t00)
    ! 首先计算非零元素数量
    nnz = 0
    
    DO i = 1, M
        leftrange  = max(1, i - 2*Span)
        rightrange = min(N, i + 2*Span)
        DO j = leftrange, rightrange !1, N
            IF (ABS(H2(i,j)) > 1.0D-12) THEN
                nnz = nnz + 1
            END IF
        END DO
    END DO
    Call CPU_TIME(t11)
    ! 分配稀疏矩阵存储空间
    ALLOCATE(row_ptr(M+1))
    ALLOCATE(col_ind(nnz))
    ALLOCATE(values(nnz))

    ! 构建CSR格式
    nnz = 0
    row_ptr(1) = 1
    DO i = 1, M
        leftrange  = max(1, i - 2*Span)
        rightrange = min(N, i + 2*Span)
        DO j = leftrange, rightrange !1, N
        !DO j = 1, N
            IF (ABS(H2(i,j)) > 1.0D-12) THEN
                nnz = nnz + 1
                values(nnz) = H2(i,j)
                col_ind(nnz) = j
            END IF
        END DO
        row_ptr(i+1) = nnz + 1
    END DO
    Call CPU_TIME(t22)
    ! 使用CSR格式进行矩阵乘法
    HPSI = (0.0D0, 0.0D0)
    DO i = 1, M
        row_start = row_ptr(i)
        row_end = row_ptr(i+1) - 1
        DO k = row_start, row_end
            j = col_ind(k)
            temp = values(k)
            HPSI(i,:) = HPSI(i,:) + temp * PSI(j,:)
        END DO
    END DO


    Call CPU_TIME(t33)
    !Write(*,*)(t11-t00)/(t33-t00), (t22-t11)/(t33-t00), (t33-t22)/(t33-t00)
    !stop
    ! 释放内存
    DEALLOCATE(row_ptr, col_ind, values)

    RETURN
End Subroutine

SUBROUTINE CSOLVE_EQUATION(N, HPSI, HPSI0, CC)
      IMPLICIT NONE
      INTEGER N, INFO
      COMPLEX*16 HPSI(N,N), HPSI0(N), CC(N)

      INTEGER IPIV(N)
      COMPLEX*16 A(N,N), B(N)
      INTEGER I, J

      DO 20 J = 1, N
          DO 10 I = 1, N
              A(I,J) = HPSI(I,J)
   10     CONTINUE
   20 CONTINUE

      DO 30 I = 1, N
          B(I) = -HPSI0(I)
   30 CONTINUE

      CALL ZGESV(N, 1, A, N, IPIV, B, N, INFO)

      DO 40 I = 1, N
          CC(I) = B(I)
   40 CONTINUE

      RETURN
End Subroutine

Subroutine myCOULFG(angL,eta,rho,myf,mydf,myg,mydg,iv)
  double precision  angL,eta, rho, myf, mydf,myg,mydg
  integer iv,lmax 
  double precision XLMIN, XLMAX
  double precision XX, ETA1
  integer  MODE1, KFN, IFAIL
  double precision , allocatable, dimension(:):: FC,GC,FCP,GCP
  intent(in) angL
  intent(in) eta
  intent(in) rho
  intent(out) myf
  intent(out) mydf
  intent(out) myg
  intent(out) mydg
  intent(out) iv     
  XLMIN=angL
  XLMAX=angL
  lmax=IDINT(XLMAX)+1
  allocate (FC(lmax),GC(lmax),FCP(lmax),GCP(lmax))
  MODE1=1
  KFN=0 
  XX=rho
  ETA1=eta
  call COULFG(XX,ETA1,XLMIN,XLMAX,FC,GC,FCP,GCP,MODE1,KFN,IFAIL)
  myf=FC(lmax)
  mydf=FCP(lmax)
  myg=GC(lmax)
  mydg=GCP(lmax) 
  iv=IFAIL
  return 
End Subroutine

SUBROUTINE COULFG(XX, ETA1, XLMIN, XLMAX, FC, GC, FCP, GCP, MODE1, KFN, IFAIL)
    IMPLICIT NONE
    ! Argument type declarations
    DOUBLE PRECISION, INTENT(IN) :: XX, ETA1, XLMIN, XLMAX
    DOUBLE PRECISION, DIMENSION(*), INTENT(OUT) :: FC, GC, FCP, GCP
    INTEGER, INTENT(IN) :: MODE1, KFN
    INTEGER, INTENT(OUT) :: IFAIL

    ! Local variable declarations
    DOUBLE PRECISION :: ACCUR, ACC, ACC4, ACCH, ETA, GJWKB, PACCQ, X, XLM, E2MM1, DELL, XLL, XI, FCL, PK, PX, F, D, C
    DOUBLE PRECISION :: PK1, EK, RK2, TK, DF, FPL, XL, RL, EL, SL, FCL1, TA, WI, P, Q, AR, AI, BR, BI, DR, DI, DP, DQ
    DOUBLE PRECISION :: A, B, GAM, W, ALPHA, BETA, FCM, GCL, GPL, FJWKB, GCL1
    INTEGER :: MODE, NFP, NPQ, IEXP, M1, LXTRA, L1, LP, MAXL, L, LP1
    LOGICAL :: ETANE0, XLTURN

    ! Common block
    COMMON /STEE/ PACCQ, NFP, NPQ, IEXP, M1

    ! Constants
    DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, TEN2 = 1.0D2, ABORT = 2.0D4
    DOUBLE PRECISION, PARAMETER :: HALF = 0.5D0, TM30 = 1.0D-30, BIG = 1.0D+100
    DOUBLE PRECISION, PARAMETER :: RT2DPI = 0.7978845608286535587989211986876373D0

    ! Initialize accuracy and other variables
    ACCUR = 1.0D-16
    MODE = 1
    IF (MODE1 == 2 .OR. MODE1 == 3) MODE = MODE1
    IFAIL = 0
    IEXP = 1
    NPQ = 0
    ETA = ETA1
    GJWKB = ZERO
    PACCQ = ONE
    IF (KFN /= 0) ETA = ZERO
    ETANE0 = (ETA /= ZERO)
    ACC = ACCUR * 10.0D0
    ACC4 = ACC * TEN2 * TEN2
    ACCH = SQRT(ACC)

    ! Check for small XX
    IF (XX <= ACCH) THEN
        IFAIL = -1
        RETURN
    END IF

    X = XX
    XLM = XLMIN
    IF (KFN == 2) XLM = XLM - HALF
    IF (XLM <= -ONE .OR. XLMAX < XLMIN) THEN
        IFAIL = -2
        RETURN
    END IF

    E2MM1 = ETA*ETA + XLM*XLM + XLM
    XLTURN = (X*(X - TWO*ETA) < XLM*XLM + XLM)
    DELL = XLMAX - XLMIN + ACC
    LXTRA = INT(DELL)
    XLL = XLM + DBLE(LXTRA)

    M1 = MAX(INT(XLMIN + ACC), 0) + 1
    L1 = M1 + LXTRA

    XI = ONE / X
    FCL = ONE
    PK = XLL + ONE
    PX = PK + ABORT
    F = ETA/PK + PK*XI
    IF (ABS(F) < TM30) F = TM30
    D = ZERO
    C = F

    DO
        PK1 = PK + ONE
        EK = ETA / PK
        RK2 = ONE + EK*EK
        TK = (PK + PK1)*(XI + EK/PK1)
        D = TK - RK2 * D
        C = TK - RK2 / C
        IF (ABS(C) < TM30) C = TM30
        IF (ABS(D) < TM30) D = TM30
        D = ONE / D
        DF = D * C
        F = F * DF
        IF (D < ZERO) FCL = -FCL
        PK = PK1
        IF (PK > PX) THEN
            IFAIL = -3
            RETURN
        END IF
        IF (ABS(DF-ONE) < ACC) EXIT
    END DO

    NFP = PK - XLL - 1
    IF (LXTRA /= 0) THEN
        FCL = FCL / BIG
        FPL = FCL * F
        IF (MODE == 1) FCP(L1) = FPL
        FC(L1) = FCL
        XL = XLL
        RL = ONE
        EL = ZERO
        DO LP = 1, LXTRA
            IF (ETANE0) EL = ETA / XL
            IF (ETANE0) RL = SQRT(ONE + EL*EL)
            SL = EL + XL*XI
            L = L1 - LP
            FCL1 = (FCL * SL + FPL) / RL
            FPL = FCL1 * SL - FCL * RL
            FCL = FCL1
            FC(L) = FCL
            IF (MODE == 1) FCP(L) = FPL
            IF (MODE /= 3 .AND. ETANE0) GC(L+1) = RL
            IF (ABS(FCL) > BIG) THEN
                DO LP1 = L, M1+LXTRA
                    IF (MODE == 1) FCP(LP1) = FCP(LP1) * 1.0D-20
                    FC(LP1) = FC(LP1) * 1.0D-20
                END DO
                FCL = FC(L)
                FPL = FPL * 1.0D-20
            END IF
            XL = XL - ONE
        END DO
        IF (FCL == ZERO) FCL = ACC
        F = FPL / FCL
    END IF

    IF (XLTURN) CALL JWKB(X, ETA, MAX(XLM, ZERO), FJWKB, GJWKB, IEXP)
    IF (IEXP > 1 .OR. GJWKB > ONE/(ACCH*TEN2)) THEN
        W = FJWKB
        GAM = GJWKB * W
        P = F
        Q = ONE
    ELSE
        XLTURN = .FALSE.
        TA = TWO * ABORT
        PK = ZERO
        WI = ETA + ETA
        P = ZERO
        Q = ONE - ETA*XI
        AR = -E2MM1
        AI = ETA
        BR = TWO*(X - ETA)
        BI = TWO
        DR = BR/(BR*BR + BI*BI)
        DI = -BI/(BR*BR + BI*BI)
        DP = -XI*(AR*DI + AI*DR)
        DQ = XI*(AR*DR - AI*DI)
        DO
            P = P + DP
            Q = Q + DQ
            PK = PK + TWO
            AR = AR + PK
            AI = AI + WI
            BI = BI + TWO
            D = AR*DR - AI*DI + BR
            DI = AI*DR + AR*DI + BI
            C = ONE/(D*D + DI*DI)
            DR = C * D
            DI = -C * DI
            A = BR*DR - BI*DI - ONE
            B = BI*DR + BR*DI
            C = DP*A - DQ*B
            DQ = DP*B + DQ*A
            DP = C
            IF (PK > TA) THEN
                IFAIL = -4
                RETURN
            END IF
            IF (ABS(DP) + ABS(DQ) < (ABS(P) + ABS(Q)) * ACC) EXIT
        END DO
        NPQ = INT(PK / TWO)
        PACCQ = HALF * ACC / MIN(ABS(Q), ONE)
        IF (ABS(P) > ABS(Q)) PACCQ = PACCQ * ABS(P)
        GAM = (F - P) / Q
        IF (Q <= ACC4 * ABS(P)) THEN
            IFAIL = -5
            RETURN
        END IF
        W = ONE / SQRT((F - P)*GAM + Q)
    END IF

    ALPHA = ZERO
    IF (KFN == 1) ALPHA = XI
    IF (KFN == 2) ALPHA = XI * HALF
    BETA = ONE
    IF (KFN == 1) BETA = XI
    IF (KFN == 2) BETA = SQRT(XI) * RT2DPI
    FCM = SIGN(W, FCL) * BETA
    FC(M1) = FCM
    IF (MODE /= 3) THEN
        IF (.NOT. XLTURN) THEN
            GCL = FCM * GAM
        ELSE
            GCL = GJWKB * BETA
        END IF
        IF (KFN /= 0) GCL = -GCL
        GC(M1) = GCL
        GPL = GCL*(P - Q/GAM) - ALPHA*GCL
        IF (MODE /= 2) THEN
            GCP(M1) = GPL
            FCP(M1) = FCM*(F - ALPHA)
        END IF
    END IF

    IF (LXTRA == 0) RETURN

    W = BETA * W / ABS(FCL)
    MAXL = L1 - 1
    DO L = M1, MAXL
        IF (MODE /= 3) THEN
            XL = XL + ONE
            IF (ETANE0) EL = ETA / XL
            IF (ETANE0) RL = GC(L+1)
            SL = EL + XL*XI
            GCL1 = ((SL - ALPHA)*GCL - GPL) / RL
            IF (ABS(GCL1) > BIG) THEN
                IFAIL = L1 - L
                RETURN
            END IF
            GPL = RL*GCL - (SL + ALPHA)*GCL1
            GCL = GCL1
            GC(L+1) = GCL1
            IF (MODE /= 2) THEN
                GCP(L+1) = GPL
                FCP(L+1) = W*(FCP(L+1) - ALPHA*FC(L+1))
            END IF
        END IF
        FC(L+1) = W * FC(L+1)
    END DO

END SUBROUTINE COULFG

subroutine JWKB(XX, ETA1, XL, FJWKB, GJWKB, IEXP)
    implicit none
    ! ===== 输入输出参数 =====
    real(8), intent(in)  :: XX, ETA1, XL
    real(8), intent(out) :: FJWKB, GJWKB
    integer, intent(out) :: IEXP

    ! ===== 内部变量 =====
    real(8) :: ZERO, HALF, ONE, SIX, TEN
    real(8) :: DZ, RL35
    real(8) :: aloge, X, ETA, GH2, XLL1, HLL, HL, SL, RL2, GH
    real(8) :: PHI, PHI10

    ! ===== 常量赋值 =====
    ZERO = 0.0d0
    HALF = 0.5d0
    ONE  = 1.0d0
    SIX  = 6.0d0
    TEN  = 10.0d0
    DZ   = 0.0d0
    RL35 = 35.0d0

    aloge = log(TEN)

    ! ===== 核心计算 =====
    X   = XX
    ETA = ETA1

    GH2  = X * (2.0d0*ETA - X)
    XLL1 = max(XL*XL + XL, DZ)

    if (GH2 + XLL1 <= ZERO) then
        FJWKB = 0.0d0
        GJWKB = 0.0d0
        IEXP  = 0
        return
    end if

    HLL = XLL1 + SIX / RL35
    HL  = sqrt(HLL)
    SL  = ETA/HL + HL/X
    RL2 = ONE + ETA*ETA/HLL
    GH  = sqrt(GH2 + HLL) / X

    PHI = X*GH - HALF * ( HL*log((GH+SL)**2 / RL2) - log(GH) )
    if (ETA /= ZERO) PHI = PHI - ETA*atan2(X*GH, X - ETA)

    PHI10 = -PHI * aloge
    IEXP  = int(PHI10)

    if (IEXP > 70) then
        GJWKB = TEN**(PHI10 - dble(IEXP))
    else
        GJWKB = exp(-PHI)
        IEXP  = 0
    end if

    FJWKB = HALF / (GH * GJWKB)

end subroutine JWKB

Subroutine WHIT(HETA, R, XK, E, LL, F, FD, IE)
    implicit none
    ! ====== 输入输出参数 ======
    integer, intent(in) :: LL
    integer, intent(inout) :: IE
    real(8), intent(in) :: HETA, R, XK, E
    real(8), intent(out) :: F(LL+1), FD(LL+1)

    ! ====== 内部变量 ======
    integer :: L, LP1, LM, LMP1, IS, H
    integer :: I, M, M1, M2, JS, IFEQL
    real(8) :: fpmax, fpminl
    real(8) :: EE, AK, ETA, RHO, RHOA, PJE, A, B, C, D
    real(8) :: T(12), S(7)

    fpmax = 1.0d290
    L = LL + 1

    EE = -1.0d0
    AK = XK
    ETA = HETA
    LP1 = L + 1
    RHO = AK * R

    S(:) = 0.0d0

    ! ====== 设置 LM ======
    if (L <= 50) then
        LM = 60
    else
        LM = L + 10
    end if
    LMP1 = LM + 1
    IS = 7

    PJE = 30.0d0 * RHO + 1.0d0
    H = max(int(PJE), 4)
    H = RHO / H

    RHOA = 10.0d0 * (ETA + 1.0d0)
    if (RHOA <= RHO) then
        IFEQL = 1
        RHOA = RHO
    else
        IFEQL = 0
    end if

    PJE = RHOA / H + 0.5d0
    RHOA = H * int(PJE)

    if (IFEQL == 1) then
        if (RHOA <= RHO + 1.5d0 * H) then
            RHOA = RHO + 2.0d0 * H
        end if
    end if

    if (EE < 0.0d0) then
        ! ----------- 级数展开部分 ----------
        C = 1.0d0 / RHOA
        A = 1.0d0
        B = 1.0d0 - C * ETA
        F(1) = A
        FD(1) = B
        do M = 1, 26
            D = 0.5d0 * (ETA + dble(M - 1)) * (ETA + dble(M)) * C / dble(M)
            A = -A * D
            B = -B * D - A * C
            F(1) = F(1) + A
            FD(1) = FD(1) + B
        end do

        A = -ETA * log(2.0d0 * RHOA) - RHOA
        fpminl = -log(fpmax)
        if (IE == 0 .and. A < fpminl) IE = int(fpminl - A)
        A = exp(A + dble(IE))

        F(1) = A * F(1)
        FD(1) = A * FD(1) * (-1.0d0 - 2.0d0 * ETA / RHOA)

        if (IFEQL == 1) then
            S(IS) = F(1)
            if (IS == 7) then
                IS = 6
                RHOA = RHOA + H
                ! 再做一次展开
                C = 1.0d0 / RHOA
                A = 1.0d0
                B = 1.0d0 - C * ETA
                F(1) = A
                FD(1) = B
                do M = 1, 26
                    D = 0.5d0 * (ETA + dble(M - 1)) * (ETA + dble(M)) * C / dble(M)
                    A = -A * D
                    B = -B * D - A * C
                    F(1) = F(1) + A
                    FD(1) = FD(1) + B
                end do
                A = -ETA * log(2.0d0 * RHOA) - RHOA
                fpminl = -log(fpmax)
                if (IE == 0 .and. A < fpminl) IE = int(fpminl - A)
                A = exp(A + dble(IE))
                F(1) = A * F(1)
                FD(1) = A * FD(1) * (-1.0d0 - 2.0d0 * ETA / RHOA)
            end if
        else
            F(1) = T(1)
            FD(1) = T(2)
        end if

    else
        stop "WHIT1: EE >= 0 not implemented"
    end if

    ! ----------- 递推关系 ----------
    C = 1.0d0 / RHO
    do M = 1, L-1
        A = ETA / dble(M)
        B = A + C * dble(M)
        F(M+1) = (B * F(M) - FD(M)) / (A + 1.0d0)
        FD(M+1) = (A - 1.0d0) * F(M) - B * F(M+1)
    end do

    do M = 1, L
        FD(M) = AK * FD(M)
    end do

    return
End Subroutine

Subroutine Ec_training()
  Use ccfull_initialization_mod
  Implicit None
  Integer ::  Iparameter, i, ipara1, ipara2, j
  
  Real(8) :: P

  wf_base = (0.0D0, 0.0D0)
  Iparameter = 0
  Do ipara1 = 1, num_traning_set1
    Do ipara2 = 1, num_traning_set2

      Iparameter = Iparameter + 1
      Call Prameterize_H(training_set1(ipara1), training_set2(ipara2), Hamilton, Ec_psi, Ec_psi0)
      !Call DiscreteBasis(P)
      Call Numerov_Ec(P)
      wf_base(:, Iparameter) = wf
      Do j = Nx_cc + 1, Nx2
        wf_base(j, Iparameter) = (0.0D0, 0.0D0)
      End Do

    End Do
  End Do
End Subroutine

Subroutine Ec_emulator(P_ec)
  Use ccfull_initialization_mod
  Implicit None
  Integer ::  Iparameter, i, j
  Real(8) ::  P_ec, k, kk
  Real(8) ::  t00,t11,t22,t33,t44

  !Write(*,*)'emulator',Nx1, Nx2, NumEmulator+Nlevel
  Call CPU_TIME(t00)
 !construct the total wave function
  Do i = 1, 2*Nlevel
    Do j = 1, Nlevel
      wf_base(Nx_cc+i, NumEmulator+j) = Ec_psi(Nx_cc+i, Nx_cc+j)
    End Do 
  End Do
  Call CPU_TIME(t11)
  !|H|psi>
  !Call matmul_cplx(Hamilton, wf_base, wf_base2, Nx1, Nx2, NumEmulator+Nlevel)
  Call matmul_cplx_sparse(Hamilton, wf_base, wf_base2, Nx1, Nx2, NumEmulator+Nlevel, Nlevel)
  Call CPU_TIME(t22)
  !|H|psi0>
  !Call matmul_cplx(Hamilton, Ec_psi0, Ec_psi0_2, Nx1, Nx2, 1)
  Call matmul_cplx_sparse(Hamilton, Ec_psi0, Ec_psi0_2, Nx1, Nx2, 1, Nlevel)
  Call CPU_TIME(t33)
  !Write(*,*)(t11-t00)/(t33-t00), (t22-t11)/(t33-t00), (t33-t22)/(t33-t00)
  !|psi>  ---> <psi|
  Do i = 1, Nx1
    Do j = 1, NumEmulator+Nlevel
        wf_base3_dagger(j, i) = dconjg(wf_base2(i, j))
    End Do
  End Do

  !<psi|H|psi>
  Call matmul_cplx(wf_base3_dagger, wf_base2, cc_ec, NumEmulator+Nlevel, Nx1, NumEmulator+Nlevel)


  !<psi|H|psi0>
  Call matmul_cplx(wf_base3_dagger, Ec_psi0_2, dd_ec, NumEmulator+Nlevel, Nx1, 1)
  Call matinv_lapack(NumEmulator+Nlevel,cc_ec,cc_ec_inv)

  Call matmul_cplx(-cc_ec_inv, dd_ec, dd2_ec, NumEmulator+Nlevel, NumEmulator+Nlevel, 1)
  

  
  !Write(*,*)'emulator',Nx1, Nx2, NumEmulator+Nlevel
  a_tra = (0.D0, 0.D0)
  Do i = 1, NumEmulator
    Do j = 1, Nlevel
      a_tra(j) = a_tra(j) + dd2_ec(i)*wf_base(j,i)
    End Do
  End Do

  P_ec = 0.

  Do i = 1, Nlevel
    If(Echannel(i)<0)Exit
    k = sqrt((2.d0* ReduceMass/Hbar**2 * Echannel(i)))
    kk= sqrt(2.d0 * ReduceMass/hbar**2 * E)
    P_ec = P_ec + (abs(a_tra(i)))**2*k/kk
  End Do

End Subroutine


Program CCFull_Main
  Use ccfull_initialization_mod
  Implicit None
  Integer ::  i, j, ir, ie, il
  Integer ::  Iestep, Iparameter, ipara1, ipara2
  Real(8) ::  P, P_ec, s0
  Real(8) ::  V
  Real(8) ::  t0,t1,t2,t3,t4, t00, t11, t22
  Real(8) ::  dsigma_t, dsigma_ec
 
  Real(8), allocatable :: P_t_0(:), sigma_t(:)
  Real(8), allocatable :: P_Ec_0(:), sigma_Ec(:)
  
  character(len=20)   :: str_Beta2T, str_Beta4T, str_Apro, str_Atar
  character(len=20)   :: SystemName
  character(len=100)  :: Filename

  external V

  Call Initialize_CCFull()
  Write(*,*)'Nlevle = ', Nlevel


  Allocate(P_t_0(Num_testing), sigma_t(Num_testing))
  Allocate(P_Ec_0(Num_testing), sigma_Ec(Num_testing))

  Iestep = int((Emax - Emin)/dE)

  Do ie = 0, Iestep

    E = Emin + dE * ie

    P_t_0 = 0.0D0
    sigma_t = 0.0D0
    P_Ec_0 = 0.0D0
    sigma_Ec = 0.0D0

    Call CPU_TIME(t0)
    Do il = 0, 200
      Call CPU_TIME(t00)
      L_i = il
      
      Call PotShape(L_i, R_barrier_l, V_barrier_l, curv_l, R_bottom_l, V_bottom_l)
      If( V(R_bottom_l, L_i) > E .or. R_bottom_l < 0 )Then
          P = 0.
          P_ec = 0.
          Cycle  
      End If

      !Builde the emulator

      Call Ec_training()
      Call CPU_TIME(t11)
     
      !Write(*,*)E, L_i, 'traning time:', t11 - t00
      !Use emulator to calculate
      Iparameter = 0
      
      Do ipara1 = 1, num_testing_set1
        Do ipara2 = 1, num_testing_set2
          Iparameter = Iparameter + 1
          !s0 = sigma_t(Iparameter)
          Call Prameterize_H(testing_set1(ipara1), testing_set2(ipara2),  Hamilton, Ec_psi, Ec_psi0)

          Call CPU_TIME(t2)
          !Call DiscreteBasis(P)
          Call Numerov_Ec(P) !numerov
          dsigma_t = (2.0D0 * L_i + 1.0D0) * P * pi * Hbar**2 / (2.0D0 * ReduceMass * E) * 10.0D0
          
          sigma_t(Iparameter) = sigma_t(Iparameter) + dsigma_t
                              
          
          Call CPU_TIME(t3)


          Call Ec_emulator(P_ec)
          dsigma_ec = (2.0D0 * L_i + 1.0D0) * P_ec * pi * Hbar**2 / (2.0D0 * ReduceMass * E) * 10.0D0
          sigma_ec(Iparameter) = sigma_ec(Iparameter) + dsigma_ec
                              


          Call CPU_TIME(t4)

          !If(L_i == 0 .and. Iparameter == 1)Write(*,*)E, L_i, 'Error:', abs(P-P_ec)/P, 'accelarator ratio:', (t3-t2)/(t4-t3)
          !If(L_i == 0 .and. Iparameter == 1)Write(*,*)E, L_i,  (t3-t2), (t4-t3)
          Write(999,*)E, L_i, Iparameter, P, P_ec, abs(P-P_ec)/P, (t3-t2), (t4-t3)

          If (L_i == 0) Then
              P_t_0(Iparameter) = P
              P_Ec_0(Iparameter) = P_ec
          End If
          !If(sigma_t(Iparameter) - s0 < s0*1.d-4)Exit

        End Do 
      End Do
      
    Call CPU_TIME(t22)
    WRITE(*, '("E = ", F8.2, ", L_i = ", I0, " Training time = ", F8.2, " Prediction time = ", F8.2)')E, L_i, t11-t00, t22-t11
    

    End Do

    Call CPU_TIME(t1)
    Write(*,*)E, 'Each Energy time:', t1-t0, 's'


    !Output Result
    Iparameter = 0
    Do ipara1 = 1, num_testing_set1
      Do ipara2 = 1, num_testing_set2
        Iparameter = Iparameter + 1
        write(str_Apro, '(I0)') Int(Apro)
        write(str_Atar, '(I0)') Int(Atar)
        write(str_Beta2T, '(F6.3)') testing_set1(ipara1)
        write(str_Beta4T, '(F6.3)') testing_set2(ipara2)
        write(SystemName, '(A3,A,A1,A3,A2)')TRIM(ADJUSTL(str_Apro)), TRIM(ADJUSTL(Element(Zpro))), &
                                          '+', TRIM(ADJUSTL(str_Atar)), TRIM(ADJUSTL(Element(Ztar)))
        Filename = "CCfull_sigma_" // TRIM(ADJUSTL(SystemName)) // "_Ec_" // &
                  "Beta2T_" // TRIM(ADJUSTL(str_Beta2T)) // "_" // &
                  "Beta4T_" // TRIM(ADJUSTL(str_Beta4T)) // ".dat"
        Open(20, File=Filename,  Status = 'unknown',position='append',action='write')
        write(20,'(F8.3,5E15.5)')E, P_t_0(Iparameter), P_Ec_0(Iparameter), sigma_t(Iparameter), sigma_ec(Iparameter),&
                                                     abs(sigma_t(Iparameter)-sigma_ec(Iparameter))/sigma_t(Iparameter)
        Close(20)

      End Do
    End Do

  End Do

  Print *, "Program completed successfully."

End Program CCFull_Main


      !Call Prameterize_H(0.205d0,Hamilton, Ec_psi, Ec_psi0)
      !Call Numerov_Ec(P1)
      !Do i = 0, R_iterat-1
      !  Write(50,*)R_min+i*dr, Real(wf(nlevel*i+1)),  AIMAG(wf(nlevel*i+1)), Real(wf(nlevel*i+2)),  AIMAG(wf(nlevel*i+2))
      !End Do 
      !Call DiscreteBasis(P2)
      !Do i = 0, R_iterat-1
      !  Write(40,*)R_min+i*dr, Real(wf(nlevel*i+1)),  AIMAG(wf(nlevel*i+1)), Real(wf(nlevel*i+2)),  AIMAG(wf(nlevel*i+2))
      !End Do 
    
      !Write(*,*)E, L_i, P1, P2, abs(P1-P2)/P1
      !stop
      !if (E == 60 .and. L_i == 0)Then
      !  Do i = 0, R_iterat+1
      !    r = R_min + dr * i
      !    write(40, '(F10.5, 1X, 100(1X, F12.5, 1X,  F12.5))') r, (REAL(Wf(i,j)), AIMAG(Wf(i,j)), j=1, Nlevel)
      !  End Do
      !End If


!       55.00000    0.97449E-02        5.87031
!       56.00000        0.05489        5.94333
!       57.00000        0.28583        6.05134
!       58.00000        1.36500        6.19272
!       59.00000        5.84375        6.40451
!       60.00000       20.59856        6.86092
!       61.00000       52.14435        7.81887
!       62.00000       94.62477        9.18913
!       63.00000      139.58988       10.65032
!       64.00000      185.55960       11.98384
!       65.00000      234.04527       13.13045
!       66.00000      283.93527       14.18620
!       67.00000      333.25670       15.21143
!       68.00000      381.19965       16.20589
!       69.00000      427.60179       17.16365
!       70.00000      472.46037       18.08247



! E =   55.000 MeV,   Sigma =     0.97448E-02 mb, <L> =  5.87031 hbar, <P_0> =     0.22629E-03
! E =   56.000 MeV,   Sigma =         0.05489 mb, <L> =  5.94333 hbar, <P_0> =     0.12592E-02
! E =   57.000 MeV,   Sigma =         0.28583 mb, <L> =  6.05134 hbar, <P_0> =     0.64144E-02
! E =   58.000 MeV,   Sigma =         1.36500 mb, <L> =  6.19272 hbar, <P_0> =     0.29721E-01
! E =   59.000 MeV,   Sigma =         5.84375 mb, <L> =  6.40451 hbar, <P_0> =     0.11981E+00
! E =   60.000 MeV,   Sigma =        20.59856 mb, <L> =  6.86092 hbar, <P_0> =     0.35491E+00
! E =   61.000 MeV,   Sigma =        52.14434 mb, <L> =  7.81887 hbar, <P_0> =     0.62839E+00
! E =   62.000 MeV,   Sigma =        94.62476 mb, <L> =  9.18913 hbar, <P_0> =     0.75138E+00
! E =   63.000 MeV,   Sigma =       139.58987 mb, <L> = 10.65032 hbar, <P_0> =     0.79530E+00
! E =   64.000 MeV,   Sigma =       185.55959 mb, <L> = 11.98384 hbar, <P_0> =     0.86342E+00
! E =   65.000 MeV,   Sigma =       234.04527 mb, <L> = 13.13045 hbar, <P_0> =     0.94486E+00
! E =   66.000 MeV,   Sigma =       283.93526 mb, <L> = 14.18620 hbar, <P_0> =     0.98364E+00
! E =   67.000 MeV,   Sigma =       333.26115 mb, <L> = 15.21129 hbar, <P_0> =     0.99507E+00
! E =   68.000 MeV,   Sigma =       381.21016 mb, <L> = 16.20563 hbar, <P_0> =     0.99820E+00
! E =   69.000 MeV,   Sigma =       427.61803 mb, <L> = 17.16333 hbar, <P_0> =     0.99921E+00
! E =   70.000 MeV,   Sigma =       472.48080 mb, <L> = 18.08211 hbar, <P_0> =     0.99962E+00


! E =   55.000 MeV,   Sigma =     0.99270E-02 mb, <L> =  5.86833 hbar, <P_0> =     0.22967E-03
! E =   56.000 MeV,   Sigma =         0.05545 mb, <L> =  5.95481 hbar, <P_0> =     0.12643E-02
! E =   57.000 MeV,   Sigma =         0.28575 mb, <L> =  6.08041 hbar, <P_0> =     0.63493E-02
! E =   58.000 MeV,   Sigma =         1.35843 mb, <L> =  6.21020 hbar, <P_0> =     0.29472E-01
! E =   59.000 MeV,   Sigma =         5.80359 mb, <L> =  6.40377 hbar, <P_0> =     0.11981E+00
! E =   60.000 MeV,   Sigma =        20.61289 mb, <L> =  6.84429 hbar, <P_0> =     0.35802E+00
! E =   61.000 MeV,   Sigma =        52.45950 mb, <L> =  7.79227 hbar, <P_0> =     0.63535E+00
! E =   62.000 MeV,   Sigma =        95.05488 mb, <L> =  9.17641 hbar, <P_0> =     0.75517E+00
! E =   63.000 MeV,   Sigma =       140.09738 mb, <L> = 10.64040 hbar, <P_0> =     0.79767E+00
! E =   64.000 MeV,   Sigma =       186.08861 mb, <L> = 11.97409 hbar, <P_0> =     0.86742E+00
! E =   65.000 MeV,   Sigma =       234.57189 mb, <L> = 13.12039 hbar, <P_0> =     0.94765E+00
! E =   66.000 MeV,   Sigma =       284.44688 mb, <L> = 14.17473 hbar, <P_0> =     0.98670E+00
! E =   67.000 MeV,   Sigma =       333.72636 mb, <L> = 15.19986 hbar, <P_0> =     0.99845E+00
! E =   68.000 MeV,   Sigma =       381.62127 mb, <L> = 16.19269 hbar, <P_0> =     0.10017E+01
! E =   69.000 MeV,   Sigma =       427.95219 mb, <L> = 17.14905 hbar, <P_0> =     0.10030E+01
! E =   70.000 MeV,   Sigma =       472.74120 mb, <L> = 18.06709 hbar, <P_0> =     0.10032E+01


! nlevel=           3
!       50.00000    0.13216E-03        5.67434
!       51.00000    0.11981E-02        5.70974
!       52.00000        0.01004        5.78231
!       53.00000        0.07784        5.86061
!       54.00000        0.54862        5.96724
!       55.00000        3.15947        6.33560
!       56.00000       10.99611        7.47440
!       57.00000       22.48305        9.18772
!       58.00000       35.44297       10.80500
!       59.00000       52.44661       11.78712
!       60.00000       78.76923       12.12881
!       61.00000      115.61499       12.41005
!       62.00000      163.07554       12.77286
!       63.00000      215.39252       13.42451
!       64.00000      267.38084       14.28697
!       65.00000      317.69128       15.22206
!       66.00000      366.19428       16.16368
!       67.00000      412.91931       17.08422
!       68.00000      457.93722       17.97450
!       69.00000      501.29524       18.83117
!       70.00000      543.04197       19.65412
!       71.00000      583.22750       20.44471
!       72.00000      621.93634       21.20563
!       73.00000      659.23615       21.93896
!       74.00000      695.12471       22.64453
!       75.00000      729.90346       23.33079