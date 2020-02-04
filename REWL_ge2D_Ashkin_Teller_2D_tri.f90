!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!Definitions of kind value to integer variables
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
MODULE integer_kind
  IMPLICIT NONE

  INTEGER, PARAMETER :: K15=selected_int_kind(9)!selected_int_kind(15)
  INTEGER, PARAMETER :: K4B=selected_int_kind(9)
END MODULE integer_kind
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!Variables for the pseudo-random numbers generator
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
MODULE raset1
  USE integer_kind
  IMPLICIT NONE

  REAL(8), SAVE :: u(97), cc, cd, cm
  INTEGER(K15), SAVE :: i97, j97
END MODULE raset1
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!Variables for the Wang-Landau Algorithm
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
MODULE WL_variables
  USE integer_kind
  IMPLICIT NONE
  INTEGER(K4B), PARAMETER :: mcs_bin = 1000, mcs_bin_term = 10, ireplica=1000, max_neighbors = 4
  REAL(8), PARAMETER :: lnf_min=1.0D-9,accept_min=0.6D0,p=0.7,overlap_ratio=0.75D0,&
                      & dE1=8.0D0, dE2=4.0D0, dE1_inv=1.0D0/dE1, dE2_inv=1.0D0/dE2
  INTEGER(K4B) :: ierr,numprocs,myid,mcs,n_mcs,energy_comm,energy_id,numranks,iWalker,&
                & numbers_count, nWalkers, nWindows, nEbins, cont, kWindow, mWindows
  INTEGER(k15) :: delta_range1, delta_range2, rangeE1_min, rangeE1_max, rangeE2_min, rangeE2_max
  INTEGER(K4B) :: iseed1,iseed2,accept_counter,exchange_counter, number_windows, even_windows
  INTEGER(K4B),ALLOCATABLE,DIMENSION(:) :: vizXdir,vizXesq,vizYdir,vizYesq,counter_s,counter_t,&
                                          & my_neighbor_id,position_aux
  INTEGER(K4B),ALLOCATABLE,DIMENSION(:) :: number_neighbor_windows,neighbor_windows
  REAL(8) :: pi,mcs_inv,accept_ratio,lnf,dU,E1_min,E1_max,E2_min,E2_max,phi_optimum,dble_nWalkers_1,&
            & replica_ratio,Emin_N, Emax_N = 0.0D0
  INTEGER(K15),ALLOCATABLE,DIMENSION(:) :: he
  REAL(8),ALLOCATABLE,DIMENSION(:) :: lng, random_numbers
  LOGICAL,ALLOCATABLE,DIMENSION(:) :: my_neighbor
  LOGICAL :: flat, accept

END MODULE WL_variables
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!Variables for the spin model
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
MODULE spins_variables
  USE integer_kind
  IMPLICIT NONE
  INTEGER(K4B), PARAMETER :: L=16, L2 = L*L
  REAL(8),PARAMETER :: L2_dble=DBLE(L2)-1.0D0, L2_inv = 1.0D0/DBLE(L2)
  INTEGER(K4B) :: delta_sx, delta_tx
  REAL(8) :: Mx,My,Mz,theta_opt,Sx_new,Sy_new,Sz_new,delta_U1,delta_U2!, A, B
  !New variable type for spin configuration
  TYPE spins_config
    REAL(8) :: U1, U2
    INTEGER(K4B),DIMENSION(1:L2) :: Sx, Tx
  END TYPE spins_config
  TYPE (spins_config) :: current, trial
  INTEGER(K15),DIMENSION(1:3) :: E_current, E_trial

END MODULE spins_variables
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!Variables for the MPI library
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
MODULE MPI_new_variables
USE integer_kind
  IMPLICIT NONE
  include "mpif.h"

  INTEGER(K4B), PARAMETER  :: cont_blocks=4
  INTEGER :: mpistatus(MPI_STATUS_SIZE), MPI_config_spins
  INTEGER,DIMENSION(:) :: block_length(1:cont_blocks), address(1:cont_blocks+1), typelist(1:cont_blocks)
  INTEGER(kind=MPI_ADDRESS_KIND),DIMENSION(1:cont_blocks) :: displacement

END MODULE MPI_new_variables
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!##################################################################################################!
!Wang-Landau Method for xy-model
!##################################################################################################!
PROGRAM WangLandau_xy_model
  USE integer_kind
  USE WL_variables
  USE spins_variables
  USE MPI_new_variables
  IMPLICIT NONE
  INTEGER(K15) :: i,j,k
  CHARACTER(50) :: lnge!,xmakemol_file
913 FORMAT('oREWL_AshkinTeller_2D_L=',I2,'seed',I2,'_a.out')
!914 FORMAT('config_xy_2D_L='I3'seed'I2'.xyz')

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
!MPI initialization
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  CALL MPI_INIT(ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
!Parameters and Variables initialization
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
!  WRITE(xmakemol_file,914) L,myid+1
!  OPEN(UNIT=20, FILE=xmakemol_file)
  CALL initial()
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
!Wang-Landau Algorithm
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  he = 0; lng = 1.0D0; lnf = 1.0D0 !Initialization of the Variables
  DO WHILE (lnf .GT. lnf_min) 
    CALL sweep()    ! Perform mcs Monte Carlo sweeps (mcs = mcs_bin*nEbins)
    CALL flathist() !Check if histogram is flat
  END DO
  !. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .!
  !Simulation Output (log of the density of states)
  !. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .!
  WRITE(lnge,913) L,myid+1
  OPEN(UNIT=10, FILE=lnge)
  DO i=RangeE1_min,RangeE1_max
    DO j=RangeE2_min,RangeE2_max
      k = delta_range2*i + j
      IF(lng(k) > 1.0D0) WRITE(10,*) DBLE(i)*dE1, DBLE(j)*dE2, lng(k)
    END DO
  END DO
  CLOSE(10)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
!MPI finalization
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
DEALLOCATE(he,lng,vizXdir,vizXesq,vizYdir,vizYesq,counter_s,counter_t,&
           & position_aux,neighbor_windows,number_neighbor_windows)
  CALL MPI_FINALIZE(ierr)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  STOP
!##################################################################################################!
CONTAINS
!------------------------------------------------------------------------------!
!Initialization (initial configuration, constants, table of energy and space Ex)
!------------------------------------------------------------------------------!
SUBROUTINE initial()
USE integer_kind
USE WL_variables
USE spins_variables
USE MPI_new_variables
IMPLICIT NONE
  INTEGER(K4B) :: i,j,k,ii,kk!,ii
  REAL(8) :: aux1,aux2,aux4,DeltaEp1,DeltaEp2!,norm,phi
  INTEGER :: count, count_rate, count_max, iWindows, iaux, jaux
  REAL(8), ALLOCATABLE,DIMENSION(:) :: E1_min_aux,E1_max_aux,E2_min_aux,E2_max_aux
  INTEGER(k4B),DIMENSION(12) :: iseed_aux
  INTEGER:: rangeE1_min_aux, rangeE1_max_aux, rangeE2_min_aux, rangeE2_max_aux
  INTEGER(k4B),ALLOCATABLE,DIMENSION(:) :: jWindow
  REAL :: seed_aux
  LOGICAL :: found, found_core

!  REAL(8),ALLOCATABLE :: lngE_med(:)
!  LOGICAL :: found_E, found_E_core
!  CHARACTER(30) :: estimative
!666 FORMAT('Estimative_xy_2D_L=',I2,'.in')

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
!Reading input file
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  IF (myid==0) THEN
    OPEN(UNIT=10, FILE='input_REWL_AshkinTeller_2D.dat')
    READ(10,*) nWalkers, mWindows
    CLOSE(10)
    number_windows = 0
    DO i=0,mWindows-1
      nWindows = 1 + 2*i
      number_windows = number_windows + nWindows
    END DO
    PRINT*,'number of windows of energy',nWindows,'X',mWindows, '=>',&
           & number_windows
    PRINT*,'number of walkers',nWalkers
  END IF
  CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
  IF(myid==0) PRINT*,'The input file was read.'
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  !Broadcasting the initial data
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  CALL MPI_Bcast(nWalkers,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_Bcast(nWindows,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_Bcast(mWindows,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  CALL MPI_Bcast(number_windows,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  IF(myid==0)  PRINT*,'The data was broadcasted to all processors.'
  CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  !Checking if the number of processors was passed correctL.
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  IF(numprocs .NE. number_windows*nWalkers) THEN
    IF (myid==0) THEN
      PRINT*,'PROGRAM ABORT!'
      PRINT*,'wrong number of processors',numprocs,'it should be',number_windows*nWalkers
    END IF
    CALL MPI_FINALIZE(ierr)
    STOP
  END IF
  dble_nWalkers_1 = DBLE(nWalkers-1)
  IF(myid==0) PRINT*,'The number of processors was checked.'
  CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  !New communicator (For walkers inside a given energy window)
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  kWindow = myid/nWalkers
  iWalker = MOD(myid,nWalkers)
  CALL MPI_Comm_split(MPI_COMM_WORLD,kWindow,iWalker,energy_comm,ierr)
  CALL MPI_COMM_SIZE(energy_comm,numranks,ierr)
  CALL MPI_COMM_RANK(energy_comm,energy_id,ierr)
  IF(myid==0) PRINT*,'The new communicator was created.'
  CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  !New MPI Data Type (Configuration of the System State)
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  typelist(1) = MPI_DOUBLE_PRECISION
  typelist(2) = MPI_DOUBLE_PRECISION
  typelist(3) = MPI_INTEGER
  typelist(4) = MPI_INTEGER
  block_length(1) = 1
  block_length(2) = 1
  block_length(3) = L2
  block_length(4) = L2
  CALL MPI_ADDRESS(current,address(1),ierr)
  CALL MPI_ADDRESS(current%U1,address(2),ierr)
  CALL MPI_ADDRESS(current%U2,address(3),ierr)
  CALL MPI_ADDRESS(current%Sx,address(4),ierr)
  CALL MPI_ADDRESS(current%Tx,address(5),ierr)
  DO i=1,cont_blocks
    displacement(i) = address(i+1) - address(1)
  END DO
  CALL MPI_TYPE_CREATE_STRUCT(cont_blocks,block_length,displacement,typelist,MPI_config_spins,ierr)
  CALL MPI_TYPE_COMMIT(MPI_config_spins,ierr)
  IF(myid==0)  PRINT*,'The new mpi data type was created.'
  CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  !Random Number Generator inicialization
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  DO i=1,12
    DO j=0,myid+i
      CALL SYSTEM_CLOCK(count, count_rate, count_max)
      iseed_aux(i) = count
    END DO
    iseed_aux(i) = iseed_aux(i)*myid*i + iseed_aux(i)*i
  END DO
  CALL random_seed(PUT=iseed_aux)
  DO i=1,numprocs + 10*myid
    CALL random_number(seed_aux)
    iseed1 = NINT(31328.D0*seed_aux)
    CALL random_number(seed_aux)
    iseed2 = NINT(30081.D0*seed_aux)
  END DO
  PRINT*, 'Processor ' ,myid, '. Seeds for the Marsaglia random number generator', iseed1,iseed2
  CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
  CALL rmarin(iseed1,iseed2)
  IF (myid == 0) PRINT*, 'Gerador de números aleatórios inicializado '
  CALL MPI_Barrier(MPI_COMM_WORLD,ierr)

  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  !Alocando variaveis e tabelando parametros
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
!  IF (myid ==0) CALL CPU_TIME(t_inic)    !Inicializa a contagem do tempo de computacao
  !Defining the lattice (simple cubic (SC) lattice).
  !A surface with periodic boundary condition in xy-plane.
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  ALLOCATE (vizXdir(1:L2),vizXesq(1:L2),vizYdir(1:L2),vizYesq(1:2*L2),&
           & counter_s(1:L2),counter_t(1:L2))
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  !Tabelando os primeiros vizinhos
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  cont = 0
  DO i=1,L
    DO j=1,L
      cont = cont + 1
      counter_s(cont) = cont
      counter_t(cont) = cont
      !Define a posicao dos vizinhos na direcao x
      vizXdir(cont) = cont + 1
      vizXesq(cont) = cont - 1
      IF(j .EQ. 1) vizXesq(cont) = cont - 1 + L
      IF(j .EQ. L) vizXdir(cont) = cont + 1 - L
      !Define a posicao dos vizinhos na direcao y
      vizYdir(cont) = cont + L
      vizYesq(cont) = cont - L
      IF(i .EQ. 1) vizYesq(cont) = cont - L + L2
      IF(i .EQ. L) vizYdir(cont) = cont + L - L2
    END DO
  END DO
  current%sx = 1
  current%tx = 1
  CALL energy(current,E_current)

  IF(myid==0)  PRINT*,'A square lattice was built.'
  CALL MPI_Barrier(MPI_COMM_WORLD,ierr)

  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  !Defining some variables and allocating memory:
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  E1_min = -4.0D0*L2       !Minimum energy
  E1_max = 4.0D0*L2        !Maximum energy
  E2_min = -2.0D0*L2       !Minimum energy
  E2_max = 2.0D0*L2        !Maximum energy

  rangeE1_min = NINT(E1_min*dE1_inv)
  rangeE1_max = NINT(E1_max*dE1_inv)
  rangeE2_min = NINT(E2_min*dE2_inv)
  rangeE2_max = NINT(E2_max*dE2_inv)
  delta_range1 = rangeE1_max-rangeE1_min + 1
  delta_range2 = rangeE2_max-rangeE2_min + 1

  ALLOCATE(he(delta_range2*rangeE1_min+rangeE2_min:delta_range2*rangeE1_max+rangeE2_max),&
           & lng(delta_range2*rangeE1_min+rangeE2_min:delta_range2*rangeE1_max+rangeE2_max))
           he = 0; lng = 1.0D0; lnf = 1.0D0 !Initialization of the Variables

  nEbins = NINT((E1_max-E1_min)*dE1_inv*overlap_ratio*(E2_max-E2_min)*dE2_inv*overlap_ratio)
  nEbins = nEbins/number_windows
  mcs = mcs_bin_term*nEbins
  mcs_inv = 1.0D0/mcs
  CALL MPI_Barrier(MPI_COMM_WORLD,ierr)

  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  !Dividing the range of energy into 'nWindows' windos of the same size:
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!

  ALLOCATE (E1_min_aux(0:numprocs-1),E1_max_aux(0:numprocs-1),E2_min_aux(0:numprocs-1),&
            & E2_max_aux(0:numprocs-1),jWindow(0:number_windows))
  jWindow = -1
  DeltaEp2 = NINT((E2_max - E2_min)/(DBLE(mWindows) - DBLE(mWindows-1)*overlap_ratio))
  aux1 = E2_max
  aux2 = aux1 - DeltaEp2
  jWindow = 0
  DO ii=0,nWalkers-1
    E2_max_aux(ii) = aux1
    E2_min_aux(ii) = aux2
  END DO
  k = 1
  DO i=1,mWindows-1
    aux1 = aux2 + DeltaEp2*overlap_ratio
    aux2 = aux1 - DeltaEp2
    iWindows = 1 + 2*i
    DO j=0,iWindows-1
      jWindow(k) = i
      DO ii=0,nWalkers-1
        kk = k*nWalkers + ii
        E2_max_aux(kk) = aux1
        E2_min_aux(kk) = aux2
      END DO
       k = k + 1
    END DO
  END DO

  DeltaEp1 = NINT((E1_max - E1_min)/(DBLE(nWindows) - DBLE(nWindows-1)*overlap_ratio))
  aux4 = E1_min
  kk = numprocs-1
  DO i=mWindows-1,0,-1
    iWindows = 1 + 2*i
    aux1 = aux4
    aux2 = aux1 + DeltaEp1
    DO j=0,iWindows-1
      DO ii=0,nWalkers-1
        E1_min_aux(kk) = aux1
        E1_max_aux(kk) = aux2
        kk = kk - 1
      END DO
      aux1 = aux2 - DeltaEp1*overlap_ratio
      aux2 = aux1 + DeltaEp1
    END DO
    aux4 = aux4 + DeltaEp1*(1.0D0-overlap_ratio)
  END DO


ALLOCATE (my_neighbor_id(0:max_neighbors*numprocs*nWalkers-1), number_neighbor_windows(0:number_windows-1))
ALLOCATE (neighbor_windows(0:max_neighbors*number_windows-1),&
          & my_neighbor(0:numprocs*numprocs-1), position_aux(0:numprocs-1))

  number_neighbor_windows = 0
  k = 0
  kk = 0
  DO i=0,mWindows-1
    kk = kk + 2
    DO j=0,2*i
      IF(jWindow(k) == jWindow(k+1)) THEN
        iaux = k+1
        neighbor_windows(k*max_neighbors+number_neighbor_windows(k)) = iaux
        number_neighbor_windows(k) = number_neighbor_windows(k) + 1
        neighbor_windows(iaux*max_neighbors+number_neighbor_windows(iaux)) = k
        number_neighbor_windows(iaux) = number_neighbor_windows(iaux) + 1
      ELSE
        iaux = -1
      END IF
      IF(i < mWindows-1) THEN
        jaux = k + kk
        neighbor_windows(k*max_neighbors+number_neighbor_windows(k)) = jaux
        number_neighbor_windows(k) = number_neighbor_windows(k) + 1
        neighbor_windows(jaux*max_neighbors+number_neighbor_windows(jaux)) = k
        number_neighbor_windows(jaux) = number_neighbor_windows(jaux) + 1
      ELSE
        jaux = -1
      END IF
       k = k + 1
    END DO
  END DO

  my_neighbor = .FALSE.
  DO i=0,numprocs-1
    k=i/nWalkers
    DO j=0,number_neighbor_windows(k)-1
      DO kk=0,nWalkers-1
        my_neighbor(i*numprocs + neighbor_windows(4*k+j)*nWalkers + kk) = .TRUE.
      END DO
    END DO
  END DO

  my_neighbor_id = -1
  DO i=0,numprocs-1
    k = 0
    DO j=0,numprocs-1
      IF(my_neighbor(i*numprocs+j)) THEN
        kk = MOD(j,nWalkers)
        my_neighbor_id(i*max_neighbors*nWalkers+k*nWalkers+kk) = j
        IF (MOD(kk+1,nWalkers)==0) k = k + 1
      END IF
      IF (MOD(k,number_neighbor_windows(i/nWalkers))==0) k = 0
    END DO
  END DO

!IF (myid==0) THEN
!  DO i=0,numprocs-1
!    kk=i/nWalkers
!    DO j=0,number_neighbor_windows(kk)-1
!      DO k=0,nWalkers-1
!        PRINT*,i,my_neighbor_id(i*4*nWalkers+j*nWalkers + k)
!      END DO
!    END DO
!  END DO
!END IF

  CALL MPI_Barrier(MPI_COMM_WORLD,ierr)

  rangeE1_min_aux = NINT(E1_min_aux(myid)*dE1_inv)
  rangeE1_max_aux = NINT(E1_max_aux(myid)*dE1_inv)
  rangeE2_min_aux = NINT(E2_min_aux(myid)*dE2_inv)
  rangeE2_max_aux = NINT(E2_max_aux(myid)*dE2_inv)

  DO i=rangeE1_min,rangeE1_min_aux
    DO j=rangeE2_min,rangeE2_max
      k = i*delta_range2 + j
      lng(k) = lng(k) + 1.0D6 + 1.0D2*(i-rangeE1_min_aux)*(i-rangeE1_min_aux)
    END DO
  END DO

  DO i=rangeE1_max_aux,rangeE1_max
    DO j=rangeE2_min,rangeE2_max
      k = i*delta_range2 + j
      lng(k) = lng(k) + 1.0D6 + 1.0D2*(i-rangeE1_max_aux)*(i-rangeE1_max_aux)
    END DO
  END DO

  DO i=rangeE1_min,rangeE1_max
    DO j=rangeE2_min,rangeE2_min_aux
      k = i*delta_range2 + j
      lng(k) = lng(k) + 1.0D6 + 1.0D2*(j-rangeE2_min_aux)*(j-rangeE2_min_aux)
    END DO
  END DO

  DO i=rangeE1_min,rangeE1_max
    DO j=rangeE2_max_aux,rangeE2_max
      k = i*delta_range2 + j
      lng(k) = lng(k) + 1.0D6 + 1.0D2*(j-rangeE2_max_aux)*(j-rangeE2_max_aux)
    END DO
  END DO

  DO
    found_core = .FALSE.
    CALL sweep()
     IF((E_current(1) > rangeE1_min_aux) .AND. (E_current(1) < rangeE1_max_aux) &
       & .AND. (E_current(2) > rangeE2_min_aux) .AND. (E_current(2) < rangeE2_max_aux)) found_core = .TRUE.
     PRINT*,myid,E_current(1),rangeE1_min_aux,rangeE1_max_aux,E_current(2),rangeE2_min_aux,rangeE2_max_aux, found_core
     IF (myid == 0) THEN !Processor 0 will process this check-up.
       found = found_core !The global flatness status recieves the flatness status of the walker 0.
       DO i=1,numprocs-1
         !recieve all core flatness status results
         CALL MPI_Recv(found_core,1,MPI_LOGICAL,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,mpistatus,ierr)
         IF (.NOT.found_core) found = .FALSE. !If one walker is not flat the global flatness status would be false
       END DO
     ELSE
       !Send the result of the core flatness status to processor 0
       CALL MPI_Send(found_core,1,MPI_LOGICAL,0,myid,MPI_COMM_WORLD,ierr)
     END IF
     CALL MPI_Bcast(found,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
     IF(found) EXIT
  END DO

17  CALL energy(current,E_current)

  CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
  IF(myid==0) PRINT*,'The configuration inside the energy window was found.'
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  !Redefining the variables and reallocating the memory:
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  E1_min = E1_min_aux(myid)
  E1_max = E1_max_aux(myid)
  E2_min = E2_min_aux(myid)
  E2_max = E2_max_aux(myid)

  rangeE1_min = NINT(E1_min*dE1_inv)
  rangeE1_max = NINT(E1_max*dE1_inv)
  rangeE2_min = NINT(E2_min*dE2_inv)
  rangeE2_max = NINT(E2_max*dE2_inv)
  delta_range1 = rangeE1_max-rangeE1_min + 1
  delta_range2 = rangeE2_max-rangeE2_min + 1

  PRINT*,'Processor ' ,myid, 'E1_min', E1_min,'E1_max', E1_max,&
  &'E2_min', E2_min,'E2_max', E2_max
  CALL MPI_Barrier(MPI_COMM_WORLD,ierr)

  mcs = mcs_bin_term*nEbins
  mcs_inv = 1.0D0/mcs
  IF(myid==0)  PRINT*,'$$$$$$$$$$$$$$$$ Monte Carlo Sweeps = ', mcs, '$$$$$$$$$$$$$$$$'

  DEALLOCATE (he,lng)
  ALLOCATE(he(delta_range2*rangeE1_min+rangeE2_min:delta_range2*rangeE1_max+rangeE2_max),&
          & lng(delta_range2*rangeE1_min+rangeE2_min:delta_range2*rangeE1_max+rangeE2_max))
  he = 0; lng = 1.0D0; lnf = 1.0D0 !Initialization of the Variables

  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  !Finishing the initialization subroutine:
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
  IF(myid ==0) PRINT*,'The initialization was done.'

  RETURN
END SUBROUTINE initial

!------------------------------------------------------------------------------!
!Performing mcs Monte Carlo sweeps:
!------------------------------------------------------------------------------!
SUBROUTINE sweep()
  USE integer_kind
  USE WL_variables
  USE spins_variables
  IMPLICIT NONE
  INTEGER(K4B) :: j, imcs, icount_replica!, i, mcs_aux, n_numbers
  INTEGER(K15) :: k
!  mcs_aux = mcs/n_mcs
!  n_numbers = mcs_aux*3*L2
!  ALLOCATE(random_numbers(1:n_numbers))
  accept_counter = 0
  exchange_counter = 0
  imcs = 0
  icount_replica = 0
  DO k=1,mcs
    CALL energy(current,E_current)
    trial = current
    E_trial = E_current
    numbers_count = 0
!    CALL vranmar(n_numbers,random_numbers)
    imcs = imcs + 1
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    !Single spin rotation
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    DO j = 1,L2
        cont = counter_s(j)
        delta_sx = -2*current%sx(cont)
        delta_tx = 0
        CALL delta_energy(cont)
        CALL wangLandau_step()
        IF(accept) THEN
          current%U1 = trial%U1
          current%U2 = trial%U2
          E_current = E_trial
          current%Sx(cont) = -current%Sx(cont)
        END IF
        cont = counter_t(j)
        delta_sx = 0
        delta_tx = -2*current%tx(cont)
        CALL delta_energy(cont)
        CALL wangLandau_step()
        IF(accept) THEN
          current%U1 = trial%U1
          current%U2 = trial%U2
          current%Tx(cont) = -current%Tx(cont)
          E_current = E_trial
        END IF
    END DO
    CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    !REPLICA EXCHANGE
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    trial = current
    E_trial = E_current
    IF (MOD(imcs,ireplica) == 0) THEN
      icount_replica = icount_replica + 1
      CALL replica_exchange()
    END IF
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  END DO
  accept_ratio = 0.5D0*DBLE(accept_counter)*L2_inv*mcs_inv
  replica_ratio = DBLE(exchange_counter)/DBLE(icount_replica)

  RETURN
END SUBROUTINE sweep

!------------------------------------------------------------------------------!
!Calcula a variacao da energia e da magnetizacao
!------------------------------------------------------------------------------!
SUBROUTINE delta_energy(i)
USE WL_variables
USE spins_variables
USE integer_kind
IMPLICIT NONE
INTEGER(K4B), INTENT(IN) :: i
INTEGER(K15) :: ivx1, ivx2, ivx3
INTEGER(K4B) :: idir,iesq,jdir,jesq

   idir = vizXdir(i)
   iesq = vizXesq(i)
   jdir = vizYdir(i)
   jesq = vizYesq(i)
  ! Soma das componentes dos primeiros vizinhos do sitio 'i'
  ivx1 = current%sx(idir) + current%sx(iesq) + &
      & current%sx(jdir) + current%sx(jesq)
  ivx2 = current%tx(idir) + current%tx(iesq) + &
      & current%tx(jdir) + current%tx(jesq)
  ivx3 = current%sx(idir)*current%tx(idir) + &
      & current%sx(iesq)*current%tx(iesq) + &
      & current%sx(jdir)*current%tx(jdir) + &
      & current%sx(jesq)*current%tx(jesq)
  delta_U1 = -DBLE(delta_sx*ivx1 + delta_tx*ivx2)
  delta_U2 = -DBLE(delta_sx*current%tx(i)*ivx3 + delta_tx*current%sx(i)*ivx3)
  trial%U1 = delta_U1 + current%U1
  trial%U2 = delta_U2 + current%U2
  E_trial(1) = NINT(trial%U1*dE1_inv)
  E_trial(2) = NINT(trial%U2*dE2_inv)
  E_trial(3) = delta_range2*E_trial(1) + E_trial(2)

RETURN
END SUBROUTINE delta_energy

!------------------------------------------------------------------------------!
!Calcula a energia e a magnetizacao do sistema
!------------------------------------------------------------------------------!
SUBROUTINE energy(config_aux,E_config_aux)
USE integer_kind
USE spins_variables
USE WL_variables
IMPLICIT NONE
INTEGER(K4B)        :: k, idir, jdir, kk, counter_aux, ivx1, ivx2, ivx3
INTEGER(K15), DIMENSION(1:3), INTENT(OUT) :: E_config_aux
TYPE (spins_config), INTENT(INOUT) :: config_aux
DOUBLE PRECISION:: rand

  !Inicializando as variaveis
  config_aux%U1 = 0
  config_aux%U2 = 0
  DO k=1,L2
    idir = vizXdir(k)
    jdir = vizYdir(k)
    ! soma das componentes dos primeiros vizinhos
    ivx1 = config_aux%sx(idir) + config_aux%sx(jdir)
    ivx2 = config_aux%tx(idir) + config_aux%tx(jdir)
    ivx3 = config_aux%sx(idir)*config_aux%tx(idir) + &
        & config_aux%sx(jdir)*config_aux%tx(jdir)
    !Contribuicao aa energia do sitio k
    config_aux%U1 = config_aux%U1 - DBLE(config_aux%sx(k)*ivx1 + config_aux%tx(k)*ivx2)
    config_aux%U2 = config_aux%U2 - DBLE(config_aux%sx(k)*config_aux%tx(k)*ivx3)
    !Embaralhando ordem dos passos de Monte Carlo
    ! Rede S
    rand = ranmar()
    kk = NINT(L2_dble*rand+1.0D0)
    counter_aux = counter_s(k)
    counter_s(k) = counter_s(kk)
    counter_s(kk) = counter_aux
    ! Rede T
    rand = ranmar()
    kk = NINT(L2_dble*rand+1.0D0)
    counter_aux = counter_t(k)
    counter_t(k) = counter_t(kk)
    counter_t(kk) = counter_aux
  END DO
  E_config_aux(1) = NINT(config_aux%U1*dE1_inv)
  E_config_aux(2) = NINT(config_aux%U2*dE2_inv)
  E_config_aux(3) = delta_range2*E_config_aux(1) + E_config_aux(2)

  RETURN
END SUBROUTINE energy

!------------------------------------------------------------------------------!
!Wang-Landau Algorithm
!------------------------------------------------------------------------------!
SUBROUTINE wangLandau_step()
  USE WL_variables
  USE spins_variables
  IMPLICIT NONE
  REAL(8) :: aux1,aux2

  accept = .FALSE.
  !- - - - - - - - - - - - - - - - - - - - - - - - - - !
  !Checking if the new state is in the range of energy:
  !- - - - - - - - - - - - - - - - - - - - - - - - - - !
  IF (E_trial(1) > rangeE1_max) THEN
!    Print*, 'E_trial > rangeE_max'
    GO TO 15
  ELSE
    IF (E_trial(1) < rangeE1_min) THEN
!      Print*, 'E_trial < rangeE_min'
      GO TO 15
    ELSE
      IF (E_trial(2) > rangeE2_max) THEN
!    Print*, 'E_trial > rangeE_max'
        GO TO 15
      ELSE
        IF (E_trial(2) < rangeE2_min) THEN
!      Print*, 'E_trial < rangeE_min'
          GO TO 15
        END IF
      END IF
    END IF
  END IF
  !- - - - - - - - - - - - - - - - - - - - - - - - -!
  !Wang Landau sweep
  !- - - - - - - - - - - - - - - - - - - - - - - - -!
  IF (lng(E_trial(3)) < lng(E_current(3))) THEN
    E_current = E_trial
    accept_counter = accept_counter + 1
    accept = .TRUE.
!    Print*, 'accepted lng(E_trial) < lng(E_current)'
  ELSE
!    numbers_count = numbers_count + 1
!    aux1 = random_numbers(numbers_count)
    aux1 = ranmar()
    aux2 = DEXP(lng(E_current(3)) - lng(E_trial(3)))
    IF (aux1 < aux2) THEN
      E_current = E_trial
      accept_counter = accept_counter + 1
      accept = .TRUE.
!      Print*, 'accepted g(E_trial > g(E_current)'
!    ELSE
!      Print*, 'REJECTED'
    END IF
  END IF
  !- - - - - - - - - - - - - - - - - - - - - - - - -!
  !Updating the histogram and the density of states:
  !- - - - - - - - - - - - - - - - - - - - - - - - -!
15 lng(E_current(3)) = lng(E_current(3)) + lnf
   he(E_current(3)) = he(E_current(3)) + 1

  RETURN
END SUBROUTINE wangLandau_step
!------------------------------------------------------------------------------!
!Replica exchange subroutine
!------------------------------------------------------------------------------!
SUBROUTINE replica_exchange()
  USE integer_kind
  USE WL_variables
  USE spins_variables
  USE MPI_new_variables
  IMPLICIT NONE
  INTEGER(K4B) :: i, j, k, jaux, iaux
  INTEGER(K4B),ALLOCATABLE,DIMENSION(:)  :: partner
  REAL(8) :: aux1,aux2,aux3,aux4
  REAL(8),DIMENSION(1:2) :: ge_aux
  LOGICAL :: accept_exchange, range_bin, range_bin_partner
  LOGICAL :: main_proc, secondary_proc
  LOGICAL,ALLOCATABLE,DIMENSION(:) :: picked

  ALLOCATE(partner(0:numprocs-1))
  partner = -1
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  !Defining the walkers partners for the replica exchange:
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  IF (myid == 0) THEN
    ALLOCATE(picked(0:numprocs-1))
      !After that, the values inside a window will be shuffled.
    picked = .TRUE.
    DO i=0,numprocs-1
      iaux = i/nWalkers
      j = NINT(DBLE(number_neighbor_windows(iaux)-1)*ranmar())
      k = NINT(DBLE(nWalkers-1)*ranmar())
      jaux = my_neighbor_id(i*4*nWalkers + j*nWalkers + k)
      IF (picked(i) .AND. picked(jaux)) THEN
        partner(i) = jaux
        partner(jaux) = i
        picked(i) = .FALSE.
        picked(jaux) = .FALSE.
      END IF
    END DO
    DEALLOCATE(picked)
  END IF
  CALL MPI_barrier(MPI_COMM_WORLD,ierr)
  CALL MPI_Bcast(partner,numprocs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  main_proc = .FALSE.
  secondary_proc = .FALSE.
  IF(partner(myid) .GE. 0) THEN
    IF(myid < partner(myid)) THEN
      main_proc = .TRUE.
    ELSE
      secondary_proc = .TRUE.
    END IF
!  ELSE
!    GO TO 17
  END IF

!  PRINT*,myid,partner(myid),main_proc,secondary_proc
 CALL MPI_barrier(MPI_COMM_WORLD,ierr)
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  !REPLICA EXCANGE
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
  range_bin = .FALSE.
  IF (main_proc) THEN
    !Sending/Recieving the energy bin flag to/from the secondary processor.
    CALL MPI_Sendrecv(E_current,3,MPI_INTEGER,partner(myid),myid,&
          & E_trial,3,MPI_INTEGER,partner(myid),partner(myid)+6*numprocs,MPI_COMM_WORLD,mpistatus,ierr)
    !Checking if the energy of the state of the secondary processor is in the range of energy of the main processor.
    IF ((E_trial(1) > rangeE1_min) .AND. (E_trial(1) < rangeE1_max) &
       & .AND. (E_trial(2) > rangeE2_min) .AND. (E_trial(2) < rangeE2_max)) range_bin = .TRUE.
    !Sending/Recieving the result of the previous check to/from the secondary processor.
    CALL MPI_Sendrecv(range_bin,1,MPI_LOGICAL,partner(myid),myid+7*numprocs,&
         & range_bin_partner,1,MPI_LOGICAL,partner(myid),partner(myid)+numprocs,MPI_COMM_WORLD,mpistatus,ierr)
    !If the walkers are in the same range we can calculate the acceptance ratio (P_acc = min[1,gi(X)gj(Y)/gi(Y)gj(X)]
    IF (range_bin .AND. range_bin_partner) THEN
      E_trial(3) = delta_range2*E_trial(1) + E_trial(2)
      !Reciving the density of states from the secondary processor {ge_aux(1) = ln[gj(Y)] and ge_aux(2) = ln[gj(X))]}
      CALL MPI_Recv(ge_aux,2,MPI_DOUBLE_PRECISION,partner(myid),myid+2*numprocs,MPI_COMM_WORLD,mpistatus,ierr)
      aux1 = lng(E_current(3)) + ge_aux(1)
      aux2 = lng(E_trial(3)) + ge_aux(2)
      !Wang-Landau Algorithm
      IF (aux1 < aux2) THEN ! P_acc=min[gi(X)gj(Y)/gi(Y)gj(X)] [It will be accept if gi(X)gj(Y) < gi(Y)gj(X)]
         accept_exchange = .TRUE. !The exchange is accepted - acceptance status = True
      ELSE
!         numbers_count = numbers_count + 1
        aux3 = ranmar()!random_numbers(numbers_count)
        aux4 = DEXP(aux1 - aux2)
        IF (aux3 < aux4) THEN ! P_acc=min[gi(X)gj(Y)/gi(Y)gj(X)]  [It will be accept if the random number < gi(X)gj(Y)/gi(Y)gj(X)]
          accept_exchange = .TRUE. !The exchange is accepted - acceptance status = True
        ELSE
          accept_exchange = .FALSE. !The exchange is not accepted - acceptance status = False
        END IF
      END IF
      !Sending the acceptance status to the secondary processor.
      CALL MPI_Send(accept_exchange,1,MPI_LOGICAL,partner(myid),myid+3*numprocs,MPI_COMM_WORLD,ierr)
      !If it is accepted the subroutine will exchange the configurations of the main and secondary processor.
      IF(accept_exchange) THEN
        exchange_counter = exchange_counter + 1 !Count how many times the exchange was accepted.
        E_current = E_trial
        !Exchanging the configuration.
        CALL MPI_Sendrecv_replace(current,1,MPI_config_spins,partner(myid),myid+5*numprocs,&
              & partner(myid),partner(myid)+4*numprocs,MPI_COMM_WORLD,mpistatus,ierr)
      END IF
    END IF
  ELSE
    IF(secondary_proc) THEN
      !Sending/Recieving the energy bin flag to/from the main processor.
      CALL MPI_Sendrecv(E_current,3,MPI_INTEGER,partner(myid),myid+6*numprocs,&
           & E_trial,3,MPI_INTEGER,partner(myid),partner(myid),MPI_COMM_WORLD,mpistatus,ierr)
      !Checking if the energy of the state of the main processor is in the range of energy of the secondary processor.
      IF ((E_trial(1) > rangeE1_min) .AND. (E_trial(1) < rangeE1_max) &
         & .AND. (E_trial(2) > rangeE2_min) .AND. (E_trial(2) < rangeE2_max)) range_bin = .TRUE.
      !Sending/Recieving the result of the previous check to/from the main processor.
      CALL MPI_Sendrecv(range_bin,1,MPI_LOGICAL,partner(myid),myid+numprocs,&
           & range_bin_partner,1,MPI_LOGICAL,partner(myid),partner(myid)+7*numprocs,MPI_COMM_WORLD,mpistatus,ierr)
      IF (range_bin .AND. range_bin_partner) THEN
        E_trial(3) = delta_range2*E_trial(1) + E_trial(2)
        !Sending the density of states to the main processor {ge_aux(1) = ln[gj(Y)] and ge_aux(2) = ln[gj(X))]}
        ge_aux(1) = lng(E_current(3))
        ge_aux(2) = lng(E_trial(3))
        CALL MPI_Send(ge_aux,2,MPI_DOUBLE_PRECISION,partner(myid),partner(myid)+2*numprocs,MPI_COMM_WORLD,ierr)
        !Reciving the acceptance status from the main processor.
        CALL MPI_Recv(accept_exchange,1,MPI_LOGICAL,partner(myid),partner(myid)+3*numprocs,MPI_COMM_WORLD,mpistatus,ierr)
        !If it is accepted the subroutine will exchange the configurations of the main and secondary processor.
        IF(accept_exchange) THEN
          exchange_counter = exchange_counter + 1 !Count how many times the exchange was accepted.
          E_current = E_trial
          !Exchanging the configuration.
          CALL MPI_Sendrecv_replace(current,1,MPI_config_spins,partner(myid),myid+4*numprocs,&
               & partner(myid),partner(myid)+5*numprocs,MPI_COMM_WORLD,mpistatus,ierr)
        END IF
      END IF
    ELSE
    END IF
  END IF

  lng(E_current(3)) = lng(E_current(3)) + lnf !Density of states
  he(E_current(3)) = he(E_current(3)) + 1 !Histogram

  DEALLOCATE(partner)

  RETURN
END SUBROUTINE replica_exchange

!------------------------------------------------------------------------------!
!Verifica se o histograma esta flat
!------------------------------------------------------------------------------!
SUBROUTINE flathist()
  USE integer_kind
  USE WL_variables
  USE MPI_new_variables
  IMPLICIT NONE
  INTEGER(K15) :: i
  INTEGER(K4B) :: kk
  REAL(8) :: lngE_med(delta_range2*rangeE1_min+rangeE2_min:delta_range2*rangeE1_max+rangeE2_max),&
           & media,mini,aux,p2
  LOGICAL :: flat_core

  !--------------------------------------------------!
  !check if the walker myid is flat
  !--------------------------------------------------!
  flat_core = .FALSE.
  kk = 0
  media = 0.0D0
  mini = 1.0D9
  DO i = delta_range2*rangeE1_min+rangeE2_min, delta_range2*rangeE1_max+rangeE2_max
    IF (lng(i) > 1.0D2) THEN
      kk = kk + 1
      aux = DBLE(he(i))
      media = media + aux
      IF (aux .LT. mini) mini = aux
    END IF
  END DO
  media=media/DBLE(kk)
  p2=mini/media
  IF (p2 > p) flat_core = .TRUE.
  PRINT*,'id', myid, 'flatness', p2, flat_core,'acceptance', accept_ratio, replica_ratio
  !--------------------------------------------------!
  !Checking if all walkers are flat.
  !--------------------------------------------------!
  CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
  IF (myid == 0) THEN !Processor 0 will process this check-up.
    flat = flat_core !The global flatness status recieves the flatness status of the walker 0.
    DO i=1,numprocs-1
      !recieve all core flatness status results
      CALL MPI_Recv(flat_core,1,MPI_LOGICAL,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,mpistatus,ierr) 
      IF (.NOT.flat_core) flat = .FALSE. !If one walker is not flat the global flatness status would be false
    END DO
  ELSE
    !Send the result of the core flatness status to processor 0
    CALL MPI_Send(flat_core,1,MPI_LOGICAL,0,myid,MPI_COMM_WORLD,ierr)
  END IF
  !Send the result of flatness of all walkers to all processors
  CALL MPI_Bcast(flat,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
  !------------------------------------------------------------------------------------------------------------------------!
  !If all walkers are flat we have to distribute the average of the density of states of/to all walkers inside a subwindow
  !------------------------------------------------------------------------------------------------------------------------!
  IF (flat) THEN
    IF(myid==0) PRINT*,'############## FLAT HISTOGRAM #################',lnf
    lnf = 0.5D0*lnf !Reduce the modification fator f
    he = 0          !Reset the histogram to zero
    flat = .FALSE.  !Reset flatness variable status
    IF (lnf > lnf_min) THEN
      IF (energy_id == 0) THEN !sub_processor 0 will calculate de average
        lngE_med = lng
        DO i=1,nWalkers-1
          !recieve all density of states of the sub-window
          CALL MPI_Recv(lng,delta_range1*delta_range2,&
                       &MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG,energy_comm,mpistatus,ierr)
          lngE_med = lngE_med + lng
        END DO
        lngE_med = lngE_med/DBLE(nWalkers)
        lng = lngE_med
      ELSE
        !send density of the state to sub processor 0
        CALL MPI_Send(lng,delta_range1*delta_range2,&
                     &MPI_DOUBLE_PRECISION,0,energy_id,energy_comm,ierr)
      END IF
      !distribute the density of states to all walkers inside the subwindow
      CALL MPI_Bcast(lng,delta_range1*delta_range2,&
                    &MPI_DOUBLE_PRECISION,0,energy_comm,ierr)
    END IF
  ELSE
    IF(myid==0) PRINT*,'****NOT FLAT HISTOGRAM *****', lnf
  END IF

  RETURN
END SUBROUTINE flathist
!------------------------------------------------------------------------------!
!Pseudorandom numbers generator
!------------------------------------------------------------------------------!
SUBROUTINE rmarin(ij, kl)
!  This subroutine and the next function generate random numbers. See
!  the comments for SA for more information. The onL changes from the
!  orginal code is that (1) the test to make sure that RMARIN runs first
!  was taken out since SA assures that this is done (this test didn't
!  compile under IBM's VS Fortran) and (2) typing ivec as integer was
!  taken out since ivec isn't used. With these exceptions, all following
!  lines are original.

! This is the initialization routine for the random number generator
!     RANMAR()
! NOTE: The seed variables can have values between:    0 <= IJ <= 31328
!                                                      0 <= KL <= 30081
USE integer_kind
USE raset1
IMPLICIT NONE
INTEGER(K4B), INTENT(IN) :: ij, kl
INTEGER(K4B) :: i, j, k, l, ii, jj, m
REAL(8) :: s, t

IF( ij < 0  .OR.  ij > 31328  .OR. kl < 0  .OR.  kl > 30081 ) THEN
  WRITE(*, '(A)') ' The first random number seed must have a value ',  &
               'between 0 AND 31328'
  WRITE(*, '(A)') ' The second seed must have a value between 0 and 30081'
  STOP
END IF

i = MOD(ij/177, 177) + 2
j = MOD(ij, 177) + 2
k = MOD(kl/169, 178) + 1
l = MOD(kl, 169)
DO ii = 1, 97
  s = 0.0D0
  t = 0.5D0
  DO jj = 1, 24
    m = MOD(MOD(i*j, 179)*k, 179)
    i = j
    j = k
    k = m
    l = MOD(53*l + 1, 169)
    IF (MOD(l*m, 64) >= 32) THEN
      s = s + t
    END IF
    t = 0.5D0*t
  END DO
  u(ii) = s
END DO
cc = 362436.0D0/16777216.0D0
cd = 7654321.0D0/16777216.0D0
cm = 16777213.0D0/16777216.0D0
i97 = 97
j97 = 33

RETURN
END SUBROUTINE rmarin

!------------------------------------------------------------------------------!
! This is the random number generator proposed by George Marsaglia
! in Florida State University Report: FSU-SCRI-87-50
!------------------------------------------------------------------------------!
FUNCTION ranmar() RESULT(fn_val)
USE raset1
IMPLICIT NONE
REAL(8) :: fn_val
! Local variable
REAL(8):: uni

  uni = u(i97) - u(j97)
  IF( uni < 0.0D0 ) uni = uni + 1.0D0
  u(i97) = uni
  i97 = i97 - 1
  IF(i97 == 0) i97 = 97
  j97 = j97 - 1
  IF(j97 == 0) j97 = 97
  cc = cc - cd
  IF( cc < 0.0D0 ) cc = cc + cm
  uni = uni - cc
  IF( uni < 0.0D0 ) uni = uni + 1.0D0
!  IF( uni == 0.0D0 ) uni = 2.0D-38
  fn_val = uni

  RETURN
END FUNCTION ranmar

!------------------------------------------------------------------------------!
!The same generator as above but for a vector of n random numbers 
!------------------------------------------------------------------------------!
!SUBROUTINE vranmar(n, fn_val)
!USE integer_kind
!USE raset1
!IMPLICIT NONE
!  INTEGER(K4B), INTENT(IN) :: n
!  REAL(8), DIMENSION(1:n),INTENT(OUT) :: fn_val
!
!  REAL(8) :: uni
!  INTEGER(K4B) :: i
!
!  DO i = 1, n
!    uni = u(i97) - u(j97)
!    IF( uni < 0.0D0 ) uni = uni + 1.0D0
!    u(i97) = uni
!    i97 = i97 - 1
!    IF(i97 == 0) i97 = 97
!    j97 = j97 - 1
!    IF(j97 == 0) j97 = 97
!    cc = cc - cd
!    IF( cc < 0.0D0 ) cc = cc + cm
!    uni = uni - cc
!    IF( uni < 0.0D0 ) uni = uni + 1.0D0
!!    IF( uni == 0.0D0 ) uni = 2.0D-38
!    fn_val(i) = uni
!  END DO
!
!  RETURN
!END SUBROUTINE vranmar

!------------------------------------------------------------------------------!
!Initial random number generator (numerical recipes)
!------------------------------------------------------------------------------!
DOUBLE PRECISION FUNCTION ran_init(idum)
  USE integer_kind
  IMPLICIT NONE
  INTEGER(K4B), INTENT(INOUT)  :: idum
  INTEGER(K4B), PARAMETER      :: IA=16807,IM=2147483647,IQ=127773,IR=2836
  REAL(8), SAVE                :: am
  INTEGER(K4B), SAVE           :: ix=-1,iy=-1,k

  IF (idum <= 0 .OR. iy < 0) THEN
    am = nearest(1.0,-1.0)/DBLE(IM)
    iy = ior(ieor(888889999,abs(idum)),1)
    ix = ieor(777755555,abs(idum))
    idum = abs(idum) + 1
  END IF
  ix = ieor(ix,ishft(ix,13))
  ix = ieor(ix,ishft(ix,-17))
  ix = ieor(ix,ishft(ix,5))
  k = iy/IQ
  iy = IA*(iy - k*IQ) - IR*k
  IF (iy < 0) iy = iy + IM
  ran_init = am*ior(iand(IM,ieor(ix,iy)),1)

  RETURN
END FUNCTION ran_init
!------------------------------------------------------------------------------!

END PROGRAM WangLandau_xy_model
!##################################################################################################!
