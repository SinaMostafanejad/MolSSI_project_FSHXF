PROGRAM calculate_BO
  !
  USE kinds
  USE m_potential,  ONLY: define_pot_type, extpotential, potentials, extgradient
  USE m_solver,     ONLY: grid_diag_spjd
  USE m_math,       ONLY: integrate
  !
  IMPLICIT NONE
  !
  REAL(dp) :: pos_inp
  CHARACTER(LEN=256) :: model_inp
  INTEGER :: nstates
  !
  INTEGER, PARAMETER :: nx= 401
  INTEGER, PARAMETER :: lr= 200
  REAL(dp), PARAMETER :: xmin = -20.0_dp
  REAL(dp), PARAMETER :: xmax = 20.0_dp
  !
  REAL(dp) :: xe(nx), vext(nx), grad(1,nx),las_der(1,nx), laserfield
  REAL(dp), ALLOCATABLE :: ev(:), nac(:,:,:), evec(:,:), force(:,:)
  REAL(dp), ALLOCATABLE :: evec_l(:,:), evec_r(:,:), norm_lr(:,:)
  REAL(dp), ALLOCATABLE :: laser_coup(:,:,:)
  REAL(dp) :: dx
  !
  INTEGER :: i, noe, ist, jst
  REAL(dp) :: integral, denom
  REAL(dp) :: integral_las
  !
  REAL(dp) :: xn(nx), pot(nx,2), nacs(nx)
  REAL(dp) :: force_num(nx)
  !
  REAL(dp) :: norm
  !
  REAL(dp) :: time
  REAL(dp) :: step
  !
  logical :: exist
  !
  WRITE(6,*) 'Model Hamiltonian'
  READ(5,*) model_inp
  WRITE(6,*) 'Nuclear position?'
  READ(5,*) pos_inp
  WRITE(6,*) '# of BO states?'
  READ(5,*) nstates
  WRITE(6,*) 'Time'
  READ(5,*) time
  WRITE(6,*) 'Step size'
  READ(5,*) step
  !
  step=step*0.02418884344 !step size in atomic units
  !
  ALLOCATE(ev(nstates),evec(nx,nstates),evec_l(lr,nstates),evec_r(lr,nstates),nac(1,nstates,nstates),force(1,nstates),laser_coup(1,nstates,nstates))
  ALLOCATE(norm_lr(2,2))
  ! ... Calculate BO stuff
  !
  ! ... Set up grid
  dx = (xmax-xmin)/DBLE(nx-1)
  DO i = 1, nx
    xe(i) = xmin + (i-1)*dx
  END DO
  ! ... Set up potential
  CALL define_pot_type(potentials)
  ! ... external potentials
  CALL extpotential(TRIM(model_inp),pos_inp,nx,xe,time,step,vext,laserfield)
  ! ... write laser field
  open(33,file='laser.data',status='unknown',access='append')
  write(33,*) laserfield, time
  close(33)
  ! ... solve eigenvalue equation
  CALL grid_diag_spjd(nx,dx,vext,ev,evec,noe,.TRUE.)
  ! ... NAC calculation
  ! ... potential gradient
  CALL extgradient(TRIM(model_inp),pos_inp,nx,xe,grad,las_der,time,step)
  norm = integrate(nx,evec(:,1)*evec(:,1))*dx
  evec(:,1) = evec(:,1)/SQRT(norm)
  norm = integrate(nx,evec(:,1)*evec(:,1))*dx
  WRITE(6,*) 'NORM 1: ', norm
  norm = integrate(nx,evec(:,2)*evec(:,2))*dx
  evec(:,2) = evec(:,2)/SQRT(norm)
  norm = integrate(nx,evec(:,2)*evec(:,2))*dx
  !
  !Localization probabilities
  do i=1,lr !loading electronic wf defined onnegative and positive z values
    evec_l(i,:)=evec(i,:)
    evec_r(i,:)=evec(201+i,:)
  enddo
  norm_lr(1,1) = integrate(lr,evec_l(:,1)*evec_l(:,1))*dx 
  norm_lr(1,2) = integrate(lr,evec_l(:,2)*evec_l(:,2))*dx
  norm_lr(2,1) = integrate(lr,evec_r(:,1)*evec_r(:,1))*dx
  norm_lr(2,2) = integrate(lr,evec_r(:,2)*evec_r(:,2))*dx
  !
  WRITE(6,*) 'NORM 2: ', norm
  DO ist = 1, nstates
    force(1,ist) = -integrate(nx,evec(:,ist)*evec(:,ist)*grad(1,:))*dx
    DO jst = ist+1, nstates
!      integral = integrate(nx,CONJG(evec(:,ist))*evec(:,jst)*grad(1,:))*dx
      integral = integrate(nx,evec(:,ist)*evec(:,jst)*grad(1,:))*dx
      denom = ev(jst)-ev(ist)

      integral_las=integrate(nx,evec(:,ist)*evec(:,jst)*las_der(1,:))*dx

      IF(ABS(denom) .LT. 1.0d-12 ) THEN
        nac(1,ist,jst) = 1.0d4
        laser_coup(1,ist,jst) = 1.0d4
      ELSE 
        nac(1,ist,jst) = integral / denom
        laser_coup(1,ist,jst) = integral_las / denom
      END IF
    nac(1,jst,ist) = - nac(1,ist,jst)
    laser_coup(1,jst,ist) = -laser_coup(1,ist,jst)
    END DO
  END DO
  !
  OPEN(15,FILE="ENERGY.DAT") 
  OPEN(16,FILE="FORCE.DAT") 
  OPEN(17,FILE="NAC.DAT")
  OPEN(18,FILE="LASER_COUP.dat")
  !
  DO ist=1, nstates
    WRITE(15,*) ev(ist)
    WRITE(16,"(i1)") ist
!    WRITE(16,*) force(1,ist)
    WRITE(16,'(ES22.15)') force(1,ist)
    DO jst=1, nstates
      WRITE(17,"(i1,1x,i1)") ist, jst
      WRITE(17,*) nac(1,ist,jst)
      WRITE(18,"(i1,1x,i1)") ist, jst
      WRITE(18,*) laser_coup(1,ist,jst)
    END DO
  END DO
  !
  inquire(file="NORM_LR.dat", exist=exist)
  if (exist) then
    OPEN(19,FILE="NORM_LR.dat", status="old", position="append", action="write")
  else
    OPEN(19,FILE="NORM_LR.dat", status="new", action="write")
  endif
  WRITE(19,*) norm_lr(1,1), norm_lr(2,1), norm_lr(1,2), norm_lr(2,2)
  CLOSE(15)
  CLOSE(16)
  CLOSE(17)
  CLOSE(18)
  CLOSE(19)
!  OPEN(99,FILE="BO.DAT")
!  force_num = 0.0_dp
!  DO i = 1, nx
!    !
!    xn(i) = -7.0_dp + (i-1)*DBLE(14.0_dp/(nx-1))
!    ! ... Set up potential
!    CALL define_pot_type(potentials)
!    ! ... external potentials
!    CALL extpotential(TRIM(model_inp),xn(i),nx,xe,vext)
!    ! ... solve eigenvalue equation
!    CALL grid_diag_spjd(nx,dx,vext,pot(i,:),evec,noe,.TRUE.)
!    ! ... NAC calculation
!    ! ... potential gradient
!    CALL extgradient(TRIM(model_inp),xn(i),nx,xe,grad)
!    norm = integrate(nx,evec(:,1)*evec(:,1))*dx
!    evec(:,1) = evec(:,1)/SQRT(norm)
!    norm = integrate(nx,evec(:,1)*evec(:,1))*dx
!    WRITE(6,*) 'NORM 1: ', norm
!    norm = integrate(nx,evec(:,2)*evec(:,2))*dx
!    evec(:,2) = evec(:,2)/SQRT(norm)
!    norm = integrate(nx,evec(:,2)*evec(:,2))*dx
!    WRITE(6,*) 'NORM 2: ', norm
!    force(1,2) = - integrate(nx,evec(:,2)*evec(:,2)*grad(2,:))*dx
!    integral = integrate(nx,evec(:,1)*evec(:,2)*grad(1,:))*dx
!    denom = pot(i,2)-pot(i,1)
!    !
!    IF(i /= 1) force_num(i) = -(pot(i,2)-pot(i-1,2))/(DBLE(14.0_dp/(nx-1)))
!    !
!    IF(ABS(denom) .LT. 1.0d-12 ) THEN
!      nacs(i) = 1.0d4
!    ELSE 
!      nacs(i) = integral / denom
!    END IF
!    IF (i > 1) THEN
!      IF(nacs(i-1)*nacs(i)<0.0_dp) nacs(i)=-nacs(i)
!    END IF
!    WRITE(99,"(7(1x,f13.8))") xn(i), pot(i,1), pot(i,2), nacs(i), force(1,2), force_num(i), force_num(i)/force(1,2)
!  END DO
!  !
!  CLOSE(99)
  !
  !
  DEALLOCATE(ev,evec,nac,laser_coup)
  !
END PROGRAM calculate_BO
