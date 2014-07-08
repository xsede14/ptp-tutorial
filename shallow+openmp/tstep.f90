! Translation of tstep.c into Fortran 2003 (J. Overbey, 17 Sept 2011)
! Refactored and added OpenMP (July 2014)
! Use the Introduce Implicit None refactoring to add declataions for i, j
!************************************************************************
!                                                                       *
! Commonwealth Scientific and Industrial Research Organisation (CSIRO)  *
!        - Division of Information Technology        (DIT)              *
!        - Division of Atmospheric Research        (DAR)                *
!                                                                       *
! Shallow water weather model - Distributed Memory Version              *
!                                                                       *
! Finite difference model of shallow water equations based on :-        *
! "The dynamics of finite difference models of the shallow water        *
! equations" by R. Sadourney, JAS, 32, 1975.                            *
! Code from:-                                                           *
! "An introduction to three-dimensional climate modelling"              *
! by Washington and Parkinson                                           *
!                                                                       *
! Programmers   = David Abramson        (DIT) rcoda@koel.co.rmit.oz     *
!               = Paul Whiting          (DIT) rcopw@koel.co.rmit.oz     *
!               = Martin Dix            (DAR) mrd@koel.co.rmit.oz       *
! Language      = BSD c using Argonne NL macros                         *
! O/S           = Unix System V                                         *
! H/W           = Encore Multimax 320                                   *
!                                                                       *
!************************************************************************

subroutine tstep(m,n,alpha,jstart,jend,cpold,cuold,cvold,cp,cu,cv,cpnew,cunew,cvnew,cdpdt,cdudt,cdvdt,firststep,tdt) bind(c)
  use, intrinsic :: ISO_C_BINDING

  integer(kind=C_INT), value :: m, n
  real(kind=C_FLOAT), value :: alpha
  integer(kind=C_INT), value :: jstart,jend
  type(C_PTR), value :: cpold; real(kind=C_FLOAT), pointer :: pold(:,:)
  type(C_PTR), value :: cuold; real(kind=C_FLOAT), pointer :: uold(:,:)
  type(C_PTR), value :: cvold; real(kind=C_FLOAT), pointer :: vold(:,:)
  type(C_PTR), value :: cp;    real(kind=C_FLOAT), pointer ::    p(:,:)
  type(C_PTR), value :: cu;    real(kind=C_FLOAT), pointer ::    u(:,:)
  type(C_PTR), value :: cv;    real(kind=C_FLOAT), pointer ::    v(:,:)
  type(C_PTR), value :: cpnew; real(kind=C_FLOAT), pointer :: pnew(:,:)
  type(C_PTR), value :: cunew; real(kind=C_FLOAT), pointer :: unew(:,:)
  type(C_PTR), value :: cvnew; real(kind=C_FLOAT), pointer :: vnew(:,:)
  type(C_PTR), value :: cdpdt; real(kind=C_FLOAT), pointer :: dpdt(:,:)
  type(C_PTR), value :: cdudt; real(kind=C_FLOAT), pointer :: dudt(:,:)
  type(C_PTR), value :: cdvdt; real(kind=C_FLOAT), pointer :: dvdt(:,:)
  integer(kind=C_INT), value :: firststep
  real(kind=C_FLOAT), value  :: tdt

  call c_f_pointer(cpold, pold, shape=[m, n])
  call c_f_pointer(cuold, uold, shape=[m, n])
  call c_f_pointer(cvold, vold, shape=[m, n])
  call c_f_pointer(cp,    p,    shape=[m, n])
  call c_f_pointer(cu,    u,    shape=[m, n])
  call c_f_pointer(cv,    v,    shape=[m, n])
  call c_f_pointer(cpnew, pnew, shape=[m, n])
  call c_f_pointer(cunew, unew, shape=[m, n])
  call c_f_pointer(cvnew, vnew, shape=[m, n])
  call c_f_pointer(cdpdt, dpdt, shape=[m, n])
  call c_f_pointer(cdudt, dudt, shape=[m, n])
  call c_f_pointer(cdvdt, dvdt, shape=[m, n])

  !$omp parallel do private(i)
  do j = jstart+1, jend+1
    do i = 1, m
      pnew(i,j) = pold(i,j) + tdt*dpdt(i,j)
      unew(i,j) = uold(i,j) + tdt*dudt(i,j)
      vnew(i,j) = vold(i,j) + tdt*dvdt(i,j)
    end do

    ! Don't apply time filter on first step
    if ( firststep == 0 ) then
      do i = 1, m
        pold(i,j) = p(i,j)+alpha*(pnew(i,j)-2._c_float*p(i,j)+pold(i,j))
        uold(i,j) = u(i,j)+alpha*(unew(i,j)-2._c_float*u(i,j)+uold(i,j))
        vold(i,j) = v(i,j)+alpha*(vnew(i,j)-2._c_float*v(i,j)+vold(i,j))
      end do
    end if

    do i = 1, m
      p(i,j) = pnew(i,j)
      u(i,j) = unew(i,j)
      v(i,j) = vnew(i,j)
    end do
  end do
  !$omp end parallel do
end subroutine
