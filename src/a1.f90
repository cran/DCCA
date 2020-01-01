subroutine inv(n, A, Ainv)
  !-------------------------------------------------------------------------------------------
  ! Returns the inverse of a matrix calculated by finding the LU
  ! decomposition.  Depends on LAPACK.
  !
  !  DGETRF,  DGETRI
  !-------------------------------------------------------------------------------------------
  integer, parameter :: dp = kind(1.d0)
  integer :: n
  real(dp) :: A(n, n)
  real(dp) :: Ainv(n, n)
  real(dp) :: work(n)  ! work array for LAPACK
  integer :: ipiv(n)   ! pivot indices
  integer ::  info

  ! Store A in Ainv to prevent it from being overwritten by LAPACK
  Ainv = A

  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  call DGETRF(n, n, Ainv, n, ipiv, info)

  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  call DGETRI(n, Ainv, n, ipiv, work, n, info)
  return
end subroutine inv


subroutine cumsum(y, r, n)
  !-----------------------------------------
  !
  !  Calculate the integrated process R_t
  !
  !------------------------------------------
  implicit none
  integer, parameter :: dp = kind(1.d0)
  integer :: n
  real(dp) :: y(n), r(n)
  integer :: i
  r(1) = y(1)
  do i = 2,n
     r(i) = r(i-1) + y(i)
  end do
  return
end subroutine cumsum


subroutine Jmatrix(n,J)
  !-----------------------------------------
  !
  !  Matrix J
  !
  !------------------------------------------
  implicit none
  integer, parameter :: dp = kind(1.d0)
  integer :: n
  real(dp) :: J(n,n)
  integer :: r, s
  J = 0.d0
  do r = 1,n
     do s = 1,r
        J(r,s) = 1.d0
     end do
  end do
  return
end subroutine Jmatrix


subroutine Pmatrix1(m, P)
  !-----------------------------------------
  !
  !  Matrix P: projection matrix
  !  used only when v = 0
  !
  !------------------------------------------
  implicit none
  integer, parameter :: dp = kind(1.d0)
  integer :: m
  real(dp) :: P(m+1, m+1)
  integer :: r,s
  do r = 1, ceiling((m+1.0)/2.0)
     do s = r, (m+2-r)
        P(r,s) = (4.d0 - 6.d0*(r+s-1.0)/m +  12.d0*r*s/(m*m + 2.d0*m))/(m + 1.d0)
        P(s,r) = P(r,s)
        P(m+2-r, m+2-s) = P(r,s)
        P(m+2-s, m+2-r) = P(r,s)
     end do
  end do
  return
end subroutine Pmatrix1


subroutine Pmatrix(m, v, P)
  !-----------------------------------------
  !
  !  Matrix P: projection matrix
  !
  !------------------------------------------
  implicit none
  integer, parameter :: dp = kind(1.d0)
  integer :: m, v
  real(dp) :: P(m+1, m+1), D(m+1, v+2)
  real(dp) :: DD(v+2, v+2), DDI(v+2,v+2)
  integer :: r,s
  D(:,1) = 1.d0
  do r = 2,(v+2)
     D(:,r) = (/(1.d0*s**(r-1), s = 1,(m+1))/)
  end do
  DD = matmul(transpose(D), D)
  call inv(v+2, DD, DDI)
  P = matmul(matmul(D, DDI), transpose(D))
  return
end subroutine Pmatrix


subroutine Qm(m, P, Q)
  !------------------------------------------------
  !
  !  Matrix Q: matrix I - P
  !  This subroutine requires P
  !
  !------------------------------------------------
  implicit none
  integer, parameter :: dp = kind(1.d0)
  integer :: m
  real(dp) :: P(m+1, m+1), Q(m+1, m+1)
  integer :: r
  Q = -P
  do r = 1, m+1
     Q(r,r) = 1.d0 - P(r,r)
  end do
  return
end subroutine Qm

subroutine Qmatrix(m, v, Q)
  !-----------------------------------------
  !
  !  Matrix Q: matrix I - P
  !
  !------------------------------------------
  implicit none
  integer, parameter :: dp = kind(1.d0)
  integer :: m, v
  real(dp) :: P(m+1, m+1), Q(m+1, m+1)
  call Pmatrix(m, v, P)
  call Qm(m, P, Q)
  return
end subroutine Qmatrix


subroutine Km(m,J,Q,K)
  !-----------------------------------------
  !
  !  Matrix K = J'QJ
  !  This subrotuine requires J and Q
  !
  !------------------------------------------
  implicit none
  integer, parameter :: dp = kind(1.d0)
  integer :: m
  real(dp) :: J(m+1, m+1), Q(m+1, m+1), K(m+1, m+1)
  K = matmul(matmul(transpose(J), Q),J)
  K(1,:) = 0.d0
  K(:,1) = 0.d0
  return
end subroutine Km



subroutine Kmatrix(m,v,K)
  !-----------------------------------------
  !
  !  Matrix K = J'QJ, for a given m and v
  !
  !------------------------------------------
  implicit none
  integer, parameter :: dp = kind(1.d0)
  integer :: m, v
  real(dp) :: J(m+1, m+1), Q(m+1, m+1), K(m+1, m+1)
  call Qmatrix(m, v, Q)
  call Jmatrix(m+1,J)
  call Km(m, J, Q, K)
  return
end subroutine Kmatrix



subroutine Fmi(m, Ei1, Ei2, fdmi)
  !--------------------------------------------------------------
  !
  !     calculates fdfa_i(m) and/or fdcca_i(m)
  !
  !--------------------------------------------------------------
  implicit none
  integer, parameter :: dp = kind(1.d0)
  integer :: m
  real(dp) :: Ei1(m+1), Ei2(m+1),  fdmi
  fdmi = 1.d0/m*sum(Ei1*Ei2)
  return
end subroutine Fmi


subroutine Fm(m, n, overlap, Q, R1n, R2n, f1, fd1m, f2, fd2m, f12, fd12m, rhom)
  !---------------------------------------------------------------------------------
  !
  !  Generic function to calculate  Fdfa(m)  and  Fdcca(m)
  !
  !  Input:
  !     m: integer indicating the size of the window
  !     n:  sample size
  !     overlap: 0 = NO; 1 = Yes
  !     R1n and R2n: integrated time series
  !     f1, f2, f12: indicate if DFA and DCCA must be calculated
  !                         0 = NO; 1 = YES
  !    fd1m, fd2m, dd12m: F2_DFA AND F_DCCA
  !    rhom: cross correlation coefficient
  !---------------------------------------------------------------------------------
  implicit none
  integer, parameter :: dp = kind(1.d0)
  integer :: m, n, overlap, f1, f2, f12
  real(dp) :: R1n(max(1,n*max(f1,f12))), R2n(max(1,n*max(f2,f12))), Q(m+1,m+1)
  real(dp) :: fd1m, fd2m, fd12m, rhom
  integer :: i, j
  real(dp) :: fd1(max(((n-m)*overlap + n/(m+1)*(1-overlap))*f1,1))
  real(dp) :: fd2(max(((n-m)*overlap + n/(m+1)*(1-overlap))*f2,1))
  real(dp) :: fd12(max(((n-m)*overlap + n/(m+1)*(1-overlap))*f12,1))
  real(dp), allocatable :: E1i(:), E2i(:)
  integer ::  nw

  fd1 = 0.d0
  fd2 = 0.d0
  fd12 = 0.d0
  rhom = 0.d0

  if(allocated(E1i)) deallocate(E1i)
  if(allocated(E2i)) deallocate(E2i)
  if(f1 == 1 .or. f12 == 1) allocate(E1i(m+1))
  if(f2 == 1 .or. f12 == 1) allocate(E2i(m+1))
	
  nw = (n-m)*overlap + n/(m+1)*(1-overlap)

  do j = 1, nw
     i = j
     if(overlap == 0) i = (j-1)*(m+1) + 1
     ! error E
     if(f1 == 1 .or. f12 == 1) E1i = matmul(Q, R1n(i:(i+m)))
     if(f2 == 1 .or. f12 == 1) E2i = matmul(Q, R2n(i:(i+m)))
     ! fdfa(i,m) and fdcca(i,m)
     if(f1 == 1) call Fmi(m, E1i, E1i, fd1(j))
     if(f2 == 1) call Fmi(m, E2i, E2i, fd2(j))
     if(f12 == 1) call Fmi(m, E1i, E2i, fd12(j))
  end do

  if(f1 == 1) fd1m = 1.d0/nw*sum(fd1)
  if(f2 == 1) fd2m = 1.d0/nw*sum(fd2)
  if(f12 == 1) fd12m = 1.d0/nw*sum(fd12)
  if(f1*f2*f12 == 1) rhom = fd12m/sqrt(fd1m*fd2m)
  return
end subroutine Fm


subroutine dfadcca(lm, m, v, overlap, n, Y1n, Y2n, f1, fd1m, f2, fd2m, f12, fd12m, rhom)
  implicit none
  integer, parameter :: dp = kind(1.d0)
  integer :: lm, n, overlap, f1, f2, f12, v
  integer :: m(lm)
  real(dp) :: Y1n(max(1, n*max(f1,f12))), Y2n(max(1, n*max(f2,f12)))
  real(dp) :: rhom(max(lm*f1*f2*f12,1))
  real(dp) :: fd1m(max(lm*f1,1)), fd2m(max(lm*f2,1)), fd12m(max(lm*f12,1))
  integer :: im, im1, im2, im12, imr, n1, n2
  real(dp), allocatable :: Q(:,:), R1n(:), R2n(:)

  n1 = max(1, n*max(f1,f12))
  n2 = max(1, n*max(f2,f12))
  if(allocated(R1n)) deallocate(R1n)
  if(allocated(R2n)) deallocate(R2n)  
  allocate(R1n(n1), R2n(n2))
  R1n = 0.d0
  R2n = 0.d0
  call cumsum(Y1n, R1n, n1)
  call cumsum(Y2n, R2n, n2)
  
  do im = 1, lm
     if(allocated(Q)) deallocate(Q)
     allocate(Q(m(im)+1, m(im)+1))
     call Qmatrix(m(im), v, Q )
     im1 = max(im*f1,1)
     im2 = max(im*f2,1)
     im12 = max(im*f12,1)
     imr  = max(im*f1*f2*f12,1)
     call Fm(m(im), n, overlap, Q, R1n, R2n, f1, fd1m(im1), f2, fd2m(im2), f12, fd12m(im12), rhom(imr))
  end do
  return
end subroutine dfadcca


subroutine trace(n, A, tr)
  !-----------------------------------------
  !
  !  Trace function
  !
  !-----------------------------------------
  implicit none
  integer, parameter :: dp = kind(1.d0)
  integer :: n, i
  real(dp) :: A(n,n), tr
  tr = A(1,1)
  do i = 2,n
     tr = tr + A(i,i)
  end do
  return
end subroutine trace

subroutine Em(m, K, G, E)
  !-----------------------------------------
  !
  !  E(F). Depends on K.
  !
  !-----------------------------------------
  implicit none
  integer, parameter :: dp = kind(1.d0)
  integer :: m
  real(dp) :: K(m+1, m+1), G(m+1, m+1), E
  call trace(m+1, matmul(K,G), E)
  E = 1.d0/m*E
  return
end subroutine Em

subroutine Expectm(m,v,G,E)
  !-----------------------------------------
  !
  !  E(F). Depends on v.
  !     K is calculated by the subroutine
  !
  !-----------------------------------------
  implicit none
  integer, parameter :: dp = kind(1.d0)
  integer :: m,v
  real(dp) :: K(m+1, m+1), G(m+1, m+1), E
  call Kmatrix(m,v,K)
  call Em(m, K, G, E)
  return
end subroutine Expectm


subroutine Expectms(lm, m, v, c1, G1, c2, G2, c12, G12, E1, E2, E12, rhos)
  implicit none
  integer, parameter :: dp = kind(1.d0)
  integer :: lm, m(lm), c1, c2, c12, v
  real(dp) :: G1(m(lm)*c1+1, m(lm)*c1+1), E1(max(lm*c1,1))
  real(dp) :: G2(m(lm)*c2+1, m(lm)*c2+1), E2(max(lm*c2,1))
  real(dp) :: G12(m(lm)*c12+1, m(lm)*c12+1), E12(max(lm*c12,1))
  real(dp) :: rhos(max(lm*c1*c2*c12, 1))
  real(dp), allocatable ::  K(:,:)
  integer :: i

  E1 = 0.d0
  E2 = 0.d0
  E12 = 0.d0

  do i = 1, lm
     if(allocated(k)) deallocate(k)
     allocate(k(m(i)+1, m(i)+1))
     call Kmatrix(m(i), v, K)
     if(c1 == 1) call trace(m(i)+1, matmul(K,G1(1:(m(i)+1), 1:(m(i)+1))), E1(i))
     if(c2 == 1) call trace(m(i)+1, matmul(K,G2(1:(m(i)+1), 1:(m(i)+1))), E2(i))
     if(c12 == 1) call trace(m(i)+1, matmul(K,G12(1:(m(i)+1), 1:(m(i)+1))), E12(i))
  end do

  if(c1 == 1) E1 = 1.d0/m*E1
  if(c2 == 1) E2 = 1.d0/m*E2
  if(c12 == 1) E12 = 1.d0/m*E12
  if(c1*c2*c12 == 1) rhos = E12/sqrt(E1*E2)
  return
end subroutine Expectms



subroutine Kroenecker(nlA, ncA, A, nlB, ncB, B, AB)
  !-----------------------------------------
  !
  !  Kroencker Product
  !
  !-----------------------------------------
  implicit none
  integer, parameter :: dp = kind(1.d0)
  integer :: nlA, nlB, ncA, ncB
  real(dp) :: A(nlA, ncA), B(nlB, ncB), AB(nlA*nlB, ncA*ncB)
  integer :: r,s

  do r = 1, nlA
     do s = 1,ncA
        !Place a block of B values.
        AB(((r-1)*nlB + 1):(r*nlB),((s-1)*ncB + 1):(s*ncB)) = A(r,s)*B
     end do
  end do

  return
end subroutine Kroenecker
	
subroutine Kkronm(m,h,overlap,K,Kkron)
  !-----------------------------------------
  !
  !  Matrix Kkron
  !  J and Q must be provided.
  !
  !  (m+1+h)*overlap + (m+1)*(h+1)*(1-overlap)
  !  is the same as
  !  (m+1)*(h+1) - m*h*overlap
  !------------------------------------------
  implicit none
  integer, parameter :: dp = kind(1.d0)
  integer :: m, h, overlap
  real(dp) :: K(m+1,m+1)
  real(dp) :: Kkron(((m+1)*(h+1) - m*h*overlap)*(m+1),((m+1)*(h+1) - m*h*overlap)*(m+1))
  real(dp) :: Kh((m+1)*(h+1) - m*h*overlap, (m+1)*(h+1) - m*h*overlap)
  integer :: mm

  mm = (m+1)*(h+1) - m*h*overlap
  Kh = 0.d0
  Kh((mm-m):mm,(mm-m):mm) = K
  call Kroenecker(m+1,m+1,K,mm,mm,Kh,Kkron)
  return
end subroutine Kkronm


subroutine Kkronmatrix(m,h,v,overlap,Kkron)
  !-----------------------------------------
  !
  !  Matrix Kkron
  !  J and Q are calculated
  !
  !  (m+1+h)*overlap + (m+1)*(h+1)*(1-overlap)
  !  is the same as
  !  (m+1)*(h+1) - m*h*overlap
  !------------------------------------------
  implicit none
  integer, parameter :: dp = kind(1.d0)
  integer :: m, h, v, overlap
  real(dp) :: K(m+1, m+1)
  real(dp) :: Kkron(((m+1)*(h+1) - m*h*overlap)*(m+1),((m+1)*(h+1) - m*h*overlap)*(m+1))
  call Kmatrix(m,v,K)
  call Kkronm(m, h, overlap, K, Kkron)
  return
end subroutine Kkronmatrix


subroutine Covamh(m, h, overlap, K, Kkron, G1, G2, G12, nrc, Cum, df, cova)
  !-----------------------------------------
  !
  !  Covariance of DFA and DCCA
  !
  !  (m+1+h)*overlap + (m+1)*(h+1)*(1-overlap)
  !  is the same as
  !  mm = (m+1)*(h+1) - m*h*overlap
  !
  !  K = (m+1) by (m+1)
  !  Kh = (m+1)mm by (m+1)mm
  !  Kkron = K %x% Kh
  !------------------------------------------
  implicit none
  integer, parameter :: dp = kind(1.d0)
  integer :: m, h, overlap, df, nrc
  real(dp) :: K(m+1, m+1)
  real(dp) :: Kkron(max(1,nrc),max(1,nrc))
  real(dp) :: G1(m+1,m+1), G2(m+1,m+1)
  real(dp) :: G12((m+1)*(h+1) - m*h*overlap, (m+1)*(h+1) - m*h*overlap)
  real(dp) :: Cum(max(nrc,1),max(nrc,1))
  real(dp) :: cova, p1, p2, p3
  integer :: mm
  mm = (m+1)*(h+1) - m*h*overlap
  p1 = 0.d0
  if(nrc > 0) call trace((m+1)*mm, matmul(Kkron,Cum), p1)
  p2 = 0.d0
  call trace(m+1, matmul(matmul(K,G1),matmul(K,G2)), p2)
  p3 = p2
  if(df == 0) call trace(m+1, matmul(matmul(K,G12(1:(m+1),(mm-m):mm)),matmul(K,G12((mm-m):mm, 1:(m+1)))), p3)
  cova = 1.d0/m**2*(p1+p2+p3)
  return
end subroutine Covamh


subroutine CovF2dfa(lm, m, v, lh, h, overlap, nr1, nc1, G1, nrc, Cum, df, cova)
  !-----------------------------------------
  !
  !  Covariance of DFA and DCCA
  !
  !  (m+1+h)*overlap + (m+1)*(h+1)*(1-overlap)
  !  is the same as
  !  mm = (m+1)*(h+1) - m*h*overlap
  !------------------------------------------
  implicit none
  integer, parameter :: dp = kind(1.d0)
  integer :: lm, lh, v, overlap, df
  integer :: m(lm), h(lh)
  integer :: nr1, nc1, nrc
  real(dp) :: cova(lm, lh)
  real(dp) :: G1(nr1,nc1)
  real(dp) :: Cum(max(nrc,1), max(nrc,1))
  real(dp), allocatable :: K(:,:), Kkron(:,:), Kh(:,:)
  real(dp), allocatable :: Q(:,:), J(:,:)
  integer :: mm, im, ih, mi, hi, nk

  cova = 0.d0

  ! in case the cumulants are zero:
  if(allocated(Kkron)) deallocate(Kkron)
  allocate(kkron(1,1))
  kkron = 0.d0
  nk = 1

  mm =  (maxval(m)+1)*(maxval(h)+1) - maxval(m)*maxval(h)*overlap
  if(allocated(J)) deallocate(J)
  allocate(J(mm,mm))
  call Jmatrix(mm, J)

  do im = 1,lm
     mi = m(im)
     if(allocated(Q)) deallocate(Q)
     allocate(Q(mi+1, mi+1))
     call Qmatrix(mi,v,Q)

     if(allocated(K)) deallocate(K)
     allocate(K(mi+1, mi+1))
     call Km(mi,J(1:(mi+1),1:(mi+1)),Q, K)

     do ih = 1,lh
        hi = h(ih)
        mm =  (mi+1)*(hi+1) - mi*hi*overlap        

        if(nrc > 0) then
           if(allocated(Kh)) deallocate(Kh)
           allocate(Kh(mm,mm))
           Kh = matmul(matmul(transpose(J((mm-mi):mm,1:mm)),Q),J((mm-mi):mm,1:mm))
           nk = mm*(mi+1)
           if(allocated(Kkron)) deallocate(Kkron)
           allocate(kkron(nk,nk))
           call Kroenecker(mi+1,mi+1,K,mm,mm,Kh,Kkron)
        end if

        call Covamh(mi, hi, overlap, K, Kkron, G1(1:(mi+1),(mm-mi):mm),&
             G1((mm-mi):mm, 1:(mi+1)), G1(1:mm,1:mm), nrc, Cum(1:nk,1:nk),&
             df,cova(im,ih))
     end do
  end do
  return
end subroutine CovF2dfa

subroutine CovFdcca(lm, m, v, lh, h, overlap, nr1, nc1, G1, nr2, nc2, G2, &
     nr12, nc12, G12, nrc, Cum, df, cova)
  !-----------------------------------------
  !
  !  Covariance of DFA and DCCA
  !
  !  (m+1+h)*overlap + (m+1)*(h+1)*(1-overlap)
  !  is the same as
  !  mm = (m+1)*(h+1) - m*h*overlap
  !------------------------------------------
  implicit none
  integer, parameter :: dp = kind(1.d0)
  integer :: lm, lh, v, overlap, df
  integer :: m(lm), h(lh)
  integer :: nr1, nc1,nr2, nc2, nr12, nc12, nrc
  real(dp) :: cova(lm, lh)
  real(dp) :: G1(nr1,nc1), G2(nr2, nc2)
  real(dp) :: G12(nr12, nc12), Cum(max(nrc,1), max(nrc,1))
  real(dp), allocatable :: K(:,:), Kh(:,:), Kkron(:,:)
  real(dp), allocatable :: Q(:,:), J(:,:)
  integer :: mm, im, ih, mi, hi, nk

  cova = 0.d0

  ! in case the cumulants are zero:
  if(allocated(Kkron)) deallocate(Kkron)
  allocate(kkron(1,1))
  kkron = 0.d0
  nk = 1

  mm =  (maxval(m)+1)*(maxval(h)+1) - maxval(m)*maxval(h)*overlap
  if(allocated(J)) deallocate(J)
  allocate(J(mm,mm))
  call Jmatrix(mm, J)

  do im = 1,lm
     mi = m(im)
     if(allocated(Q)) deallocate(Q)
     allocate(Q(mi+1, mi+1))
     call Qmatrix(mi,v,Q)

     if(allocated(K)) deallocate(K)
     allocate(K(mi+1, mi+1))
     call Km(mi,J(1:(mi+1),1:(mi+1)),Q, K)

     do ih = 1,lh
        hi = h(ih)
        mm =  (mi+1)*(hi+1) - mi*hi*overlap

        if(nrc > 0) then
           if(allocated(Kh)) deallocate(Kh)
           allocate(Kh(mm,mm))
           Kh = matmul(matmul(transpose(J((mm-mi):mm,1:mm)),Q),J((mm-mi):mm,1:mm))
           nk = mm*(mi+1)
           if(allocated(Kkron)) deallocate(Kkron)
           allocate(kkron(nk,nk))
           call Kroenecker(mi+1,mi+1,K,mm,mm,Kh,Kkron)
        end if

        call Covamh(mi, hi, overlap, K, Kkron, G1(1:(mi+1),(mm-mi):mm),&
             G2((mm-mi):mm, 1:(mi+1)), G12(1:mm,1:mm), nrc, Cum(1:nk,1:nk),&
             df,cova(im,ih))
     end do
  end do
  return
end subroutine CovFdcca
