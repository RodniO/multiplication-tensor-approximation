
function ExampleTEN(ver,nm_int,mms,r,ret_st,dis_ret_st,maxst,dis_t,dis_start_b,stop_b,mask_c,mat_s,reg_s,reg_d_s) result(res)
  USE ModMtrx
  USE MerTwist, only : grnd !RNG
  
  Integer(4), intent(in) :: ver !verbose: how much log on screen is needed (0, 1 or 2)
  Integer(4), intent(in) :: nm_int !Integer to add in the end of file names
  Integer(4), dimension(3), intent(in) :: mms !Matrix sizes
  Integer(4), intent(in) :: r !Desired rank
  Integer(4), intent(in) :: ret_st !Restart is done at step=ret_st if error>1
  Integer(4), intent(in) :: dis_ret_st !Restart is done at step=dis_ret_st if discretization has not yet started
  Integer(4), intent(in) :: maxst !Maximum total number of steps allowed
  Double precision, intent(in) :: dis_t !Maximum time spent on trying to fit leftover (not discretized) coefficients
  Double precision, intent(in) :: dis_start_b !When error<dis_start_b is reached, discretization starts
  Double precision, intent(in) :: stop_b !When error<stop_b is reached, discretization stops
  Double precision, intent(in) :: mask_c !Only elements>mask_c are fixed to be exactly 1
  Double precision, intent(in), optional :: mat_s !Initial scale of u, v and w matrices
  Double precision, intent(in), optional :: reg_s !Initial scale of regularization after ret_st steps
  Double precision, intent(in), optional :: reg_d_s !Initial scale of regularization at discretization start
  Integer(4) :: res !2 - success, 1 - near success (small error, but not discretized), 0 - fail
  
  Integer(4) i, j, k, l !Indices
  
  Double precision, dimension(mms(1)*mms(2),mms(2)*mms(3),mms(3)*mms(1)) :: mmat, emat !Multiplication tensor and error
  Type(Mtrx) :: u, v, w !Low-rank factors. W is transposed for symmetry reasons
  Type(Mtrx) :: uu, uv, uw, u1, v1, w1 !Copies of low-rank factors
  Type(Mtrx) masku, maskv, maskw, maskzu, maskzv, maskzw !Masks for fixed elements
  Double precision scal !Scale of u, v and w
  Integer(4) ms(3) !Number of columns in u, v and w
  
  Integer(4) step !Current step
  Integer(4) addstep !Steps after adding last fixed element
  Integer(4) screenlogsteps !Number of steps between screen logs
  Double precision err, err0 !Error calculation results
  Double precision alpha, alphacoef, alphadcoef !Regularization coefficient and its scales
  Integer(4) nonzeros !Number of expected fixed elements (only nonzeros are used by default)
  Integer(4) fixed !Number of fixed elements
    
  Double precision dsecnd !CPU time evaluation from LAPACK
  Double precision tmptime !Temporary to save CPU time
  
  !TO DO
  Double precision maskcoef !REPLACE WITH JUST SEARCHING FOR MAXIMUM; REPLACE u, v, w with u(3)
  Integer(4) jkl(mms(1))
  Integer(4), parameter :: mats = 0
  Type(Mtrx) :: a(mats), b(mats), c(mats), ait(mats), bit(mats), cit(mats)
  Type(Mtrx) :: al(mats), bl(mats), cl(mats)

  Double precision x, y, z !Maximums in the arrays
  Integer(4) ij(2) !Maximums positions
  Logical disstarted !Determines, whether discretization has started
  
  ms = [mms(1)*mms(2), mms(2)*mms(3), mms(3)*mms(1)]
  
  screenlogsteps = 100
  
  scal = (2.0d0/r)**(1.5d0)/sqrt(dble(ms(1)*ms(2)*ms(3)))
  if (present(mat_s)) then
    scal = mat_s**3
  end if
  alphacoef = 0.5d0
  if (present(reg_s)) then
    alphacoef = reg_s
  end if
  alphadcoef = 2.0d0
  if (present(reg_d_s)) then
    alphadcoef = reg_d_s
  end if
  
  alpha = 0
  fixed = 0
  disstarted = .false.
  nonzeros = 1
  
  call maskzu%init(r,ms(1))
  call maskzv%init(r,ms(2))
  call maskzw%init(r,ms(3))
  
  do i = 1, mats
    call a(i)%gauss(mms(1), mms(1))
    call b(i)%gauss(mms(1), mms(1))
    call c(i)%gauss(mms(1), mms(1))
    ait(i) = .T.(.I. a(i))
    bit(i) = .T.(.I. b(i))
    cit(i) = .T.(.I. c(i))
    al(i) = kron(a(i), bit(i))
    bl(i) = kron(b(i), cit(i))
    cl(i) = kron(c(i), ait(i))
  end do
  
  mmat(:,:,:) = 0
  do i = 1, mms(1)
    do j = 1, mms(2)
      do k = 1, mms(3)
        mmat(mms(2)*(i-1)+j,mms(3)*(j-1)+k,mms(1)*(k-1)+i) = 1
      end do
    end do
  end do
  
  call u%gauss(r,ms(1))
  call v%gauss(r,ms(2))
  call w%gauss(r,ms(3))
  u = u*scal**(1.0d0/3.0d0)
  v = v*scal**(1.0d0/3.0d0)
  w = w*scal**(1.0d0/3.0d0)
  
  masku%n = 0
  
  !CAN INITIALISE MASKS, IF ONE REALLY WANTS
  !u%d(1,1:ms(1)) = [1,0,0,0,1,0,0,0,0]
  !call masku%init(r,ms(1))
  !call maskv%init(r,ms(2))
  !call maskw%init(r,ms(3))
  !masku%d(1,1:ms(1)) = [1,0,0,0,1,0,0,0,0]
  !maskzu%d(1,1:ms(1)) = [0,1,1,1,0,1,1,1,1]
  !alpha = 0.1d0**14 
  
  if (ver > 0) then
    print *, 'ITERATIONS START'
  end if

  step = 0
  addstep = 0
  tmptime = dsecnd()
  do while(step < maxst)
    step = step + 1
    addstep = addstep + 1
        
    !Alternative Least Squares
    if (masku%n > 0) then
      call ALS(u, v, w, ms, mmat, alpha, masku, maskv, maskw, maskzu, maskzv, maskzw)
    else
      call ALS(u, v, w, ms, mmat, alpha)
    end if
    
    !REGULARIZATION UPDATE
    if (mod(step,10) == 1) then
      emat(:,:,:) = mmat(:,:,:)
      do i = 1, r
        do j = 1, ms(1)
          do k = 1, ms(2)
            do l = 1, ms(3)
              emat(j,k,l) = emat(j,k,l) - u%d(i,j)*v%d(i,k)*w%d(i,l)  
            end do
          end do
        end do
      end do
      emat(:,:,:) = emat(:,:,:)**2
      err = sqrt(sum(emat))
      alpha = min(alpha, err*alphacoef)
      if (disstarted) then
        alpha = min(alpha, alphadcoef*err/sqrt(r*dble(ms(1)+ms(2)+ms(3))-fixed))
      end if
      if (err < stop_b) then
        exit
      end if
    end if
    
    !SCREEN LOG
    if ((mod(step,screenlogsteps) == 0) .and. (ver == 2)) then
      print *, 'STEP', step
      print *, 'ERROR', err
      print *, 'Alpha', alpha
      print *, 'Time passed', dsecnd() - tmptime
      tmptime = dsecnd()
    end if
    
    !ADD ELEMENTS WHEN LOCAL MINIMUM
    !if ((addstep==dis_ret_st).and.(err>dis_start_b).and.(err<1.0d0).and.(disstarted).and.(fixed==nonzeros)) then
    !  step = step - dis_ret_st/10
    !  addstep = dis_ret_st*9/10
    !  nonzeros = nonzeros + 1
    !end if
    
    !RETRY IF BAD ACCURACY
    if (((step == ret_st) .and. (err > 1.0d0)) .or. ((addstep == dis_ret_st) .and. (err > dis_start_b))) then
      !UNCOMMENT TO STOP INSTEAD OF RETRY
      !if (err < 1.0d0) then
      !  exit
      !end if
      if (ver == 2) then
        print *, 'RETRY'
      end if
      call u%gauss(r,ms(1))
      call v%gauss(r,ms(2))
      call w%gauss(r,ms(3))
      u = u*scal**(1.0d0/3.0d0)
      v = v*scal**(1.0d0/3.0d0)
      w = w*scal**(1.0d0/3.0d0)
      step = 0
      alpha = 0
      addstep = 0
      nonzeros = 1
      masku%n = 0
      maskzu%d = 0
      maskzv%d = 0
      maskzw%d = 0
      if (disstarted) then
        disstarted = .false.
        call masku%deinit()
        call maskv%deinit()
        call maskw%deinit()
      end if
    end if
    
    !INITIALISE ALPHA IF GOOD ACCURACY
    if ((step == ret_st) .and. (err <= 1.0d0)) then
      alpha = sqrt(1.0d0/(u%fnorm()**2 + v%fnorm()**2 + w%fnorm()**2))
      alpha = min(alpha, err*alphacoef)
      if (ver == 2) then
        print *, 'Alpha (initial)', alpha
      end if
    end if
    
    !START DISCRETE CONVERGENCE
    if ((err < dis_start_b) .and. (.not.disstarted) .and. (mod(step,10) == 0)) then
      disstarted = .true.
      alpha = err*alphacoef
      addstep = 0
      call masku%init(r,ms(1))
      call maskv%init(r,ms(2))
      call maskw%init(r,ms(3))
    end if
    
    !INCREASE NUMBER OF NONZEROS
    if ((disstarted) .and. (mod(step,100) == 0)) then
      if (err < dis_start_b) then
        addstep = 0
        nonzeros = nonzeros + 1
      end if
      if (err < sqrt(stop_b*dis_start_b)) then
        nonzeros = fixed+1
      end if
      call uu%copy(u)
      call uv%copy(v)
      call uw%copy(w)
      uu = uu - (u .dot. masku)
      uv = uv - (v .dot. maskv)
      uw = uw - (w .dot. maskw)
      uu%d(:,:) = abs(uu%d(:,:))
      uv%d(:,:) = abs(uv%d(:,:))
      uw%d(:,:) = abs(uw%d(:,:))
      do while (masku%fnorm()**2+maskv%fnorm()**2+maskw%fnorm()**2+maskzu%fnorm()**2+maskzv%fnorm()**2+maskzw%fnorm()**2<nonzeros)
        x = maxval(uu%d)
        y = maxval(uv%d)
        z = maxval(uw%d)
        if (max(x, max(y,z)) < mask_c) then
          exit
        end if
        if (x > max(y, z)) then
          ij = maxloc(uu%d)
          masku%d(ij(1),ij(2)) = 1.0d0
          u%d(ij(1),ij(2)) = sign(1.0d0,u%d(ij(1),ij(2)))
          uu%d(ij(1),ij(2)) = 0
        else if (y > z) then
          ij = maxloc(uv%d)
          maskv%d(ij(1),ij(2)) = 1.0d0
          v%d(ij(1),ij(2)) = sign(1.0d0,v%d(ij(1),ij(2)))
          uv%d(ij(1),ij(2)) = 0
        else
          ij = maxloc(uw%d)
          maskw%d(ij(1),ij(2)) = 1.0d0
          w%d(ij(1),ij(2)) = sign(1.0d0,w%d(ij(1),ij(2)))
          uw%d(ij(1),ij(2)) = 0
        end if
        !UNCOMMENT TO PUT ZEROS IN ZERO MASKS
        !do i = 1, r
          !do j = 1, ms(1)
            !if (abs(u%d(i,j)) < 0.01d0) then
            !  maskzu%d(i,j) = 1.0d0
            !  masku%d(i,j) = 0
            !  u%d(i,j) = 0
            !end if
          !end do
          !do j = 1, ms(2)
            !if (abs(v%d(i,j)) < 0.01d0) then
            !  maskzv%d(i,j) = 1.0d0
            !  maskv%d(i,j) = 0
            !  v%d(i,j) = 0
            !end if
          !end do
          !do j = 1, ms(3)
            !if (abs(w%d(i,j)) < 0.01d0) then
            !  maskzw%d(i,j) = 1.0d0
            !  maskw%d(i,j) = 0
            !  w%d(i,j) = 0
            !end if
          !end do
        !end do
      end do
      fixed = nint(masku%fnorm()**2+maskv%fnorm()**2+maskw%fnorm()**2+maskzu%fnorm()**2+maskzv%fnorm()**2+maskzw%fnorm()**2)
      if (err < dis_start_b) then
        emat(:,:,:) = mmat(:,:,:)
        do i = 1, r
          do j = 1, ms(1)
            do k = 1, ms(2)
              do l = 1, ms(3)
                emat(j,k,l) = emat(j,k,l) - u%d(i,j)*v%d(i,k)*w%d(i,l)  
              end do
            end do
          end do
        end do
        emat(:,:,:) = emat(:,:,:)**2
        err = sqrt(sum(emat))
        alpha = alphadcoef*err/sqrt(r*dble(ms(1)+ms(2)+ms(3))-fixed)
      end if
      if (ver == 2) then
        print *, 'NONZEROS FOUND', nint(masku%fnorm()**2), nint(maskv%fnorm()**2), nint(maskw%fnorm()**2)
        print *, 'ZEROS FOUND', nint(maskzu%fnorm()**2), nint(maskzv%fnorm()**2), nint(maskzw%fnorm()**2)
      end if
    end if
  end do
  
  if (ver > 0) then
    print *, 'ITERATIONS FINISHED'
  end if
  
  call u1%copy(u)
  call v1%copy(v)
  call w1%copy(w)
  
  !ROUND DOWN ELEMENTS > 1
  do i = 1, r
    do j = 1, ms(1)
      if (abs(u%d(i,j)) >= 1) then
        u%d(i,j) = sign(1.0d0,u%d(i,j))
      end if
    end do
    do j = 1, ms(2)
      if (abs(v%d(i,j)) >= 1) then
        v%d(i,j) = sign(1.0d0,v%d(i,j))
      end if
    end do
    do j = 1, ms(3)
      if (abs(w%d(i,j)) >= 1) then
        w%d(i,j) = sign(1.0d0,w%d(i,j))
      end if
    end do
  end do
  
  !DISCRETIZATION SEARCH
  call uu%init(r,ms(1))
  call uv%init(r,ms(2))
  call uw%init(r,ms(3))
  tmptime = dsecnd()
  err = 1
  do while (err > 0.00001d0)
    do i = 1, r
      do j = 1, ms(1)
        if (grnd() > abs(u%d(i,j))) then
          uu%d(i,j) = 0
        else
          uu%d(i,j) = sign(1.0d0,u%d(i,j))
        end if
      end do
      do j = 1, ms(2)
        if (grnd() > abs(v%d(i,j))) then
          uv%d(i,j) = 0
        else
          uv%d(i,j) = sign(1.0d0,v%d(i,j))
        end if
      end do
      do j = 1, ms(3)
        if (grnd() > abs(w%d(i,j))) then
          uw%d(i,j) = 0
        else
          uw%d(i,j) = sign(1.0d0,w%d(i,j))
        end if
      end do
    end do
    emat(:,:,:) = mmat(:,:,:)
    do i = 1, r
      do j = 1, ms(1)
        do k = 1, ms(2)
          do l = 1, ms(3)
            emat(j,k,l) = emat(j,k,l) - uu%d(i,j)*uv%d(i,k)*uw%d(i,l)  
          end do
        end do
      end do
    end do
    emat(:,:,:) = emat(:,:,:)**2
    err = sqrt(sum(emat))
    if (err > 0.00001d0) then
    do i = 1, r
      do j = 1, ms(1)
        if (grnd() > sqrt(abs(u%d(i,j)))) then
          uu%d(i,j) = 0
        else
          uu%d(i,j) = sign(1.0d0,u%d(i,j))
        end if
      end do
      do j = 1, ms(2)
        if (grnd() > sqrt(abs(v%d(i,j)))) then
          uv%d(i,j) = 0
        else
          uv%d(i,j) = sign(1.0d0,v%d(i,j))
        end if
      end do
      do j = 1, ms(3)
        if (grnd() > sqrt(abs(w%d(i,j)))) then
          uw%d(i,j) = 0
        else
          uw%d(i,j) = sign(1.0d0,w%d(i,j))
        end if
      end do
    end do
    emat(:,:,:) = mmat(:,:,:)
    do i = 1, r
      do j = 1, ms(1)
        do k = 1, ms(2)
          do l = 1, ms(3)
            emat(j,k,l) = emat(j,k,l) - uu%d(i,j)*uv%d(i,k)*uw%d(i,l)  
          end do
        end do
      end do
    end do
    emat(:,:,:) = emat(:,:,:)**2
    err = sqrt(sum(emat))
    end if
    if (dsecnd() - tmptime > dis_t) then
      exit
    end if
  end do
  
  emat(:,:,:) = mmat(:,:,:)
  do i = 1, r
    do j = 1, ms(1)
      do k = 1, ms(2)
        do l = 1, ms(3)
          emat(j,k,l) = emat(j,k,l) - u1%d(i,j)*v1%d(i,k)*w1%d(i,l)  
        end do
      end do
    end do
  end do
  emat(:,:,:) = emat(:,:,:)**2
  err0 = sqrt(sum(emat))
  if (err < 0.00001d0) then
    u%d(:,:) = uu%d(:,:)
    v%d(:,:) = uv%d(:,:)
    w%d(:,:) = uw%d(:,:)
    if (ver > 0) then
      print *, 'TOTAL NONZEROS', nint(u%fnorm()**2+v%fnorm()**2+w%fnorm()**2)
      print *, 'ERROR', err
    end if
  else
    call u%copy(u1)
    call v%copy(v1)
    call w%copy(w1)
    if (ver > 0) then
      print *, 'ERROR', err0
    end if
  end if
  
  !SAVE MATRICES
  if (err < 0.00001d0) then
    open(unit=1, file='out'//itoa(nm_int)//'.r4', form='unformatted', access='direct', recl=r*max(ms(1),ms(2),ms(3))*8)
    res = 2
  else if (err0 < stop_b) then
    open(unit=1, file='outnear'//itoa(nm_int)//'.r4', form='unformatted', access='direct', recl=r*max(ms(1),ms(2),ms(3))*8)
    res = 1
  else
    open(unit=1, file='outbad'//itoa(nm_int)//'.r4', form='unformatted', access='direct', recl=r*max(ms(1),ms(2),ms(3))*8)
    res = 0
  end if
  write(1, rec=1) u%d
  write(1, rec=2) v%d
  write(1, rec=3) w%d
  
  !PRINT MATRICES
  if (ver > 1) then
    do i = 1, r
      print *, 'MATRIX', i
      write(*,"('[', "//itoa(mms(2))//"F6.1)") u%d(i,:mms(2))
      do j = 2, mms(1)-1
        write(*,"("//itoa(mms(2))//"F6.1)") u%d(i,mms(2)*(j-1)+1:mms(2)*j)
      end do
      write(*,"("//itoa(mms(2))//"F6.1, ']')") u%d(i,mms(2)*(mms(1)-1)+1:)
      write(*,"('[', "//itoa(mms(3))//"F6.1)") v%d(i,:mms(3))
      do j = 2, mms(2)-1
        write(*,"("//itoa(mms(3))//"F6.1)") v%d(i,mms(3)*(j-1)+1:mms(3)*j)
      end do
      write(*,"("//itoa(mms(3))//"F6.1, ']')") v%d(i,mms(3)*(mms(2)-1)+1:)
      write(*,"('[', "//itoa(mms(1))//"F6.1)") w%d(i,:mms(1))
      do j = 2, mms(3)-1
        write(*,"("//itoa(mms(1))//"F6.1)") w%d(i,mms(1)*(j-1)+1:mms(1)*j)
      end do
      write(*,"("//itoa(mms(1))//"F6.1, ']')") w%d(i,mms(1)*(mms(3)-1)+1:)
    end do
  end if
  
  close(1)
end

!RESCALE, DOES NOT INCREASE TOTAL FROBENIUS NORM
subroutine rescale(u, v, w)
  USE ModMtrx
  Type(Mtrx) :: u, v, w
  Integer(4) i
  Double precision x, y, z
  
  do i = 1, u%n
    x = sum(u%d(i,:)**2)
    y = sum(v%d(i,:)**2)
    z = sum(w%d(i,:)**2)
    u%d(i,:) = u%d(i,:)*(y/x)**(1.0d0/6.0d0)*(z/x)**(1.0d0/6.0d0)
    v%d(i,:) = v%d(i,:)*(x/y)**(1.0d0/6.0d0)*(z/y)**(1.0d0/6.0d0)
    w%d(i,:) = w%d(i,:)*(x/z)**(1.0d0/6.0d0)*(y/z)**(1.0d0/6.0d0)
  end do
end

!Kronecker product
function kron(a, b) result(res)
  USE ModMtrx
  Type(Mtrx), intent(in) :: a, b
  Type(Mtrx) :: res
  
  Integer(4) i1, i2, j1, j2, m1, m2, n1, n2
  
  m1 = a%n
  m2 = b%n
  n1 = a%m
  n2 = b%m
  call res%init(m1*m2, n1*n2)

  do i2 = 1, n1
    do j2 = 1, n2
      do i1 = 1, m1
        do j1 = 1, m2
          res%d((i1-1)*m2 + j1, (i2-1)*n2 + j2) = a%d(i1, i2) * b%d(j1, j2)
        end do
      end do
    end do
  end do
end

!ALS step
subroutine ALS(u, v, w, ms, mmat, alpha, masku, maskv, maskw, maskzu, maskzv, maskzw)
  USE ModMtrx
  Type(Mtrx) :: u, v, w
  Integer(4), intent(in) :: ms(3)
  Double precision, dimension(ms(1),ms(2),ms(3)), intent(in) :: mmat
  Double precision, intent(in) :: alpha
  Type(Mtrx), intent(in), optional :: masku, maskv, maskw, maskzu, maskzv, maskzw
  
  Integer(4) i, j, k, r
  Type(Mtrx) lhs(3), rhs(3), mid(3)
  
  r = u%n
  
  if ((alpha == 0) .and. (.not.present(masku))) then
    call rhs(1)%init(1,ms(2)*ms(3))
    call mid(1)%init(r,ms(2)*ms(3))
    call rhs(2)%init(1,ms(1)*ms(3))
    call mid(2)%init(r,ms(1)*ms(3))
    call rhs(3)%init(1,ms(1)*ms(2))
    call mid(3)%init(r,ms(1)*ms(2))
  else
    call rhs(1)%init(1,ms(2)*ms(3)+r)
    call mid(1)%init(r,ms(2)*ms(3)+r)
    call rhs(2)%init(1,ms(1)*ms(3)+r)
    call mid(2)%init(r,ms(1)*ms(3)+r)
    call rhs(3)%init(1,ms(1)*ms(2)+r)
    call mid(3)%init(r,ms(1)*ms(2)+r)
  end if
  
  if (.not.present(masku)) then
    call rescale(u, v, w)
  end if
  do i = 1, ms(1)
    do j = 1, ms(2)
      do k = 1, ms(3)
        rhs(1)%d(1,ms(3)*(j-1)+k) = mmat(i,j,k)
        mid(1)%d(:,ms(3)*(j-1)+k) = v%d(:,j)*w%d(:,k)
      end do
    end do
    if (alpha .ne. 0) then
      do j = 1, r
        mid(1)%d(j,ms(2)*ms(3) + j) = alpha
      end do
    end if
    if (present(masku)) then
      do j = 1, r
        if (masku%d(j,i) == 1) then
          rhs(1)%d(1,1:ms(2)*ms(3)) = rhs(1)%d(1,1:ms(2)*ms(3)) - mid(1)%d(j,1:ms(2)*ms(3))*u%d(j,i)
          rhs(1)%d(1,ms(2)*ms(3)+j) = alpha*u%d(j,i)
          mid(1)%d(j,1:ms(2)*ms(3)) = 0
        end if
        if (maskzu%d(j,i) == 1) then
          rhs(1)%d(1,ms(2)*ms(3)+j) = 0
          mid(1)%d(j,1:ms(2)*ms(3)) = 0
        end if
      end do
    end if
    lhs(1) = rhs(1) .dI. mid(1)
    u%d(:,i) = lhs(1)%d(1,:)
  end do
  do i = 1, ms(2)
    do j = 1, ms(1)
      do k = 1, ms(3)
        rhs(2)%d(1,ms(3)*(j-1)+k) = mmat(j,i,k)
        mid(2)%d(:,ms(3)*(j-1)+k) = u%d(:,j)*w%d(:,k)
      end do
    end do
    if (alpha .ne. 0) then
      do j = 1, r
        mid(2)%d(j,ms(1)*ms(3) + j) = alpha
      end do
    end if
    if (present(masku)) then
      do j = 1, r
        if (maskv%d(j,i) == 1) then
          rhs(2)%d(1,1:ms(1)*ms(3)) = rhs(2)%d(1,1:ms(1)*ms(3)) - mid(2)%d(j,1:ms(1)*ms(3))*v%d(j,i)
          rhs(2)%d(1,ms(1)*ms(3)+j) = alpha*v%d(j,i)
          mid(2)%d(j,1:ms(1)*ms(3)) = 0
        end if
        if (maskzv%d(j,i) == 1) then
          rhs(2)%d(1,ms(1)*ms(3)+j) = 0
          mid(2)%d(j,1:ms(1)*ms(3)) = 0
        end if
      end do
    end if
    lhs(2) = rhs(2) .dI. mid(2)
    v%d(:,i) = lhs(2)%d(1,:)
  end do
  do i = 1, ms(3)
    do j = 1, ms(1)
      do k = 1, ms(2)
        rhs(3)%d(1,ms(2)*(j-1)+k) = mmat(j,k,i)
        mid(3)%d(:,ms(2)*(j-1)+k) = u%d(:,j)*v%d(:,k)
      end do
    end do
    if (alpha .ne. 0) then
      do j = 1, r
        mid(3)%d(j,ms(1)*ms(2) + j) = alpha
      end do
    end if
    if (present(masku)) then
      do j = 1, r
        if (maskw%d(j,i) == 1) then
          rhs(3)%d(1,1:ms(1)*ms(2)) = rhs(3)%d(1,1:ms(1)*ms(2)) - mid(3)%d(j,1:ms(1)*ms(2))*w%d(j,i)
          rhs(3)%d(1,ms(1)*ms(2)+j) = alpha*w%d(j,i)
          mid(3)%d(j,1:ms(1)*ms(2)) = 0
        end if
        if (maskzw%d(j,i) == 1) then
          rhs(3)%d(1,ms(1)*ms(2)+j) = 0
          mid(3)%d(j,1:ms(1)*ms(2)) = 0
        end if
      end do
    end if
    lhs(3) = rhs(3) .dI. mid(3)
    w%d(:,i) = lhs(3)%d(1,:)
  end do

end

!20% faster version, but only half precision (till 10^-7)
!Rename to ALS (and rename ALS to ALS2) to use this version
subroutine ALS2(u, v, w, ms, mmat, alpha, masku, maskv, maskw, maskzu, maskzv, maskzw)
  USE ModMtrx
  Type(Mtrx) :: u, v, w
  Integer(4), intent(in) :: ms(3)
  Double precision, dimension(ms(1),ms(2),ms(3)), intent(in) :: mmat
  Double precision, intent(in) :: alpha
  Type(Mtrx), intent(in), optional :: masku, maskv, maskw, maskzu, maskzv, maskzw
  
  Integer(4) i, j, k, r
  Type(Mtrx) mid(3)
  Type(Vector) rhs(3)
  Double precision, allocatable :: lwork(:)
  Integer(4) info
  Integer(4) lenwork
  
  r = u%n
  
  lenwork = 2*r*(max(ms(1),ms(2),ms(3))**2+r)
  Allocate(lwork(lenwork))
  
  if ((alpha == 0) .and. (.not.present(masku))) then
    call rhs(1)%init(ms(2)*ms(3))
    call mid(1)%init(r,ms(2)*ms(3))
    call rhs(2)%init(ms(1)*ms(3))
    call mid(2)%init(r,ms(1)*ms(3))
    call rhs(3)%init(ms(1)*ms(2))
    call mid(3)%init(r,ms(1)*ms(2))
  else
    call rhs(1)%init(ms(2)*ms(3)+r)
    call mid(1)%init(r,ms(2)*ms(3)+r)
    call rhs(2)%init(ms(1)*ms(3)+r)
    call mid(2)%init(r,ms(1)*ms(3)+r)
    call rhs(3)%init(ms(1)*ms(2)+r)
    call mid(3)%init(r,ms(1)*ms(2)+r)
  end if
  
  if (.not.present(masku)) then
    call rescale(u, v, w)
  end if
  do i = 1, ms(1)
    do j = 1, ms(2)
      do k = 1, ms(3)
        rhs(1)%d(ms(3)*(j-1)+k) = mmat(i,j,k)
        mid(1)%d(:,ms(3)*(j-1)+k) = v%d(:,j)*w%d(:,k)
      end do
    end do
    if (alpha .ne. 0) then
      do j = 1, r
        mid(1)%d(j,ms(2)*ms(3) + j) = alpha
      end do
    end if
    if (present(masku)) then
      do j = 1, r
        if (masku%d(j,i) == 1) then
          rhs(1)%d(1:ms(2)*ms(3)) = rhs(1)%d(1:ms(2)*ms(3)) - mid(1)%d(j,1:ms(2)*ms(3))*u%d(j,i)
          rhs(1)%d(ms(2)*ms(3)+j) = alpha*u%d(j,i)
          mid(1)%d(j,1:ms(2)*ms(3)) = 0
        end if
        if (maskzu%d(j,i) == 1) then
          rhs(1)%d(ms(2)*ms(3)+j) = 0
          mid(1)%d(j,1:ms(2)*ms(3)) = 0
        end if
      end do
    end if
    call dgels('T', mid(1)%n, mid(1)%m, 1, mid(1)%d, mid(1)%n, rhs(1)%d, rhs(1)%n, lwork, lenwork, info)
    u%d(:,i) = rhs(1)%d(1:r)
    rhs(1)%d(ms(2)*ms(3)+1:) = 0
    mid(1)%d(:,ms(2)*ms(3)+1:) = 0
  end do
  do i = 1, ms(2)
    do j = 1, ms(1)
      do k = 1, ms(3)
        rhs(2)%d(ms(3)*(j-1)+k) = mmat(j,i,k)
        mid(2)%d(:,ms(3)*(j-1)+k) = u%d(:,j)*w%d(:,k)
      end do
    end do
    if (alpha .ne. 0) then
      do j = 1, r
        mid(2)%d(j,ms(1)*ms(3) + j) = alpha
      end do
    end if
    if (present(masku)) then
      do j = 1, r
        if (maskv%d(j,i) == 1) then
          rhs(2)%d(1:ms(1)*ms(3)) = rhs(2)%d(1:ms(1)*ms(3)) - mid(2)%d(j,1:ms(1)*ms(3))*v%d(j,i)
          rhs(2)%d(ms(1)*ms(3)+j) = alpha*v%d(j,i)
          mid(2)%d(j,1:ms(1)*ms(3)) = 0
        end if
        if (maskzv%d(j,i) == 1) then
          rhs(2)%d(ms(1)*ms(3)+j) = 0
          mid(2)%d(j,1:ms(1)*ms(3)) = 0
        end if
      end do
    end if
    call dgels('T', mid(2)%n, mid(2)%m, 1, mid(2)%d, mid(2)%n, rhs(2)%d, rhs(2)%n, lwork, lenwork, info)
    v%d(:,i) = rhs(2)%d(1:r)
    rhs(2)%d(ms(1)*ms(3)+1:) = 0
    mid(2)%d(:,ms(1)*ms(3)+1:) = 0
  end do
  do i = 1, ms(3)
    do j = 1, ms(1)
      do k = 1, ms(2)
        rhs(3)%d(ms(2)*(j-1)+k) = mmat(j,k,i)
        mid(3)%d(:,ms(2)*(j-1)+k) = u%d(:,j)*v%d(:,k)
      end do
    end do
    if (alpha .ne. 0) then
      do j = 1, r
        mid(3)%d(j,ms(1)*ms(2) + j) = alpha
      end do
    end if
    if (present(masku)) then
      do j = 1, r
        if (maskw%d(j,i) == 1) then
          rhs(3)%d(1:ms(1)*ms(2)) = rhs(3)%d(1:ms(1)*ms(2)) - mid(3)%d(j,1:ms(1)*ms(2))*w%d(j,i)
          rhs(3)%d(ms(1)*ms(2)+j) = alpha*w%d(j,i)
          mid(3)%d(j,1:ms(1)*ms(2)) = 0
        end if
        if (maskzw%d(j,i) == 1) then
          rhs(3)%d(ms(1)*ms(2)+j) = 0
          mid(3)%d(j,1:ms(1)*ms(2)) = 0
        end if
      end do
    end if
    call dgels('T', mid(3)%n, mid(3)%m, 1, mid(3)%d, mid(3)%n, rhs(3)%d, rhs(3)%n, lwork, lenwork, info)
    w%d(:,i) = rhs(3)%d(1:r)
    rhs(3)%d(ms(1)*ms(2)+1:) = 0
    mid(3)%d(:,ms(1)*ms(2)+1:) = 0
  end do

  Deallocate(lwork)
  
end

!Integer to string converter
recursive function itoa(i,j) result(res)
  character(:),allocatable :: res
  integer,intent(in) :: i
  integer,intent(in),optional :: j
  character(range(i)+2) :: tmp
  if (present(j)) then
    write(tmp,'(i0.'//itoa(j)//')') i
  else
    write(tmp,'(i0)') i
  end if
  res = trim(tmp)
end function
