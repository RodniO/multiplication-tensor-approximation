
program bubr
  USE ModMtrx
  Integer(4), parameter, dimension(3) :: mms = [3,3,3] !Matrix sizes for multiplication
  Integer(4), parameter :: r = 23 !Desired rank
  Integer(4), parameter :: tries = 1 !Total number of tries
  Double precision, dimension(mms(1)*mms(2),mms(2)*mms(3),mms(3)*mms(1)) :: mmat, emat !Multiplication tensor and error tensor
  Integer(4) i, j, k, l !Indices
  Type(Mtrx) :: u, v, w !Low-rank factors. W is transposed for symmetry reasons
  Type(Mtrx) mid !Temporary matrix
  Type(Mtrx) :: mats(r*3) !Matrices for rows of u, v and w
  Integer(4) :: ms(3) !Number of rows in u, v and w
  
  ms = [mms(1)*mms(2), mms(2)*mms(3), mms(3)*mms(1)]
  
  k = 0
  l = 0
  do j = 1, tries
    i = ExampleTEN(2, j, mms, r, 100, 10000, 100000, 20.0d0, 0.05d0, 0.1d0**11, 0.6d0)
    if (i == 2) then
      k = k + 1
    end if
    if (i == 1) then
      l = l + 1
    end if
  end do
  print *, 'SUCCESSES', k
  print *, 'NEAR SUCCESSES', l
  !REMOVE STOP 0 TO MAKE TENSOR OUTPUT IN TEXT FORMAT
  stop 0
  
  mmat(:,:,:) = 0
  do i = 1, mms(1)
    do j = 1, mms(2)
      do k = 1, mms(3)
        mmat(mms(2)*(i-1)+j,mms(3)*(j-1)+k,mms(1)*(k-1)+i) = 1
      end do
    end do
  end do
  
  call u%init(r,ms(1))
  call v%init(r,ms(2))
  call w%init(r,ms(3))
  
  print *, 'DATA FROM LAST SUCCESS:'
  
  !INSERT CORRECT FILE NUMBER INSTEAD OF ?? HERE
  open(unit=1, file='out??.r4', form='unformatted', access='direct', recl=r*max(ms(1),ms(2),ms(3))*8)
  
  read(1, rec=1) u%d
  read(1, rec=2) v%d
  read(1, rec=3) w%d
  
  close(1)
  
  !1/2 is replaced with 2
  do i = 1, r
    do j = 1, ms(1)
      if (nint(2*abs(u%d(i,j))) == 1) then
        u%d(i,j) = sign(2.0d0,u%d(i,j))
      end if
    end do
    do j = 1, ms(2)
      if (nint(2*abs(v%d(i,j))) == 1) then
        v%d(i,j) = sign(2.0d0,v%d(i,j))
      end if
    end do
    do j = 1, ms(3)
      if (nint(2*abs(w%d(i,j))) == 1) then
        w%d(i,j) = sign(2.0d0,w%d(i,j))
      end if
    end do
  end do
  open(2, file='textout.txt')
  do i = 1, r
    write(2, *) nint(u%d(i,:)), nint(v%d(i,:)), nint(w%d(i,:))
  end do
  close(2)
  
  do i = 1, r
    do j = 1, ms(1)
      u%d(i,j) = nint(u%d(i,j)*2)/2.0d0
    end do
    do j = 1, ms(2)
      v%d(i,j) = nint(v%d(i,j)*2)/2.0d0
    end do
    do j = 1, ms(3)
      w%d(i,j) = nint(w%d(i,j)*2)/2.0d0
    end do
  end do
  
  print *, 'NONZEROS (IF ONLY 0,1,-1)', nint(u%fnorm()**2+v%fnorm()**2+w%fnorm()**2)
  
  print *, ''
  do i = 1, r*3
    call mats(i)%init(mms((i-1)/r+1), mms(mod((i-1)/r+1,3)+1))
  end do
  print *, 'TRACES:'
  do i = 1, r
    do j = 1, mms(1)
      do k = 1, mms(2)
        mats(i)%d(j,k) = u%d(i,mms(2)*(j-1)+k)
      end do
    end do
  end do
  do i = 1, r
    do j = 1, mms(2)
      do k = 1, mms(3)
        mats(r+i)%d(j,k) = v%d(i,mms(3)*(j-1)+k)
      end do
    end do
  end do
  do i = 1, r
    do j = 1, mms(3)
      do k = 1, mms(1)
        mats(2*r+i)%d(j,k) = w%d(i,mms(1)*(j-1)+k)
      end do
    end do
  end do
  do i = 1, r
    print *, 'MATRIX', i
    mid = mats(i)*mats(r+i)*mats(2*r+i)
    print *, mid%trace()
  end do
    
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
  print *, 'ERROR', sqrt(sum(emat))
  
  contains

include "incfiles/ExampleTEN.f90"

end
