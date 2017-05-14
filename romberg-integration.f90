program main
  implicit none
  real*8::a,b,h,macheps,sum
  real*8,allocatable::r(:,:)
  integer::s,k,n,m
  macheps=1e-6
  s=100
  a=8.d0
  b=30.d0
  allocate(r(0:s,0:s))

  r(0,0)=((b-a)/2)*(f(a)+f(b))   !set first Romberg entry
  n=1
  do
     do m=1,n
        h=(b-a)/(2**n)
        sum=0.d0
        do k=1,2**(n-1)
           sum=sum+f(a+(2*k-1)*h)
        enddo
        r(n,0)=0.5d0*r(n-1,0)+(h*sum)
        r(n,m)=(((4.d0**m)*r(n,m-1))-r(n-1,m-1))/((4.d0**m)-1)
     enddo
     if (abs(r(n,n)-r(n-1,n-1))<macheps)then         !Checks stopping condition (solution converges)
        write(6,'(a26,i3)')"Romberg n (index from 0):",n
        write(6,'(a3,f14.6)')"I:",r(n,n)
        write(6,'(a11,es12.2e2)')"Precision:",macheps
        exit
     endif
     n=n+1
  enddo


  deallocate(r)

contains

  real(kind=8) function f(z)    !Place any single variable function to be evaluated here
    implicit none                  !can remove om 
    real(kind=8),intent(in)::z
    f=(2000.d0*log(140000.d0/(140000.d0-(2100.d0*z))))-9.8d0*z
  end function f

end program main
