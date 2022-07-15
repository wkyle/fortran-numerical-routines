program main
use function
   implicit none
   real*8::a,b,macheps,suma_i
   real*8,allocatable::r(:,:)
   integer::s,k,n,m,i,j
   macheps=1e-6
write(*,*)'Enter the interval of integration from highest to lowest'
   read(*,*)a,b
write(*,*)'Enter J such that 2**(J-1) are the subintervals'
read(*,*)n
allocate(R(n,n))
R(1,1)=(b-a)*0.5*(f(a)+f(b))
write(*,*)'i=',1,"j=",1,'R(i,j)=',R(1,1)
   do i=2,n 
      suma_i=0
      do k=1,(2**(i-2)),1
      suma_i=suma_i+f(a+((b-a)/2**(i-2))*(k -0.5))
      end do
      R(i,1)=0.5*(R(i-1,1)+(((b-a)/(2**(i-2)))*suma_i))
      write(*,*)'i=',i,"j=",1,'R(i,j)=',R(i,1)
   end do
   do j=2,n
      do i=j,n
         R(i,j)=((4**(j-1))*R(i,j-1)-R(i-1,j-1))/(4**(j-1)-1)
         write(*,*)"i=",i,'j=',j,"R(i,j)=",R(i,j)
      end do
      if (abs(r(n,n)-r(n-1,n-1))<macheps)then         !Checks stopping condition (solution converges)
         write(6,'(a26,i3)')"Romberg n (index from 0):",n
         write(6,'(a3,f14.6)')"I:",R(n,n)
         write(6,'(a11,es12.2e2)')"Precision:",macheps
         
      endif
   end do
   deallocate(R)
end program main


