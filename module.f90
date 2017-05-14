module main


CONTAINS


  subroutine cspline(n,Y,ans,h)
    implicit none
    integer,intent(in)::n
    integer::i
    real(kind=8)::a(n-2),b(n-2),c(n-2),d(n-2),ans(n-2),anot,bnot,cnot,dnot
    real(kind=8)::first(n-2),second(n-2),third(n-2),fourth(n-2)
    real(kind=8),intent(in)::y(n),h

    open(19,file="ddata.dat")
    do i=1,n-2
       write(19,*)(y(i)-(2.d0*y(i+1))+y(i+2))*6.do/(h**2.d0)
    enddo
    close(19)


    open(13,file='ddata.dat')
    do i=1,n-2
       a(i)=1.d0
       b(i)=4.d0
       c(i)=1.d0
       read(13,*)d(i)
    enddo
    close(13)
    call solve_tridiag(a,b,c,d,ans,n-2)
    write(6,*)ans
    anot=ans(1)/(6.d0*h)
    bnot=0.d0
    cnot=((y(2)-y(1))/h)-(ans(1)*h/6.d0)
    dnot=y(1)
    write(6,'(A3,F14.10)')"a1=",anot
    write(6,'(A3,F14.10)')"b1=",bnot
    write(6,'(A3,F14.10)')"c1=",cnot
    write(6,'(A3,F14.10)')"d1=",dnot
    write(6,*)

    do i=2,n-2
       first(i)=(ans(i)-ans(i-1))/(6.d0*h)
       second(i)=ans(i-1)/2.d0
       third(i)=(((Y(i+1)*1.d0)-(1.d0*Y(i)))/h)-((ans(i)+(2.d0*ans(i-1)))*h/6.d0)
       fourth(i)=Y(i)*1.d0
       write(6,'(1X,A1,I1,A1,F14.10)')"a",i,"=",first(i)
       write(6,'(1X,A1,I1,A1,F14.10)')"b",i,"=",second(i)
       write(6,'(1X,A1,I1,A1,F14.10)')"c",i,"=",third(i)
       write(6,'(1X,A1,I1,A1,F14.10)')"d",i,"=",fourth(i)
       write(6,*)
    enddo
    anot=(-1.d0*ans(n-2))/(6.d0*h)
    bnot=ans(n-2)/2.d0
    cnot=((y(n)-y(n-1))/h)-(2.d0*ans(n-2)*h/6.d0)
    dnot=y(n-1)
    write(6,'(A3,F14.10)')"a8=",anot
    write(6,'(A3,F14.10)')"b8=",bnot
    write(6,'(A3,F14.10)')"c8=",cnot
    write(6,'(A3,F14.10)')"d8=",dnot

  end subroutine cspline





  subroutine csplinecub(n,Y,ans,h)
    implicit none
    integer,intent(in)::n
    integer::i
    real(kind=8)::a(n-2),b(n-2),c(n-2),d(n-2),ans(n-2),anot,bnot,cnot,dnot
    real(kind=8)::first(n-2),second(n-2),third(n-2),fourth(n-2)
    real(kind=8),intent(in)::y(n),h

    open(19,file="ddata.dat")
    do i=1,n-2
       write(19,*)(y(i)-(2.d0*y(i+1))+y(i+2))*6.do/(h**2.d0)
    enddo
    close(19)

    open(13,file='ddata.dat')
    do i=1,n-2
       a(i)=1.d0
       b(i)=4.d0
       c(i)=1.d0
       read(13,*)d(i)
    enddo
    b(1)=6.d0
    b(n-2)=6.d0
    c(1)=0.d0
    a(n-3)=0.d0
    close(13)
    call solve_tridiag(a,b,c,d,ans,n-2)
    write(6,*)ans
    anot=(-ans(1)-ans(2))/(6.d0*h)
    bnot=(2.d0*ans(1)-ans(2))/2.d0
    cnot=((y(2)-y(1))/h)-((5.d0*ans(1)-2.d0*ans(2))*h/6.d0)
    dnot=y(1)
    write(6,'(A3,F14.10)')"a1=",anot
    write(6,'(A3,F14.10)')"b1=",bnot
    write(6,'(A3,F14.10)')"c1=",cnot
    write(6,'(A3,F14.10)')"d1=",dnot
    write(6,*)

    do i=2,n-2
       first(i)=(ans(i)-ans(i-1))/(6.d0*h)
       second(i)=ans(i-1)/2.d0
       third(i)=(((Y(i+1)*1.d0)-(1.d0*Y(i)))/h)-((ans(i)+(2.d0*ans(i-1)))*h/6.d0)
       fourth(i)=Y(i)*1.d0
       write(6,'(1X,A1,I1,A1,F14.10)')"a",i,"=",first(i)
       write(6,'(1X,A1,I1,A1,F14.10)')"b",i,"=",second(i)
       write(6,'(1X,A1,I1,A1,F14.10)')"c",i,"=",third(i)
       write(6,'(1X,A1,I1,A1,F14.10)')"d",i,"=",fourth(i)
       write(6,*)
    enddo
    anot=(ans(n-2)-ans(n-3))/(6.d0*h)
    bnot=ans(n-2)/2.d0
    cnot=((y(n)-y(n-1))/h)-((4.d0*ans(n-2)-ans(n-3))*h/6.d0)
    dnot=y(n-1)
    write(6,'(A3,F14.10)')"a8=",anot
    write(6,'(A3,F14.10)')"b8=",bnot
    write(6,'(A3,F14.10)')"c8=",cnot
    write(6,'(A3,F14.10)')"d8=",dnot

  end subroutine csplinecub






  subroutine csplinepara(n,Y,ans,h)
    implicit none
    integer,intent(in)::n
    integer::i
    real(kind=8)::a(n-2),b(n-2),c(n-2),d(n-2),ans(n-2),anot,bnot,cnot,dnot
    real(kind=8)::first(n-2),second(n-2),third(n-2),fourth(n-2)
    real(kind=8),intent(in)::y(n),h

    open(19,file="ddata.dat")
    do i=1,n-2
       write(19,*)(y(i)-(2.d0*y(i+1))+y(i+2))*6.do/(h**2.d0)
    enddo
    close(19)

    open(13,file='ddata.dat')
    do i=1,n-2
       a(i)=1.d0
       b(i)=4.d0
       c(i)=1.d0
       read(13,*)d(i)
    enddo
    b(1)=5.d0
    b(n-2)=5.d0
    close(13)
    call solve_tridiag(a,b,c,d,ans,n-2)
    write(6,*)ans
    anot=0.d0
    bnot=ans(1)/2.d0
    cnot=(((y(2)*1.d0)-(1.d0*y(1)))/h)-((3.d0*ans(1))*h/6.d0)
    dnot=y(1)
    write(6,'(A3,F14.10)')"a1=",anot
    write(6,'(A3,F14.10)')"b1=",bnot
    write(6,'(A3,F14.10)')"c1=",cnot
    write(6,'(A3,F14.10)')"d1=",dnot
    write(6,*)

    do i=2,n-2
       first(i)=(ans(i)-ans(i-1))/(6.d0*h)
       second(i)=ans(i-1)/2.d0
       third(i)=(((Y(i+1)*1.d0)-(1.d0*Y(i)))/h)-((ans(i)+(2.d0*ans(i-1)))*h/6.d0)
       fourth(i)=Y(i)*1.d0
       write(6,'(1X,A1,I1,A1,F14.10)')"a",i,"=",first(i)
       write(6,'(1X,A1,I1,A1,F14.10)')"b",i,"=",second(i)
       write(6,'(1X,A1,I1,A1,F14.10)')"c",i,"=",third(i)
       write(6,'(1X,A1,I1,A1,F14.10)')"d",i,"=",fourth(i)
       write(6,*)
    enddo
    anot=0.d0
    bnot=ans(n-2)/2.d0
    cnot=((y(n)-y(n-1))/h)-(3.d0*ans(n-2)*h/6.d0)
    dnot=y(n-1)
    write(6,'(A3,F14.10)')"a8=",anot
    write(6,'(A3,F14.10)')"b8=",bnot
    write(6,'(A3,F14.10)')"c8=",cnot
    write(6,'(A3,F14.10)')"d8=",dnot

  end subroutine csplinepara




  subroutine solve_tridiag(a,b,c,d,x,n)
    implicit none
    !        a - sub-diagonal (means it is the diagonal below the main diagonal)
    !        b - the main diagonal
    !        c - sup-diagonal (means it is the diagonal above the main diagonal)
    !        d - right part
    !        x - the answer
    !        n - number of equations

    integer,intent(in) :: n
    real(kind=8),dimension(n),intent(in) :: a,b,c,d
    real(kind=8),dimension(n),intent(out) :: x
    real(kind=8),dimension(n) :: cp,dp
    real(kind=8) :: m
    integer i

    ! initialize c-prime and d-prime
    cp(1) = c(1)/b(1)
    dp(1) = d(1)/b(1)
    ! solve for vectors c-prime and d-prime
    do i = 2,n
       m = b(i)-cp(i-1)*a(i)
       cp(i) = c(i)/m
       dp(i) = (d(i)-dp(i-1)*a(i))/m
    enddo
    ! initialize x
    x(n) = dp(n)
    ! solve for x from the vectors c-prime and d-prime
    do i = n-1, 1, -1
       x(i) = dp(i)-cp(i)*x(i+1)
    end do

  end subroutine solve_tridiag




  !*******************************************************************
  !*                                                                 *
  !*  Romberg takes in the (finite) integration limits and outputs   *
  !*  the value of the integral. The integral function is stored in  *
  !*  the paired Fortran function "f(z)" and must be included with   *
  !*  subroutine romberg.                                            *
  !*                                                                 *
  !*******************************************************************


  subroutine romberg(a,b)
    implicit none
    real*8,intent(in)::a,b
    real*8::macheps,sum,h
    real*8,allocatable::r(:,:)
    integer::s,k,n,m
    macheps=1e-6
    s=100
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
    
    
  end subroutine romberg
  
  
  
  subroutine monte(sample,x1,x2,y1,y2)  !Pass monte the desired sample size and x-y limits 
    implicit none                       !of a box region.
    real*8::x,y,hit,miss,integral       !Also need to define f(x) within f90 function "f(z)"
    real*8,intent(in)::x1,x2,y1,y2
    real*8::mean,int(sample),sum,stddev
    integer::i,n
    integer,intent(in)::sample
    integral=0.d0

    do n=1,sample
       hit=0.d0
       miss=0.d0
       do i=1,100
          x=(rand()*(x2-x1))+x1
          y=(rand()*(y2-y1))+y1
          if (y.le.f(x))then
             hit=hit+1.d0
          else
             miss=miss+1.d0
          endif
       enddo
       int(n)=((abs((x2-x1)*(y2-y1))*hit)/(hit+miss))
       integral= integral+((abs((x2-x1)*(y2-y1))*hit)/(hit+miss))
    enddo
    mean=integral/(1.d0*sample)
    sum=0.d0
    do i=1,sample
       sum=sum+((int(i)-mean)**2.d0)
    enddo
    stddev=sqrt(sum/(sample*1.d0))
    write(6,*)"Monte Carlo Hit and Miss"
    write(6,'(a13,i15)')'Sample size:',sample
    write(6,'(a33,f12.6)')"Mean of the samples on integral:",mean
    write(6,'(a30,f12.6)')'Standard deviation of sample:',stddev
  end subroutine monte





  real(kind=8) function f(z)
    implicit none
    real(kind=8),intent(in)::
    f=
  end function f




!*****PREDICTOR-CORRECTOR FOR SINGLE ODE*********!


  subroutine rk4(curve_constraint,resolution,xnot,ynot,max_x,filename)
    implicit none
    real*8::k1,k2,k3,k4,step,xnot,ynot,x,y,curve_constraint,resolution,slopecheck
    real*8::initialslope,stepslope,max_x
    integer::i,n
    character(len=30)::filename
    x=xnot
    y=ynot
    open(12,file=filename)
    write(12,*)x,y
    do while(x.le.max_x)
       step=resolution

       do
          k1=step*df(x,y)
          k2=step*df(x+(step/2.d0),y+(k1/2.d0))
          k3=step*df(x+(step/2.d0),y+(k2/2.d0))
          k4=step*df(x+step,y+k3)

          initialslope=df(x,y)                                      
          stepslope=df(x+step,y+((k1+(2.d0*k2)+(2.d0*k3)+k4)/6.d0)) 
          slopecheck=abs(stepslope-initialslope)/abs(initialslope)     

          if(slopecheck.gt.curve_constraint)then                    
             step=step/2.d0                                         
             cycle                                                  
          else                                                      
             y=y+((k1+(2.d0*k2)+(2.d0*k3)+k4)/6.d0)                 
             x=x+step                                               
             exit                                                   
          endif                                                     

       enddo
       write(12,*)x,y
    enddo
    close(12)
  end subroutine rk4



!***********RK4 FOR SYSTEM OF 3 ODEs*********************


  subroutine rk4(resolution,tnot,xnot,ynot,znot,max_t,filename,sig,r,b)
    implicit none
    real*8::xk1,xk2,xk3,xk4,yk1,yk2,yk3,yk4,zk1,zk2,zk3,zk4
    real*8::max_t,t,tnot,sig,r,b,step,xnot,ynot,znot,x,y,z,resolution
    integer::i,n
    character(len=30)::filename
    step=resolution
    x=xnot
    y=ynot
    z=znot
    t=tnot
    open(12,file=filename)
    write(12,*)t,x,y,z
    do while(t.le.max_t)


!       do
          xk1=step*dfx(x,y,z,sig)
          yk1=step*dfy(x,y,z,r)
          zk1=step*dfz(x,y,z,b)

          xk2=step*dfx(x+(xk1/2.d0),y+(yk1/2.d0),z+(zk1/2.d0),sig)
          yk2=step*dfy(x+(xk1/2.d0),y+(yk1/2.d0),z+(zk1/2.d0),r)
          zk2=step*dfz(x+(xk1/2.d0),y+(yk1/2.d0),z+(zk1/2.d0),b)

          xk3=step*dfx(x+(xk2/2.d0),y+(yk2/2.d0),z+(zk2/2.d0),sig)
          yk3=step*dfy(x+(xk2/2.d0),y+(yk2/2.d0),z+(zk2/2.d0),r)
          zk3=step*dfz(x+(xk2/2.d0),y+(yk2/2.d0),z+(zk2/2.d0),b)

          xk4=step*dfx(x+xk3,y+yk3,z+zk3,sig)
          yk4=step*dfy(x+xk3,y+yk3,z+zk3,r)
          zk4=step*dfz(x+xk3,y+yk3,z+zk3,b)

          !          initialslope=df(x,y)                                      !predictor/corrector
          !          stepslope=df(x+step,y+((k1+(2.d0*k2)+(2.d0*k3)+k4)/6.d0)) !implementation
          !          slopecheck=abs(stepslope-initialslope)/abs(initialslope)     

          !          if(slopecheck.gt.curve_constraint)then                    
          !             step=step/2.d0                                         
          !             cycle                                                  
          !          else
          x=x+((xk1+(2.d0*xk2)+(2.d0*xk3)+xk4)/6.d0)
          y=y+((yk1+(2.d0*yk2)+(2.d0*yk3)+yk4)/6.d0)
          z=z+((zk1+(2.d0*zk2)+(2.d0*zk3)+zk4)/6.d0)
          t=t+step                                               
          !             exit                                                   
          !          endif                                                     

!       enddo
       write(12,*)t,x,y,z
    enddo
    close(12)
  end subroutine rk4



!**********RK4 FOR SYSTEM OF 2 ODEs*************



  subroutine rk4(resolution,tnot,vnot,thnot,max_t,filename)
    implicit none
    real*8::vk1,vk2,vk3,vk4,thk1,thk2,thk3,thk4
    real*8::max_t,t,tnot,step,vnot,thnot,v,th,resolution
    integer::i,n
    character(len=30)::filename
    step=resolution
    v=vnot
    th=thnot
    t=tnot
    open(12,file=filename)
    write(12,*)t,th,v
    do while(t.le.max_t)


!       do
          thk1=step*df1(v)
          vk1=step*df2(v,th,t)


          thk2=step*df1(v+(vk1/2.d0))
          vk2=step*df2(v+(vk1/2.d0),th+(thk1/2.d0),t+(step/2.d0))


          thk3=step*df1(v+(vk1/2.d0))
          vk3=step*df2(v+(vk2/2.d0),th+(thk2/2.d0),t+(step/2.d0))


          thk4=step*df1(v+vk3)
          vk4=step*df2(v+vk3,th+thk3,t+step)


          !          initialslope=df(x,y)                                      !predictor/corrector
          !          stepslope=df(x+step,y+((k1+(2.d0*k2)+(2.d0*k3)+k4)/6.d0)) !implementation
          !          slopecheck=abs(stepslope-initialslope)/abs(initialslope)     

          !          if(slopecheck.gt.curve_constraint)then                    
          !             step=step/2.d0                                         
          !             cycle                                                  
          !          else
          v=v+((vk1+(2.d0*vk2)+(2.d0*vk3)+vk4)/6.d0)
          th=th+((thk1+(2.d0*thk2)+(2.d0*thk3)+thk4)/6.d0)

          t=t+step                                               
          !             exit                                                   
          !          endif                                                     

!       enddo
       write(12,*)t,th,v
    enddo
    close(12)
  end subroutine rk4



!**********LAGRANGE FUNCTION****************
!xx = value you want the interpolant at
!x and y are the data points of n size


  real(kind=8) function lagrange(x,y,n,xx)
    implicit none
    integer::i,k
    integer,intent(in)::n
    real(kind=8)::product,sum
    real(kind=8),intent(in)::x(n),y(n),xx
    sum=0.d0
    do i=1,n
       product=y(i)
       do k=1,n
          if(i.ne.k)then
             product=product*(xx-x(k))/(x(i)-x(k))
          endif
       enddo
       sum=sum+product
    enddo
    lagrange=sum

  end function lagrange





  
  !***NEWTON-RAPHSON METHOD --- Include the functions!!!***
  i=2
  x(1)=20.d0
  delta=1.d0
  open(12,file='raphson2.dat')
  do 

     f=f(x(i-1))
     fprime=fprime(x(i-1))
     
     c=f/fprime
     x(i)=x(i-1)-c
     delta=abs(x(i)-x(i-1))
     
     if(delta.le.accuracy)then
        write(12,*)i-1,x(i),0.d0,delta
        exit
     endif
     write(12,*)i-1,x(i),0.d0,delta
     i=i+1
     if(i.gt.96)then
        write(6,*)"Failure to converge to solution."
        stop
     endif
  enddo
close(12)

  
  !***SECANT METHOD --- Include the functions!!!***
  delta=1.d0
  i=3
  x(1)=20.d0
  x(2)=30.d0
  open(13,file='secant2.dat')
  do
     f=f(x(i-1))
     s=(f(x(i-1))-f(x(i-2)))/(x(i-1)-x(i-2))
          
     c=f/s
     x(i)=x(i-1)-c
     delta=abs(x(i)-x(i-1))
     if(delta.le.accuracy)then
        write(13,*)i-2,x(i),0.d0,delta
        exit
     endif
     write(13,*)i-2,x(i),0.d0,delta
     i=i+1
     if(i.gt.96)then
        write(6,*)"Failure to converge to solution."
        stop
     endif
  enddo




end module main
