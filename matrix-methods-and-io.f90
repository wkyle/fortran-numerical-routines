!This is all code snippets for the final exam
!including subroutines, functions, and formatting
!
!*****************************************************

!*******************
!FACTORIAL
!*******************

CONTAINS
  
  real*8 function factorial(n)
    integer,intent(in)::n
    real*8::sofar
    integer::k
    
    if (n .eq. 0 .or. n .eq. 1) then
       factorial=1
    else
       sofar=1.
       do k=n,1,-1
          sofar=sofar*k
       enddo
       factorial=anint(sofar,8)
    endif
    
  end function factorial

!*******************************************
!GAUSSIAN-JORDAN ELIMINATION MATRIX SOLVER
!*******************************************

  subroutine naivegauss(A,n,ST)
    implicit none
    integer::i,j,k
    integer,intent(in)::n
    real*8,intent(in)::A(n+1,n),ST(n+1,n)
    real*8::soln(n),sum
    real*8::O(n+1,n),P(n+1,n)
    O=A
    P=ST
    
!****************************************                                           
!Initializes top row of storage matrix                                                
!****************************************                                             
    do i=1,n+1
       P(i,1)=O(i,1)
    enddo
    
!**********************************************************************************
!This is the loop that reduces to row echelon form, or an upper triangular matrix
!**********************************************************************************
    do k=1,n-1
       do j=k+1,n
          do i=k,n+1
             P(i,j)=O(i,j)-(O(i,k)/O(k,k)*O(k,j))
          enddo
       enddo
       O=P
    enddo

!************************************************************************
!This section performs back substitution and writes the solution
!to the screen/file in proper row-column form
!************************************************************************
    soln(n)=O(n+1,n)/O(n,n)
    write(6,*)'The naive Gauss solution is...'
    write(6,'(A1,I1,A1,F20.16)')'x',n,'=',soln(n)
    do j=n-1,1,-1
       sum=0
       do i=j+1,n
          sum=sum+(O(i,j)*soln(i))
       enddo
       soln(j)=(O(n+1,j)-sum)/O(j,j)
       write(6,'(A1,I1,A1,F20.16)')'x',j,'=',soln(j)
    enddo
    
  end subroutine naivegauss
  
!********************************************************
!LU DECOMPOSITION MATRIX SOLVER
!********************************************************


  subroutine lud(A,n,ST)
    implicit none
    integer,intent(in)::n
    integer::i,j,k
    real*8::sum,U(n,n),L(n,n),y(n),soln(n),O(n+1,n),P(n+1,n)
    real*8,intent(in)::A(n+1,n),ST(n+1,n)
    
    O=A
    P=ST
    
!****************************************                                           
!Initializes top row of storage matrix                                                
!****************************************                                             
    do i=1,n+1
       P(i,1)=O(i,1)
    enddo
    
!**********************************************************************************
!This is the loop that reduces square matrix to a lower triangular matrix
!**********************************************************************************
    L=0
    do i=1,n
       L(i,i)=1
    enddo
    do k=1,n-1
       do j=k+1,n
          do i=k,n+1
             P(i,j)=O(i,j)-(O(i,k)/O(k,k)*O(k,j))
          enddo
          L(k,j)=O(k,j)/O(k,k)
       enddo
       O=P
    enddo
    do j=1,n
       do i=1,n
          U(i,j)=O(i,j)
       enddo
    enddo

!***********************************************************************
!Forward substitution (only for square matrix with normalized diagonal)
!***********************************************************************
    y(1)=O(n+1,1)
    
    do j=2,n
       sum=0
       do i=1,j-1
          sum=sum+L(i,j)*y(i)
       enddo
       y(j)=O(n+1,j)-sum
    enddo
    
!****************************************************************
!Backsubstitution 
!****************************************************************
    soln(n)=O(n+1,n)/U(n,n)
    write(6,*)'The LU decomposed solution is...'
    write(6,'(A1,I1,A1,F20.16)')'x',n,'=',soln(n)
    do j=n-1,1,-1
       sum=0
       do i=j+1,n
          sum=sum+(U(i,j)*soln(i))
       enddo
       soln(j)=(O(n+1,j)-sum)/U(j,j)
       write(6,'(A1,I1,A1,F20.16)')'x',j,'=',soln(j)
    enddo
    
    
    
  end subroutine lud

!********************************************
!ALLOCATABLE ARRAYS 
!********************************************

  real(kind=8), allocatable::A(:,:),ST(:,:),soln(:),L(:,:),U(:,:),B(:,:)
  integer::n
  write(6,*)"How big?"
  read(5,*)n
  allocate (A(n+1,n),ST(n+1,n),soln(n),L(n,n),U(n,n),B(n,n))

!*****************************************
!RANDOM NUMBER and COMPLEX NUMBERS
!*****************************************

integer :: i
real :: x,y
complex :: z
call random_seed()
call random_number(x)
z=cmplx(x,y)
enddo

!******************************************
!FORMATTING OUTPUT EXAMPLE
!******************************************

write(6,'(A5,4ES10.3E2,2x,2F9.4)')x,y,int

Exponential:  rESw.d and rESw.dEe
Real:         rFw.d
Character:    rAw
Integer:      rIw
Space:        rX

!*************************************
!MODULE EXAMPLE
!*************************************

When compiling main: gfortran module1.f90 module2.f90 main.f90 -o main.exe

In main: USE module_name_of_program

module summer
implicit none
real::num1,num2,sum

CONTAINS

  subroutine adder (sum,num1,num2)
    real,intent(in)::num1,num2
    real,intent(out)::sum
    
    sum=num1+num2
    
  end subroutine adder

  real function adding (num1,num2)
    implicit none
    real::adding
    real,intent(in)::num1,num2
    adding=num1+num2
  end function adding

end module summer

!*****************************************
!CALCULATING PI
!*****************************************
pi=acos(-1.d0)
