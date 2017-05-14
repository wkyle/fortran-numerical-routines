program main
  implicit none
  integer,parameter::nx=50, ny=50
  integer::i,j,it,maxit
  real*4,parameter::pi=acos(-1.d0), dx=1.d0/50.d0, dy=dx, &
       uleft=0.d0, utop=0.d0, ubottom=0.d0, uright=0.d0
  character(len=64),parameter::file1="gs1.dat"
  real*4::u(0:nx+1,0:ny+1),store(0:nx+1,0:ny+1),maxr,tempr

maxit=10000

  u=0.d0
  do i=0,ny
     u(0,i)=uleft
  enddo
  do i=0,ny
     u(nx+1,i)=uright
  enddo
  do i=0,nx
     u(i,0)=ubottom
  enddo
  do i=0,nx
     u(i,ny+1)=utop
  enddo
  store=1.d20

  do it=1,maxit
     maxr=-1.d0
     do i=1,nx
        do j=1,ny
           u(i,j)=.25d0*(u(i-1,j)+u(i+1,j)+u(i,j-1)+u(i,j+1))-(10.d0*dx*dx/4.d0)

           tempr=abs(u(i,j)-store(i,j))
           if(tempr.gt.maxr)then
              maxr=tempr
write(6,*)"changed maxr to",maxr
           endif
        enddo
     enddo
     if(maxr.lt.1.d-3)then
        write(6,*)"# iterations=",it
        exit
     else
        store=u
     endif
  enddo


  open(12,file=file1)
  do j=1,ny
     write(12,*)(u(i,j), i=1,nx)
  enddo
  close(12)

end program main
