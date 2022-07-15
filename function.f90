
module function
contains 
function f(x)   
real*8,intent(in)::x
f=(2000.d0*log(140000.d0/(140000.d0-(2100.d0*x))))-9.8d0*x
return
end function f
end module 