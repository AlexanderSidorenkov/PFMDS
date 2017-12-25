module cut_off_poly
implicit none

contains

function f_cut(r,R1,R2)
	real f_cut
	real :: r,R1,R2
	
	if(r>R1 .and. r<R2) then
		f_cut = 1.+(-10.*(r-R1)**3*(R2-R1)**2+15.*(r-R1)**4*(R2-R1)-6.*(r-R1)**5)/(R2-R1)**5
	elseif(r<R1) then
		f_cut = 1.
	else
		f_cut = 0.
	endif
end function f_cut

function df_cut(r,R1,R2)
	real df_cut
	real :: r,R1,R2

	if(r>R1 .and. r<R2) then
		df_cut = (-30.*(r-R1)**2*(R2-R1)**2+60.*(r-R1)**3*(R2-R1)-30.*(r-R1)**4)/(R2-R1)**5/r
	else
		df_cut = 0.
	endif	
end function df_cut

subroutine f_dfr_cut(f,dfr,r,R1,R2)
	real :: f,dfr,r,R1,R2,x,x2!,s,ss

	!x = (r-R1)/(R2-R1)
	!x2 = x*x
	!s = ((R2-r)/abs(R2-r)+1)/2
	!ss = ((r-R1)/abs(r-R1)*(R2-r)/abs(R2-r)+1)/2
	
	!if (s>1.1 .or. ss>1.1) print*,x,s,ss
	!if (s/=s) print*,x,s,ss
	!if (ss/=ss) print*,x,s,ss
	
	!f = s*(1.+x2*x*(-10.+15.*x-6.*x2))
	!dfr = ss*(x2*(-30.+60.*x-30.*x2)*r)
	
	if(r>R1 .and. r<R2) then
		x = (r-R1)/(R2-R1)
		x2 = x*x
		f = 1.+x2*x*(-10.+15.*x-6.*x2)
		dfr = x2*(-30.+60.*x-30.*x2)*r
	elseif(r<R1) then
		f = 1.
		dfr = 0.
	else
		f = 0.
		dfr = 0.
	endif	
	
end subroutine f_dfr_cut

end module cut_off_poly