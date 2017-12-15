module cut_off_function
implicit none

contains

function f_cut(r,R1,R2)
	real f_cut,r,R1,R2,pi
	pi=3.14159265358979
	if(r<R1) then
		f_cut = 1.
	elseif(r<R2) then
		f_cut = (1.+cos( pi*(r-R1)/(R2-R1) ))/2
	else
		f_cut = 0.
	endif
end function f_cut

function df_cut(r,R1,R2)
	real df_cut,r,R1,R2,pi
	pi=3.14159265358979
	if(r<R1) then
		df_cut = 0.
	elseif(r<R2) then
		df_cut = -sin( pi*(r-R1)/(R2-R1) )*pi/(R2-R1)/r/2
	else
		df_cut = 0.
	endif	
end function df_cut

end module cut_off_function