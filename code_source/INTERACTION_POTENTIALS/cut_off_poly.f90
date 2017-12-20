module cut_off_poly
implicit none

contains

function f_cut(r,R1,R2)
	real f_cut,r,R1,R2

	if(r>R1 .and. r<R2) then
		f_cut = (1.-10.*(r-R1)**3*(R2-R1)**2+15.*(r-R1)**4*(R2-R1)-6.*(r-R1))/(R2-R1)**5
	elseif(r<R1) then
		f_cut = 1.
	else
		f_cut = 0.
	endif
end function f_cut

function df_cut(r,R1,R2)
	real df_cut,r,R1,R2

	if(r>R1 .and. r<R2) then
		df_cut = (-30.*(r-R1)**2*(R2-R1)**2+60.*(r-R1)**3*(R2-R1)-6.)/(R2-R1)**5/r
	else
		df_cut = 0.
	endif	
end function df_cut

end module cut_off_poly