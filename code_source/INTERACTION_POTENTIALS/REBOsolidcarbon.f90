module REBOsolidcarbon
use cut_off_function
use md_general
implicit none

type REBOsc_parameters
	real A,Q,alpha,B(3),beta(3),T,g(6),R1,R2
end type REBOsc_parameters

contains

subroutine read_REBOsc_parameters(REBOscp,filename)
	type(REBOsc_parameters):: REBOscp
	character(*):: 	filename

	open(1,file=filename)
	read(1,*) REBOscp%A,REBOscp%Q,REBOscp%alpha
	read(1,*) REBOscp%B
	read(1,*) REBOscp%beta
	read(1,*) REBOscp%T
	read(1,*) REBOscp%g
	read(1,*) REBOscp%R1,REBOscp%R2
	close(1)
	
end subroutine read_REBOsc_parameters

subroutine REBOsc_energy(energy,nl,REBOscp)
	type(neibour_list):: nl
	type(REBOsc_parameters):: REBOscp
	integer:: i,p,q,l,j
	real:: energy,energy_priv,bsp(nl%neib_num_max,nl%N),bdh(nl%neib_num_max,nl%N),aa,ab,ac,bc,cosine
	
	energy = 0.
	energy_priv = 0.
	!$OMP PARALLEL firstprivate(energy_priv,i,p,q,l,j,aa,ab,ac,bc,cosine)
	!$OMP DO
	do i=1,nl%N
		bsp(:,i) = 0. !not very optimal but hides a bug(?) in numerical forse calc in parallel
		bdh(:,i) = 0.
		do p=1,nl%nnum(i)
			!bsp(p,i) = 0.
			!bdh(p,i) = 0.
			if (nl%moddr(p,i)<REBOscp%R2) then
				j = nl%nlist(p,i)
				do q=1,nl%nnum(i)
					if(p/=q .and. nl%moddr(q,i)<REBOscp%R2) then
						cosine = sum(nl%dr(:,p,i)*nl%dr(:,q,i))/(nl%moddr(p,i)*nl%moddr(q,i))
						bsp(p,i) = bsp(p,i)+f_cut(nl%moddr(q,i),REBOscp%R1,REBOscp%R2)*(REBOscp%g(1)+REBOscp%g(2)*cosine+&
						REBOscp%g(3)*cosine**2+REBOscp%g(4)*cosine**3+REBOscp%g(5)*cosine**4+REBOscp%g(6)*cosine**5)
						do l=1,nl%nnum(j) !mb if j>i?
							if(nl%nlist(l,j)/=i .and. nl%moddr(l,j)<REBOscp%R2) then
								aa = nl%moddr(p,i)**2
								ab = sum(nl%dr(:,p,i)*nl%dr(:,q,i))
								ac = sum(nl%dr(:,p,i)*nl%dr(:,l,j))
								bc = sum(nl%dr(:,q,i)*nl%dr(:,l,j))
								bdh(p,i) = bdh(p,i)+f_cut(nl%moddr(q,i),REBOscp%R1,REBOscp%R2)*f_cut(nl%moddr(l,j),REBOscp%R1,REBOscp%R2)*&
								(1.-(aa*bc-ab*ac)**2/(aa*nl%moddr(q,i)**2-ab**2)/(aa*nl%moddr(l,j)**2-ac**2))
							endif
						enddo
					endif
				enddo
				bsp(p,i) = (1.+bsp(p,i))**(-0.5)
				bdh(p,i) = REBOscp%T*bdh(p,i)
			endif
		enddo
	enddo
	!$OMP END DO
	!$OMP DO
	do i=1,nl%N
		do p=1,nl%nnum(i)
			if (nl%moddr(p,i)<REBOscp%R2) then
				j = nl%nlist(p,i)
				if (j>i) then
					do q=1,nl%nnum(j)
						if(nl%nlist(q,j)==i) exit
					enddo
					energy_priv = energy_priv+f_cut(nl%moddr(p,i),REBOscp%R1,REBOscp%R2)*&
					((1+REBOscp%Q/nl%moddr(p,i))*REBOscp%A*exp(-REBOscp%alpha*nl%moddr(p,i))-&
					((bsp(p,i)+bsp(q,j))/2+bdh(p,i))*&
					(REBOscp%B(1)*exp(-REBOscp%beta(1)*nl%moddr(p,i))+&
					REBOscp%B(2)*exp(-REBOscp%beta(2)*nl%moddr(p,i))+&
					REBOscp%B(3)*exp(-REBOscp%beta(3)*nl%moddr(p,i))))
				endif
			endif
		enddo
	enddo
	!$OMP END DO
	!$OMP ATOMIC
		energy = energy+energy_priv
	!$OMP END PARALLEL
	
end subroutine REBOsc_energy
	
end module REBOsolidcarbon