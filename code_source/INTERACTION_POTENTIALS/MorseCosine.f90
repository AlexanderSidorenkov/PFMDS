module MorseCosine
use cut_off_function
use md_general
implicit none

type MorseCosine_parameters
	real:: d,r,a,delt,R1,R2
	logical:: simplified
	real,allocatable:: gr_norm(:,:)
end type MorseCosine_parameters

contains

subroutine read_MorseC_parameters(MorseCp,filename)
	type(MorseCosine_parameters):: MorseCp
	character(*):: 	filename
	
	open(1,file=filename)
	read(1,*) MorseCp%d,MorseCp%r,MorseCp%a,MorseCp%delt
	read(1,*) MorseCp%R1,MorseCp%R2
	read(1,*) MorseCp%simplified
	close(1)
	
end subroutine read_MorseC_parameters

subroutine MorseC_energy(energy,nl,MorseCp)
	type(neibour_list):: nl
	type(MorseCosine_parameters):: MorseCp
	integer:: i,p
	real:: energy,energy_priv,V1,V2,V3
	
	energy = 0.
	energy_priv = 0.
	!$OMP PARALLEL firstprivate(energy_priv,i,p,V1,V2,V3)
	!$OMP DO
	do i=1,nl%N
		do p=1,nl%nnum(i)
			if (nl%moddr(p,i)<MorseCp%R2) then
				V2 = exp(-MorseCp%a*(nl%moddr(p,i)-MorseCp%r))
				V1 = V2**2
				V3 = (abs(sum(MorseCp%gr_norm(:,i)*nl%dr(:,p,i)))/nl%moddr(p,i))**MorseCp%delt
				energy_priv = energy_priv+MorseCp%d*(V1-2.*V2*V3)*f_cut(nl%moddr(p,i),MorseCp%R1,MorseCp%R2)
			endif	
		enddo
	enddo
	!$OMP END DO
	!$OMP ATOMIC
		energy = energy+energy_priv
	!$OMP END PARALLEL
	
end subroutine MorseC_energy

subroutine MorseC_forces_for_graphene(atoms,nl,nl_nn,MorseCp)
	type(particles)::	atoms
	type(neibour_list):: nl,nl_nn
	type(MorseCosine_parameters):: MorseCp
	integer:: i,p,q,l1,l2,l3,j,k,nnum_nn
	real:: drj12(3),drj31(3),drj23(3),V1,V2,V3,f_c,df_c
	
	if (nl%N/=nl_nn%N .and. nl%N/=0) then; write(*,*) 'error: nl%N/=nl_nn%N',nl%N,nl_nn%N; stop; endif
	nnum_nn = 3
	if (MorseCp%simplified) then; nnum_nn = 0; MorseCp%gr_norm = 0.; MorseCp%gr_norm(3,:) = 1.; endif!optimize
	
	!$OMP PARALLEL firstprivate(i,p,q,l1,l2,l3,j,k,drj12,drj31,drj23,V1,V2,V3,f_c,df_c)
	!$OMP DO
	do i=1,nl%N
		do p=1,nl%nnum(i)
			if (nl%moddr(p,i)<MorseCp%R2) then
				V2 = exp(-MorseCp%a*(nl%moddr(p,i)-MorseCp%r))
				V1 = V2**2
				V3 = (abs(sum(MorseCp%gr_norm(:,i)*nl%dr(:,p,i)))/nl%moddr(p,i))**MorseCp%delt
				f_c = f_cut(nl%moddr(p,i),MorseCp%R1,MorseCp%R2)
				df_c = df_cut(nl%moddr(p,i),MorseCp%R1,MorseCp%R2)
				atoms%forces(:,nl%particle_index(i)) = atoms%forces(:,nl%particle_index(i))-MorseCp%d*&
				((2.*(MorseCp%a*V1-(MorseCp%a+MorseCp%delt/nl%moddr(p,i))*V2*V3)/nl%moddr(p,i)*f_c-(V1-2.*V2*V3)*df_c)*nl%dr(:,p,i)+&
				(2.*MorseCp%delt*V2*V3/sum(MorseCp%gr_norm(:,i)*nl%dr(:,p,i))*f_c)*MorseCp%gr_norm(:,i))
			endif
		enddo
	enddo
	!$OMP END DO
	!$OMP DO
	do i=1,nl%N
		do q=1,nnum_nn
			j = nl_nn%nlist(q,i)
			do l1=1,nnum_nn
				if (nl_nn%nlist(l1,j)==i) exit
			enddo
			if(l1>3) then; write(*,*) 'l1>3'; stop; endif
			l2 = mod(l1,3)+1
			l3 = mod(l1+1,3)+1
			drj12 = nl_nn%dr(:,l2,j)-nl_nn%dr(:,l1,j)
			drj31 = nl_nn%dr(:,l1,j)-nl_nn%dr(:,l3,j)
			drj23 = nl_nn%dr(:,l3,j)-nl_nn%dr(:,l2,j)
			do p=1,nl%nnum(j)
				if (nl%moddr(p,j)<MorseCp%R2) then
					V2 = exp(-MorseCp%a*(nl%moddr(p,j)-MorseCp%r))
					V1 = V2**2
					V3 = (abs(sum(MorseCp%gr_norm(:,j)*nl%dr(:,p,j)))/nl%moddr(p,j))**MorseCp%delt
					f_c = f_cut(nl%moddr(p,j),MorseCp%R1,MorseCp%R2)
					atoms%forces(:,nl%particle_index(i)) = atoms%forces(:,nl%particle_index(i))&
					-2.*MorseCp%d*MorseCp%delt*V2*V3*f_c/(sum(drj12**2)*sum(drj31**2)-sum(drj12*drj31)**2)/sum(MorseCp%gr_norm(:,j)*nl%dr(:,p,j))*&
					MorseCp%gr_norm(:,j)*sum(nl%dr(:,p,j)*(drj12*sum(drj23*drj31)-drj31*sum(drj23*drj12)))
				endif
			enddo
		enddo
	enddo
	!$OMP END DO
	!$OMP END PARALLEL
	
end subroutine MorseC_forces_for_graphene

subroutine MorseC_forces_for_other_atoms(atoms,nl,MorseCp)
	type(particles)::	atoms
	type(neibour_list):: nl
	type(MorseCosine_parameters):: MorseCp
	integer:: i,p
	real:: gr_n(3),V1,V2,V3,f_c,df_c
	
	!$OMP PARALLEL firstprivate(i,p,gr_n,V1,V2,V3,f_c,df_c)
	!$OMP DO
	do i=1,nl%N
		do p=1,nl%nnum(i)
			if (nl%moddr(p,i)<MorseCp%R2) then
				gr_n = MorseCp%gr_norm(:,nl%nlist(p,i))
				V2 = exp(-MorseCp%a*(nl%moddr(p,i)-MorseCp%r))
				V1 = V2**2
				V3 = (abs(sum(gr_n*nl%dr(:,p,i)))/(sqrt(sum(gr_n**2))*nl%moddr(p,i)))**MorseCp%delt
				f_c = f_cut(nl%moddr(p,i),MorseCp%R1,MorseCp%R2)
				df_c = df_cut(nl%moddr(p,i),MorseCp%R1,MorseCp%R2)
				atoms%forces(:,nl%particle_index(i)) = atoms%forces(:,nl%particle_index(i))-MorseCp%d*&
				((2.*(MorseCp%a*V1-(MorseCp%a+MorseCp%delt/nl%moddr(p,i))*V2*V3)/nl%moddr(p,i)*f_c-(V1-2.*V2*V3)*df_c)*nl%dr(:,p,i)+&
				(2.*MorseCp%delt*V2*V3/sum(gr_n*nl%dr(:,p,i))*f_c)*gr_n)
			endif
		enddo
	enddo
	!$OMP END DO 
	!$OMP END PARALLEL
	
end subroutine MorseC_forces_for_other_atoms

end module MorseCosine