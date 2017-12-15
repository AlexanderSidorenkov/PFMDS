module LennardJonesCosine
use cut_off_function
use md_general
implicit none

type LennardJonesCosine_parameters
	real eps,sig,delt,R1,R2
	logical:: simplified
	real,allocatable:: gr_norm(:,:)
end type LennardJonesCosine_parameters

contains

subroutine read_LJC_parameters(LJCp,filename)
	type(LennardJonesCosine_parameters):: LJCp
	character(*):: 	filename
	
	open(1,file=filename)
	read(1,*) LJCp%eps,LJCp%sig,LJCp%delt
	read(1,*) LJCp%R1,LJCp%R2
	read(1,*) LJCp%simplified
	close(1)
	
end subroutine read_LJC_parameters

subroutine LJC_energy(energy,nl,LJCp)
	type(neibour_list):: nl
	type(LennardJonesCosine_parameters):: LJCp
	integer:: i,p
	real:: energy,energy_priv,V1,V2,V3
	
	energy = 0.
	energy_priv = 0.
	!$OMP PARALLEL firstprivate(energy_priv,i,p,V1,V2,V3)
	!$OMP DO
	do i=1,nl%N
		do p=1,nl%nnum(i)
			if (nl%moddr(p,i)<LJCp%R2) then
				V2 = (LJCp%sig/nl%moddr(p,i))**6
				V1 = V2**2
				V3 = (abs(sum(LJCp%gr_norm(:,i)*nl%dr(:,p,i)))/nl%moddr(p,i))**LJCp%delt
				energy_priv = energy_priv+4*LJCp%eps*(V1-V2*V3)*f_cut(nl%moddr(p,i),LJCp%R1,LJCp%R2)
			endif	
		enddo
	enddo
	!$OMP END DO
	!$OMP ATOMIC
		energy = energy+energy_priv
	!$OMP END PARALLEL
	
end subroutine LJC_energy

subroutine LJC_forces_for_graphene(atoms,nl,nl_nn,LJCp)
	type(particles)::	atoms
	type(neibour_list):: nl,nl_nn
	type(LennardJonesCosine_parameters):: LJCp
	integer:: i,p,q,l1,l2,l3,j,k,nnum_nn
	real:: drj12(3),drj31(3),drj23(3),V1,V2,V3,f_c,df_c
	
	if (nl%N/=nl_nn%N .and. nl%N/=0) then; write(*,*) 'error: nl%N/=nl_nn%N',nl%N,nl_nn%N; stop; endif
	nnum_nn = 3
	if (LJCp%simplified) then; nnum_nn = 0; LJCp%gr_norm = 0.; LJCp%gr_norm(3,:) = 1.; endif!optimize
	
	!$OMP PARALLEL firstprivate(i,p,q,l1,l2,l3,j,k,drj12,drj31,drj23,V1,V2,V3,f_c,df_c)
	!$OMP DO
	do i=1,nl%N
		do p=1,nl%nnum(i)
			if (nl%moddr(p,i)<LJCp%R2) then
				V2 = (LJCp%sig/nl%moddr(p,i))**6
				V1 = V2**2
				V3 = (abs(sum(LJCp%gr_norm(:,i)*nl%dr(:,p,i)))/nl%moddr(p,i))**LJCp%delt
				f_c = f_cut(nl%moddr(p,i),LJCp%R1,LJCp%R2)
				df_c = df_cut(nl%moddr(p,i),LJCp%R1,LJCp%R2)
				atoms%forces(:,nl%particle_index(i)) = atoms%forces(:,nl%particle_index(i))-4.*LJCp%eps*&
				(((12.*V1-(6.+LJCp%delt)*V2*V3)/nl%moddr(p,i)**2*f_c-(V1-V2*V3)*df_c)*nl%dr(:,p,i)+&
				(LJCp%delt*V2*V3/sum(LJCp%gr_norm(:,i)*nl%dr(:,p,i))*f_c)*LJCp%gr_norm(:,i))
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
				if (nl%moddr(p,j)<LJCp%R2) then
					V2 = (LJCp%sig/nl%moddr(p,j))**6
					V1 = V2**2
					V3 = (abs(sum(LJCp%gr_norm(:,j)*nl%dr(:,p,j)))/nl%moddr(p,j))**LJCp%delt
					f_c = f_cut(nl%moddr(p,j),LJCp%R1,LJCp%R2)
					atoms%forces(:,nl%particle_index(i)) = atoms%forces(:,nl%particle_index(i))&
					-4.*LJCp%eps*LJCp%delt*V2*V3*f_c/(sum(drj12**2)*sum(drj31**2)-sum(drj12*drj31)**2)/sum(LJCp%gr_norm(:,j)*nl%dr(:,p,j))*&
					LJCp%gr_norm(:,j)*sum(nl%dr(:,p,j)*(drj12*sum(drj23*drj31)-drj31*sum(drj23*drj12)))
				endif
			enddo
		enddo
	enddo
	!$OMP END DO
	!$OMP END PARALLEL
	
end subroutine LJC_forces_for_graphene

subroutine LJC_forces_for_other_atoms(atoms,nl,LJCp)
	type(particles)::	atoms
	type(neibour_list):: nl
	type(LennardJonesCosine_parameters):: LJCp
	integer:: i,p
	real:: gr_n(3),V1,V2,V3,f_c,df_c
	
	!$OMP PARALLEL firstprivate(i,p,gr_n,V1,V2,V3,f_c,df_c)
	!$OMP DO
	do i=1,nl%N
		do p=1,nl%nnum(i)
			if (nl%moddr(p,i)<LJCp%R2) then
				gr_n = LJCp%gr_norm(:,nl%nlist(p,i))
				V2 = (LJCp%sig/nl%moddr(p,i))**6
				V1 = V2**2
				V3 = (abs(sum(gr_n*nl%dr(:,p,i)))/nl%moddr(p,i))**LJCp%delt
				f_c = f_cut(nl%moddr(p,i),LJCp%R1,LJCp%R2)
				df_c = df_cut(nl%moddr(p,i),LJCp%R1,LJCp%R2)
				atoms%forces(:,nl%particle_index(i)) = atoms%forces(:,nl%particle_index(i))-4.*LJCp%eps*&
				(((12.*V1-(6.+LJCp%delt)*V2*V3)/nl%moddr(p,i)**2*f_c-(V1-V2*V3)*df_c)*nl%dr(:,p,i)+&
				(LJCp%delt*V2*V3/sum(gr_n*nl%dr(:,p,i))*f_c)*gr_n)
			endif
		enddo
	enddo
	!$OMP END DO 
	!$OMP END PARALLEL
	
end subroutine LJC_forces_for_other_atoms

end module LennardJonesCosine