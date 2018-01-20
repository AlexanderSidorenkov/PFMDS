module LennardJones
use cut_off_function
use md_general
implicit none

type LennardJones_parameters
	real eps,sig,R1,R2
end type LennardJones_parameters

contains

subroutine read_LJ_parameters(LJp,filename)
	type(LennardJones_parameters):: LJp
	character(*):: 	filename
	
	open(1,file=filename)
	read(1,*) LJp%eps,LJp%sig
	read(1,*) LJp%R1,LJp%R2
	close(1)
	
end subroutine read_LJ_parameters

subroutine LJ_energy(energy,nl,LJp)
	type(neighbour_list):: nl
	type(LennardJones_parameters):: LJp
	integer:: i,p
	real:: energy,energy_priv,V
	
	energy = 0.
	energy_priv = 0.
	!$OMP PARALLEL firstprivate(energy_priv,i,p,V)
	!$OMP DO
	do i=1,nl%N
		do p=1,nl%nnum(i)
			if (nl%moddr(p,i)<LJp%R2) then
				V = (LJp%sig/nl%moddr(p,i))**6
				energy_priv = energy_priv+4*LJp%eps*V*(V-1.)*f_cut(nl%moddr(p,i),LJp%R1,LJp%R2)
			endif	
		enddo
	enddo
	!$OMP END DO
	!$OMP ATOMIC
		energy = energy+energy_priv
	!$OMP END PARALLEL
	
end subroutine LJ_energy

subroutine LJ_forces(atoms,nl,LJp)
	type(particles)::	atoms
	type(neighbour_list):: nl
	type(LennardJones_parameters):: LJp
	integer:: i,p
	real:: V
	
	!$OMP PARALLEL firstprivate(i,p,V)
	!$OMP DO
	do i=1,nl%N
		do p=1,nl%nnum(i)
			if (nl%moddr(p,i)<LJp%R2) then
				V = (LJp%sig/nl%moddr(p,i))**6
				atoms%forces(:,nl%particle_index(i)) = atoms%forces(:,nl%particle_index(i))-4.*LJp%eps*&
				(V*(12.*V-6.)/nl%moddr(p,i)**2*f_cut(nl%moddr(p,i),LJp%R1,LJp%R2)-V*(V-1.)*df_cut(nl%moddr(p,i),LJp%R1,LJp%R2))*nl%dr(:,p,i)
			endif
		enddo
	enddo
	!$OMP END DO
	!$OMP END PARALLEL
	
end subroutine LJ_forces

end module LennardJones