module LennardJones_1g
use perfomance_settings
use cut_off_function
use md_general
implicit none

type LennardJones1g_parameters
	real eps,sig,R1,R2,c6,c12,c6t6,c12t12
end type LennardJones1g_parameters

contains

subroutine read_LJ1g_parameters(LJp,filename)
	type(LennardJones1g_parameters):: LJp
	character(*):: 	filename
	
	open(1,file=filename)
	read(1,*) LJp%eps,LJp%sig
	read(1,*) LJp%R1,LJp%R2
	close(1)
	LJp%c6 = 4*LJp%eps*LJp%sig**6
	LJp%c12 = 4*LJp%eps*LJp%sig**12
	LJp%c6t6 = 6*4*LJp%eps*LJp%sig**6
	LJp%c12t12 = 12*4*LJp%eps*LJp%sig**12
	
end subroutine read_LJ1g_parameters

subroutine LJ1g_energy(energy,nl,LJp)
	type(neibour_list):: nl
	type(LennardJones1g_parameters):: LJp
	integer:: i,p
	real:: energy,energy_priv,U
	
	energy = 0.
	energy_priv = 0.
	!$OMP PARALLEL firstprivate(energy_priv,i,p,U)
	!$OMP DO
	do i=1,nl%N
		do p=1,nl%nnum(i)
			if (nl%moddr(p,i)<LJp%R2 .and. nl%nlist(p,i)>i) then
				U = 1./nl%moddr(p,i)**6
				energy_priv = energy_priv+U*(LJp%c12*U-LJp%c6)*f_cut(nl%moddr(p,i),LJp%R1,LJp%R2)
			endif	
		enddo
	enddo
	!$OMP END DO
	!$OMP ATOMIC
		energy = energy+energy_priv
	!$OMP END PARALLEL
	
end subroutine LJ1g_energy

subroutine LJ1g_forces(atoms,nl,LJp)
	type(particles)::	atoms
	type(neibour_list):: nl
	type(LennardJones1g_parameters):: LJp
	integer:: i,p,k
	real:: U,F(3)
	real,allocatable:: priv_force(:,:)
	
	!$OMP PARALLEL firstprivate(i,p,U,priv_force)
	if(.not. allocated(priv_force)) allocate(priv_force(3,atoms%N))
	priv_force = 0.
	!$OMP DO schedule(dynamic,chunk_size)
	do i=1,nl%N
		do p=1,nl%nnum(i)
			if (nl%moddr(p,i)<LJp%R2 .and. nl%nlist(p,i)>i) then
				U = 1./nl%moddr(p,i)**6
				F = (U*(LJp%c12t12*U-LJp%c6t6)/nl%moddr(p,i)**2*f_cut(nl%moddr(p,i),LJp%R1,LJp%R2)&
				-U*(LJp%c12*U-LJp%c6)*df_cut(nl%moddr(p,i),LJp%R1,LJp%R2))*nl%dr(:,p,i)
				priv_force(:,nl%particle_index(i)) = priv_force(:,nl%particle_index(i))-F
				priv_force(:,nl%particle_index(nl%nlist(p,i))) = priv_force(:,nl%particle_index(nl%nlist(p,i)))+F
			endif
		enddo
	enddo
	!$OMP END DO
	do i=1,nl%N
		do k=1,3
			!$OMP ATOMIC
				atoms%forces(k,i) = atoms%forces(k,i)+priv_force(k,i)
		enddo
	enddo
	!$OMP END PARALLEL

end subroutine LJ1g_forces

end module LennardJones_1g