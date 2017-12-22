module LennardJones_1g
use perfomance_settings
use cut_off_poly
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
	!$OMP PARALLEL firstprivate(energy_priv) private(i,p,U)
	!$OMP DO schedule(dynamic,chunk_size)
	do i=1,nl%N
		do p=1,nl%nnum(i)
			if (nl%nlist(p,i)>i) exit
			!if (nl%nlist(p,i)<i) then
				U = 1./nl%moddr(p,i)**6
				energy_priv = energy_priv+U*(LJp%c12*U-LJp%c6)*f_cut(nl%moddr(p,i),LJp%R1,LJp%R2)
			!endif	
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
	integer:: i,p,k,ind,jnd
	real:: U,F,invr2
	real,allocatable:: priv_force(:,:)
	
	!$OMP PARALLEL private(i,p,k,U,F,priv_force,invr2)
	if(.not. allocated(priv_force)) allocate(priv_force(3,atoms%N))
	priv_force = 0.
	!$OMP DO schedule(dynamic,chunk_size)
	do i=1,nl%N
		do p=1,nl%nnum(i)
			if (nl%nlist(p,i)>i) exit
			!if (nl%nlist(p,i)<i) then
				invr2 = 1./nl%moddr(p,i)/nl%moddr(p,i)
				U = invr2*invr2*invr2
				F = (U*(LJp%c12t12*U-LJp%c6t6)*invr2*f_cut(nl%moddr(p,i),LJp%R1,LJp%R2)&
				-U*(LJp%c12*U-LJp%c6)*df_cut(nl%moddr(p,i),LJp%R1,LJp%R2))
				ind = nl%particle_index(i)
				do k=1,3
					priv_force(k,ind) = priv_force(k,ind)-F*nl%dr(k,p,i)
				enddo
				jnd = nl%particle_index(nl%nlist(p,i))
				do k=1,3
					priv_force(k,jnd) = priv_force(k,jnd)+F*nl%dr(k,p,i)
				enddo
			!endif
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