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
	LJp%c6 = 4.*LJp%eps*LJp%sig**6
	LJp%c12 = 4.*LJp%eps*LJp%sig**12
	LJp%c6t6 = 6.*4.*LJp%eps*LJp%sig**6
	LJp%c12t12 = 12.*4.*LJp%eps*LJp%sig**12
	
end subroutine read_LJ1g_parameters

subroutine LJ1g_energy(energy,nl,LJp)
	type(neighbour_list):: nl
	type(LennardJones1g_parameters):: LJp
	integer:: i,p,chunk_size
	real:: energy,energy_priv,U
	
	call get_chunk_size(chunk_size)
	
	energy = 0.
	
	!$OMP PARALLEL private(i,p,U,energy_priv)
	energy_priv = 0.
	!$OMP DO schedule(dynamic,chunk_size)
	do i=1,nl%N
		do p=1,nl%lessnnum(i)!do p=1,nl%nnum(i)!if (nl%nlist(p,i)>i) exit
			U = 1./(nl%moddr(p,i)*nl%moddr(p,i)*nl%moddr(p,i)*nl%moddr(p,i)*nl%moddr(p,i)*nl%moddr(p,i))
			energy_priv = energy_priv+U*(LJp%c12*U-LJp%c6)*f_cut(nl%moddr(p,i),LJp%R1,LJp%R2)
		enddo
	enddo
	!$OMP END DO
	!$OMP ATOMIC
		energy = energy+energy_priv
	!$OMP END PARALLEL
	
end subroutine LJ1g_energy

subroutine LJ1g_forces(atoms,nl,LJp)
	type(particles) ::	atoms
	type(neighbour_list) :: nl
	type(LennardJones1g_parameters) :: LJp
	integer:: i,p,k,ind,jnd,chunk_size
	real,allocatable:: priv_force(:,:),F(:),fp(:,:)
	
	call get_chunk_size(chunk_size)
	
	!$OMP PARALLEL private(i,p,k,ind,jnd,F,fp,priv_force)
	if(.not. allocated(priv_force)) allocate(priv_force(3,atoms%N))	!realloc if N changed
	if(.not. allocated(F)) allocate(F(nl%neighb_num_max))			!realloc if neighb_num_max changed
	if(.not. allocated(fp)) allocate(fp(3,nl%neighb_num_max))		!realloc if neighb_num_max changed
	!F = 0.!?
	!fp = 0.!?
	priv_force = 0.
	
	!$OMP DO SCHEDULE(dynamic,chunk_size)	 
	do i=1,nl%N
		do p=1,nl%lessnnum(i)
			F(p) = scalar_lj_force(nl%moddr(p,i),LJp%R1,LJp%R2,LJp%c12,LJp%c6,LJp%c12t12,LJp%c6t6)
		enddo
		do p=1,nl%lessnnum(i)
			fp(1,p) = F(p)*nl%dr(1,p,i)
			fp(2,p) = F(p)*nl%dr(2,p,i)
			fp(3,p) = F(p)*nl%dr(3,p,i)
		enddo
		do p=1,nl%lessnnum(i)
			ind = nl%particle_index(i)
			jnd = nl%particle_index(nl%nlist(p,i))
			priv_force(1,ind) = priv_force(1,ind)-fp(1,p)
			priv_force(2,ind) = priv_force(2,ind)-fp(2,p)
			priv_force(3,ind) = priv_force(3,ind)-fp(3,p)
			priv_force(1,jnd) = priv_force(1,jnd)+fp(1,p)
			priv_force(2,jnd) = priv_force(2,jnd)+fp(2,p)
			priv_force(3,jnd) = priv_force(3,jnd)+fp(3,p)
		enddo
	enddo
	!$OMP END DO
	
	do i=1,atoms%N
		do k=1,3
			!$OMP ATOMIC
				atoms%forces(k,i) = atoms%forces(k,i)+priv_force(k,i)
		enddo
	enddo

	!$OMP END PARALLEL
	
end subroutine LJ1g_forces

function scalar_lj_force(r,R1,R2,c12,c6,c12t12,c6t6)
	!$OMP DECLARE SIMD(scalar_lj_force) UNIFORM(R1,R2,c12,c6,c12t12,c6t6)
	real, intent(in) :: r
	real, intent(in) :: R1,R2,c12,c6,c12t12,c6t6
	real scalar_lj_force
	real :: invr2,U,fcut,dfrcut
	
	invr2 = 1./(r*r)
	U = invr2*invr2*invr2
	call f_dfr_cut(fcut,dfrcut,r,R1,R2) !remove from here? make function?
	scalar_lj_force = U*invr2*((c12t12*U-c6t6)*fcut-(c12*U-c6)*dfrcut)
	
end function

end module LennardJones_1g