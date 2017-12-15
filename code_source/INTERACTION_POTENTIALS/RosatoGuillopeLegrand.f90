module RosatoGuillopeLegrand
use cut_off_function
use md_general
implicit none

type RosatoGuillopeLegrand_parameters
	real A0,xi,p,q,r0,R1,R2
end type RosatoGuillopeLegrand_parameters

contains

subroutine read_RJL_parameters(RJLp,filename)
	type(RosatoGuillopeLegrand_parameters):: RJLp
	character(*):: 	filename
	
	open(1,file=filename)
	read(1,*) RJLp%A0,RJLp%xi,RJLp%p,RJLp%q,RJLp%r0
	read(1,*) RJLp%R1,RJLp%R2
	close(1)
	
end subroutine read_RJL_parameters

subroutine RJL_energy(energy,nl,RJLp)
	type(RosatoGuillopeLegrand_parameters):: RJLp
	type(neibour_list):: nl
	integer:: i,p
	real:: energy,energy_priv,Eb2,Er
	
	energy = 0.
	energy_priv = 0.
	!$OMP PARALLEL firstprivate(energy_priv,i,p,Eb2,Er)
	!$OMP DO
	do i=1,nl%N
		Eb2 = 0.
		Er = 0.
		do p=1,nl%nnum(i)
			if( nl%moddr(p,i)<=RJLp%R2 ) then
				Eb2 = Eb2+exp(-2.*RJLp%q*(nl%moddr(p,i)/RJLp%r0-1.))*f_cut(nl%moddr(p,i),RJLp%R1,RJLp%R2)
				Er = Er+exp(-RJLp%p*(nl%moddr(p,i)/RJLp%r0-1.))*f_cut(nl%moddr(p,i),RJLp%R1,RJLp%R2)
			endif
		enddo
		energy_priv = energy_priv+RJLp%A0*Er-RJLp%xi*sqrt(Eb2)
	enddo
	!$OMP END DO
	!$OMP ATOMIC
		energy = energy+energy_priv
	!$OMP END PARALLEL
	
end subroutine RJL_energy

subroutine RJL_forces(atoms,nl,RJLp)
	type(RosatoGuillopeLegrand_parameters):: RJLp
	type(particles)::	atoms
	type(neibour_list):: nl
	integer:: i,p
	real:: fc,dfc,Eb(nl%N),expp(nl%neib_num_max,nl%N),expq(nl%neib_num_max,nl%N)

	!$OMP PARALLEL firstprivate(i,p,fc,dfc)
	!$OMP DO
	do i=1,nl%N
		do p=1,nl%nnum(i)
			if( nl%moddr(p,i)<=RJLp%R2 ) then
				expq(p,i) = exp(-2.*RJLp%q*(nl%moddr(p,i)/RJLp%r0-1.))
				expp(p,i) = exp(-RJLp%p*(nl%moddr(p,i)/RJLp%r0-1.))
			endif
		enddo
	enddo
	!$OMP END DO
	!$OMP DO	
	do i=1,nl%N
		Eb(i) = 0.
		do p=1,nl%nnum(i)
			if( nl%moddr(p,i)<=RJLp%R2 ) then
				Eb(i) = Eb(i)+expq(p,i)*f_cut(nl%moddr(p,i),RJLp%R1,RJLp%R2)
			endif
		enddo
		Eb(i) = sqrt(Eb(i))
	enddo
	!$OMP END DO
	!$OMP DO	
	do i=1,nl%N
		do p=1,nl%nnum(i)
			if( nl%moddr(p,i)<=RJLp%R2 ) then
				fc = f_cut(nl%moddr(p,i),RJLp%R1,RJLp%R2)
				dfc = df_cut(nl%moddr(p,i),RJLp%R1,RJLp%R2)
				atoms%forces(:,nl%particle_index(i)) = atoms%forces(:,nl%particle_index(i))&
				-(2.*RJLp%A0*(RJLp%p/RJLp%r0*fc-dfc*nl%moddr(p,i))*expp(p,i)&
				-RJLp%xi*(RJLp%q/RJLp%r0*fc-dfc/2.*nl%moddr(p,i))*(1./Eb(i)+1./Eb(nl%nlist(p,i)))*expq(p,i))/nl%moddr(p,i)*nl%dr(:,p,i)
			endif
		enddo
	enddo
	!$OMP END DO
	!$OMP END PARALLEL
end subroutine RJL_forces
	
end module RosatoGuillopeLegrand