module md_integrators
use md_general
implicit none

contains

subroutine integrate_verlet_xyz_positions(atoms,group,steps,box)
	type(particles)::	atoms
	type(particle_group):: group
	type(time_steps):: steps
	type(simulation_cell):: box
	integer::		i,k,ind

	!$OMP PARALLEL firstprivate(i,ind,k)
	!$OMP DO
	do ind=1,group%N
		i = group%indexes(ind)
		do k=1,3
			atoms%positions(k,i) = atoms%positions(k,i)+atoms%velocities(k,i)*steps%ts(1)
			if (atoms%positions(k,i)>box%box_size(k)) then
				atoms%positions(k,i) = atoms%positions(k,i)-box%box_size(k)
			elseif (atoms%positions(k,i)<0.) then
				atoms%positions(k,i) = atoms%positions(k,i)+box%box_size(k)
			endif
		enddo
	enddo
	!$OMP END DO
	!$OMP END PARALLEL

	return
end subroutine integrate_verlet_xyz_positions

subroutine integrate_verlet_z_positions(atoms,group,steps,box)
	type(particles)::	atoms
	type(particle_group):: group
	type(time_steps):: steps
	type(simulation_cell):: box
	integer::		i,k,ind

	!$OMP PARALLEL firstprivate(i,ind,k)
	!$OMP DO
	do ind=1,group%N
		i = group%indexes(ind)
		k = 3
			atoms%positions(k,i) = atoms%positions(k,i)+atoms%velocities(k,i)*steps%ts(1)
			if (atoms%positions(k,i)>box%box_size(k)) then
				atoms%positions(k,i) = atoms%positions(k,i)-box%box_size(k)
			elseif (atoms%positions(k,i)<0.) then
				atoms%positions(k,i) = atoms%positions(k,i)+box%box_size(k)
			endif
	enddo
	!$OMP END DO
	!$OMP END PARALLEL

	return
end subroutine integrate_verlet_z_positions

subroutine integrate_verlet_xyz_velocities(atoms,group,steps)
	type(particles)::	atoms
	type(particle_group):: group
	type(time_steps):: steps
	real,parameter:: mass_coef=1.6605389217/1.6021765654*10.**(2)
	integer::		i,ind,k

	!$OMP PARALLEL firstprivate(i,ind,k)
	!$OMP DO
	do ind=1,group%N
		i = group%indexes(ind)
		do k=1,3
			atoms%velocities(k,i) = atoms%velocities(k,i)+atoms%forces(k,i)/atoms%masses(i)/mass_coef*steps%ts(2)
		enddo
	enddo
	!$OMP END DO
	!$OMP END PARALLEL

	return
end subroutine integrate_verlet_xyz_velocities

subroutine integrate_verlet_z_velocities(atoms,group,steps)
	type(particles)::	atoms
	type(particle_group):: group
	type(time_steps):: steps
	real,parameter:: mass_coef=1.6605389217/1.6021765654*10.**(2)
	integer::		i,ind,k

	!$OMP PARALLEL firstprivate(i,ind,k)
	!$OMP DO
	do ind=1,group%N
		i = group%indexes(ind)
		k = 3
			atoms%velocities(k,i) = atoms%velocities(k,i)+atoms%forces(k,i)/atoms%masses(i)/mass_coef*steps%ts(2)
	enddo
	!$OMP END DO
	!$OMP END PARALLEL

	return
end subroutine integrate_verlet_z_velocities

subroutine molecular_static_xyz_velocities(atoms,group)
	type(particles)::	atoms
	type(particle_group):: group
	real:: fv,ff
	integer:: i,k,ind

	!$OMP PARALLEL firstprivate(i,k,ind,fv,ff)
	!$OMP DO
	do ind=1,group%N
		i = group%indexes(ind)
		fv = sum(atoms%forces(:,i)*atoms%velocities(:,i))
		ff = sum(atoms%forces(:,i)**2)
		if (fv>0. .and. ff>10.**(-12) ) then
			do k=1,3
				atoms%velocities(k,i) = fv/ff*atoms%forces(k,i)
			enddo
		else
			atoms%velocities(:,i) = 0.
		endif
	enddo
	!$OMP END DO
	!$OMP END PARALLEL

	return
end subroutine molecular_static_xyz_velocities

subroutine molecular_static_1D_velocities(atoms,group)
	type(particles)::	atoms
	type(particle_group):: group
	real:: fv
	integer:: i,k,ind

	!$OMP PARALLEL firstprivate(i,k,ind,fv)
	!$OMP DO
	do ind=1,group%N
		i = group%indexes(ind)
		fv = sum(atoms%forces(:,i)*atoms%velocities(:,i))
		if (fv>0.) then
		else
			atoms%velocities(:,i) = 0.
		endif
	enddo
	!$OMP END DO
	!$OMP END PARALLEL

	return
end subroutine molecular_static_1D_velocities

subroutine zero_forces(atoms,group)
	type(particles)::	atoms
	type(particle_group):: group
	integer:: i,ind,k

	!$OMP PARALLEL firstprivate(i,ind,k)
	!$OMP DO
	do ind=1,group%N
		i = group%indexes(ind)
		do k=1,3
			atoms%forces(k,i) = 0.
		enddo
	enddo
	!$OMP END DO
	!$OMP END PARALLEL
	
end subroutine zero_forces

subroutine create_nose_hoover_chain(nhc,nhc_M)
	type(nose_hoover_chain):: nhc
	integer:: nhc_M
	
	allocate(nhc%x(nhc_M))
	allocate(nhc%v(nhc_M))
	allocate(nhc%q(nhc_M))
	nhc%M = nhc_M
	nhc%x = 0.
	nhc%v = 0.
	nhc%q = 0.
	nhc%e = 0.
	nhc%s = 1.
	
	return
end subroutine create_nose_hoover_chain

subroutine set_nose_hoover_chain(nhc,temp,q1,gn,l)
	type(nose_hoover_chain):: nhc
	integer:: l,i,gn
	real::	q1,temp

	if (l<1 .or. q1<0. .or. temp<0.) then; write(*,*) 'error: wrong nhc parameters'; stop; endif

	nhc%group_num = gn
	nhc%L = l
	nhc%temperature = temp
	nhc%q(1) = q1
	do i=2,nhc%M,1
		nhc%q(i) = nhc%q(1)/(3.*nhc%L)
	enddo
	
	return
end subroutine set_nose_hoover_chain

subroutine integrate_nose_hoover_chain(nhc,atoms,group,dt)
	type(particles)::	atoms
	type(particle_group):: group
	type(time_steps):: dt
	type(nose_hoover_chain):: nhc
	real:: kedif,ke,kt,b
	integer::		i

	call calculate_kinetic_energy(ke,atoms,group)
	!or !call calculate_temperature_wo_com_motion(temp,kewocomm,atoms,group)
	
	kt=1.3806488/1.6021765654*10.**(-4)*nhc%temperature
	kedif = 2.*ke-3.*nhc%L*kt

	if (nhc%M==1) then
		nhc%v(1)=nhc%v(1)+kedif/nhc%q(1)*dt%ts(3)
	else
		nhc%v(nhc%M)=nhc%v(nhc%M)+(nhc%q(nhc%M-1)*nhc%v(nhc%M-1)**2-kt)/nhc%q(nhc%M)*dt%ts(3)
		do i=nhc%M-1,2,-1
			b=exp(-nhc%v(i+1)*dt%ts(4))
			nhc%v(i)=nhc%v(i)*b**2+(nhc%q(i-1)*nhc%v(i-1)**2-kt)/nhc%q(i)*dt%ts(3)*b
		enddo
		b=exp(-nhc%v(2)*dt%ts(4))
		nhc%v(1)=nhc%v(1)*b**2+kedif/nhc%q(1)*dt%ts(3)*b
	endif

	nhc%s=exp(-nhc%v(1)*dt%ts(2))
	call scale_velocities(atoms,group,nhc%s)
	kedif = 2.*ke*nhc%s**2-3.*nhc%L*kt
	do i=1,nhc%M,1
		nhc%x(i)=nhc%x(i)+nhc%v(i)*dt%ts(2)
	enddo

	if (nhc%M==1) then
		nhc%v(1)=nhc%v(1)+kedif/nhc%q(1)*dt%ts(3)
	else
		nhc%v(1)=nhc%v(1)*b**2+kedif/nhc%q(1)*dt%ts(3)*b
		do i=2,nhc%M-1,1
			b=exp(-nhc%v(i+1)*dt%ts(4))
			nhc%v(i)=nhc%v(i)*b**2+(nhc%q(i-1)*nhc%v(i-1)**2-kt)/nhc%q(i)*dt%ts(3)*b
		enddo
		nhc%v(nhc%M)=nhc%v(nhc%M)+(nhc%q(nhc%M-1)*nhc%v(nhc%M-1)**2-kt)/nhc%q(nhc%M)*dt%ts(3)
	endif

	return
end subroutine integrate_nose_hoover_chain

subroutine calculate_nose_hoover_chain_energy(nhc)
	type(nose_hoover_chain):: nhc
	real:: kt
	integer::		i

	kt=1.3806488/1.6021765654*10.**(-4)*nhc%temperature

	nhc%e = nhc%q(1)/2*nhc%v(1)**2+3.*nhc%L*kt*nhc%x(1)
	do i=2,nhc%M,1
		nhc%e = nhc%e+nhc%q(i)/2*nhc%v(i)**2+kt*nhc%x(i)
	enddo

	return
end subroutine calculate_nose_hoover_chain_energy

end module md_integrators