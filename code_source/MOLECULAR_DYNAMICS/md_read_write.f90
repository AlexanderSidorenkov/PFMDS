module md_read_write
use md_general
implicit none

contains

subroutine read_integrator_params(integr,file_id)
	type(integrator_params):: integr
	integer:: file_id

	read(file_id,*) integr%int_name,integr%dt,integr%l,integr%period_snapshot,integr%period_log

end subroutine read_integrator_params

subroutine read_box_size(box,filename)
	type(simulation_cell):: box
	real:: boxmatrix(9)
	character(*):: 	filename
	character(len=32)::	c

	open(1,file=filename)
	read(1,*)
	read(1,*) c,boxmatrix
	close(1)
	box%box_size(1)=boxmatrix(1)
	box%box_size(2)=boxmatrix(5)
	box%box_size(3)=boxmatrix(9)
	box%half_box_size = 0.5*box%box_size
	
	return
end subroutine read_box_size

subroutine read_particles(atoms,filename)
	type(particles)::	atoms
	character(*):: 	filename
	integer::		i

	open(1,file=filename)
	read(1,*) atoms%N
	allocate(atoms%positions(3,atoms%N))
	allocate(atoms%velocities(3,atoms%N))
	allocate(atoms%forces(3,atoms%N))
	allocate(atoms%masses(atoms%N))
	allocate(atoms%atom_types(atoms%N))
	read(1,*)
	do i=1,atoms%N
		read(1,*) atoms%positions(:,i),atoms%velocities(:,i),atoms%masses(i),atoms%atom_types(i)
	enddo
	close(1)

	return
end subroutine read_particles

subroutine write_particle_group(filename,atoms,group,box)
	type(particles)::	atoms
	type(particle_group)::	group
	type(simulation_cell):: box
	character(*):: 	filename
	integer::		i,ind

	open(2,file=filename)
	write(2,*) group%N
	write(2,'(A,9f16.6,A)') 'Lattice="',box%box_size(1), 0., 0., 0., box%box_size(2), 0., 0., 0., box%box_size(3),&
			' " Properties=pos:R:3:vel:R:3:mass:R:1:species:S:1'
	do ind=1,group%N
		i = group%indexes(ind)
		write(2,'(7f27.16,A,A)') atoms%positions(:,i),atoms%velocities(:,i),atoms%masses(i),'    ',atoms%atom_types(i)
	enddo
	close(2)

	return
end subroutine write_particle_group

subroutine write_particle_group_append(filename,atoms,group,box,md_step)
	type(particles)::	atoms
	type(particle_group)::	group
	type(simulation_cell):: box
	character(*):: 	filename
	integer::		i,ind,md_step

	open(2,file=filename,action='write',position='append')
	write(2,*) group%N
	write(2,'(A,i9,A,9f10.4,A)') 'time_step: ',md_step,&
	'    Lattice="',box%box_size(1), 0., 0., 0., box%box_size(2), 0., 0., 0., box%box_size(3),&
	' " Properties=pos:R:3:vel:R:3:mass:R:1:species:S:1'
	do ind=1,group%N
		i = group%indexes(ind)
		write(2,'(7f10.4,A,A)') atoms%positions(:,i),atoms%velocities(:,i),atoms%masses(i),'    ',trim(atoms%atom_types(i))
	enddo
	close(2)
	
	return
end subroutine write_particle_group_append

end module md_read_write