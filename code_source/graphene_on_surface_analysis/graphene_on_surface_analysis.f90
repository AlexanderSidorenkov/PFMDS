module graphene_on_surface_analysis
use md_general
use md_read_write
implicit none

contains

subroutine gr_on_cu_analysis(arr1,arr2,filename,z)
	type(particles)							:: atoms
	type(particle_group)					:: groupC,groupCU
	character(len=256)						:: filename
	character(len=32)						:: type_names(3)
	real									:: z,arr1(3),arr2(3)

	call read_particles(atoms,filename)
	type_names(1) = 'C';	type_names(2) = 'C_a';		type_names(3) = 'C_b'
	call create_particle_group(groupC,type_names,atoms)
	call position_analysis(arr1(1),arr1(2),arr1(3),atoms,groupC,3,z,1000.)
	type_names(1) = 'CU';	type_names(2) = 'CU_fixed';	type_names(3) = '#'
	call create_particle_group(groupCU,type_names,atoms)
	call position_analysis(arr2(1),arr2(2),arr2(3),atoms,groupCU,3,z,1000.)
	deallocate(atoms%positions)
	deallocate(atoms%velocities)
	deallocate(atoms%forces)
	deallocate(atoms%masses)
	deallocate(atoms%atom_types)
	deallocate(groupC%indexes)
	deallocate(groupCU%indexes)
	
end subroutine 

end module graphene_on_surface_analysis