module md_simulation
use md_general
use md_interactions
implicit none 
contains
subroutine md(out_id,all_out_id,input_path,settings_filename,output_prefix,out_period,num_of_omp_treads)

type(simulation_cell)					:: cell
type(time_steps)						:: dt
type(particles)							:: atoms
type(particle_group),allocatable		:: groups(:)
type(interaction),allocatable			:: interactions(:)
type(integrator_params),allocatable		:: integrators(:)
type(nose_hoover_chain)					:: nhc
real									::	starttime,mdtime,conserved_energy,total_energy,kinetic_energy,potential_energy,prev_potential_energy,&
											ms_de,fs(3),initial_temperature,target_temperature,temperature,nhc_q1
integer									::	i,num_of_omp_treads,md_step,time_step_limit,integrators_num,integrator_index,&
											file_id,out_id,log_id,out_period,all_out_id,&
											all_atoms_group_num,termo_atoms_group_num,moving_atoms_group_num
character(len=128)						::	str,filename,logfilename,init_xyz_filename,input_path,settings_filename,output_prefix
character(len=32)						::	integrator_name,interactions_energies_format
logical 								::	new_velocities
real									::	omp_get_wtime
	
starttime = omp_get_wtime()
log_id = 108
file_id = 2017
write(out_id,'(A,i6,A)') 'RUNNING ON ',num_of_omp_treads,' OPENMP THREADS'
call omp_set_num_threads(num_of_omp_treads)
!call omp_thread_limit(num_of_omp_treads)

open(file_id,file=trim(input_path)//settings_filename);	write(out_id,'(A,A)') 'settings_filename: ',trim(settings_filename)
read(file_id,*) str,time_step_limit;					write(out_id,'(A32,i12)') str,time_step_limit
read(file_id,*) str,logfilename;						write(out_id,'(A32,A,A)') str,'	',trim(logfilename)
read(file_id,*) str,init_xyz_filename;					write(out_id,'(A32,A,A)') str,'	',trim(init_xyz_filename)
call read_box_size(cell,trim(input_path)//init_xyz_filename);	write(out_id,'(A,3f16.6)') 	'box_size: ',cell%box_size
call read_particles(atoms,trim(input_path)//init_xyz_filename);	write(out_id,'(A,i12)') 	'particles_num: ',atoms%N
read(file_id,*) str,new_velocities;						write(out_id,'(A32,l8)') str,new_velocities
call create_groups(groups,file_id,out_id,atoms)
read(file_id,*) str,moving_atoms_group_num;				write(out_id,'(A32,i12)') str,moving_atoms_group_num
read(file_id,*) str,termo_atoms_group_num;				write(out_id,'(A32,i12)') str,termo_atoms_group_num
read(file_id,*) str,all_atoms_group_num;				write(out_id,'(A32,i12)') str,all_atoms_group_num
read(file_id,*) str,integrators_num;					write(out_id,'(A32,i12)') str,integrators_num
allocate(integrators(0:integrators_num)); integrators(0)%int_name='none'; integrators(0)%dt=0.; integrators(0)%l=0
read(file_id,'(A)') str;								write(out_id,'(A32)') str
do i=1,integrators_num
	call read_integrator_params(integrators(i),file_id)
	write(out_id,'(A,A6,f10.5,4i9)') '  ',integrators(i)%int_name,&
	integrators(i)%dt,integrators(i)%l,integrators(i)%period_snapshot,integrators(i)%period_log
enddo
read(file_id,*) str,ms_de;								write(out_id,'(A32,f16.6)') str,ms_de
read(file_id,*) str,target_temperature,nhc%M,nhc_q1;	write(out_id,'(A32,f16.6,i6,f16.6)') str,target_temperature,nhc%M,nhc_q1
call create_nose_hoover_chain(nhc)
read(file_id,*) str,initial_temperature;				write(out_id,'(A32,f16.6)') str,initial_temperature
if (initial_temperature<0.) initial_temperature=0.
if (new_velocities) call set_new_temperature(atoms,groups(moving_atoms_group_num),initial_temperature)
call create_interactions(interactions,groups,file_id,out_id,input_path)
close(file_id)

write(interactions_energies_format,'("(",i0,"f20.9)")') size(interactions)
write(str,'(A,A)') trim(output_prefix),trim(logfilename)
open(log_id,file=str)
write(out_id,*) 'MD START TIME: ',omp_get_wtime()-starttime,' SECONDS'
mdtime = omp_get_wtime()
integrator_index = 0
dt%simulation_time = 0.
do md_step=0,time_step_limit

	if (md_step==sum(integrators(0:integrator_index)%l)) then
		integrator_index = integrator_index +1
		do i=integrator_index,integrators_num
			if (integrators(i)%l>0) then
				integrator_index = i; integrator_name = integrators(i)%int_name; call init_time_steps(dt,integrators(i)%dt)
				if (integrator_name=='nvt') call set_nose_hoover_chain(nhc,target_temperature,nhc_q1,groups(termo_atoms_group_num)%N)
				exit
			endif
		enddo
		if (integrator_index>integrators_num) exit
	endif
	
	call check_positions(out_id,atoms,cell)
	
	if (md_step/=0) then
		if (integrator_name=='nvms') then
			if (abs(potential_energy-prev_potential_energy)<=ms_de) exit
			prev_potential_energy = potential_energy
		endif
		if (integrator_name=='nvt') call integrate_nose_hoover_chain(nhc,atoms,groups(termo_atoms_group_num),dt)
		call integrate_verlet_velocities(atoms,groups(moving_atoms_group_num),dt)
		call integrate_verlet_positions(atoms,groups(moving_atoms_group_num),dt,cell)
	endif

	call update_interactions_neibour_lists(md_step,interactions,atoms,groups,cell)
	
	call zero_forces(atoms,groups(all_atoms_group_num))
	call calculate_forces(atoms,interactions)
	call calculate_forces_numerically(atoms,interactions)
	
	if (md_step/=0) then
		call integrate_verlet_velocities(atoms,groups(moving_atoms_group_num),dt)
		if (integrator_name=='nvt') call integrate_nose_hoover_chain(nhc,atoms,groups(termo_atoms_group_num),dt)
		if (integrator_name=='nvms') call molecular_static_velocities(atoms,groups(moving_atoms_group_num))
		dt%simulation_time = dt%simulation_time+dt%ts(1)
	endif

	if (mod(md_step,integrators(integrator_index)%period_log)==0 .or. mod(md_step,out_period)==0 &
	.or. mod(md_step,out_period)==0 .or. integrator_name=='nvms') then

		call calculate_potential_energies(interactions)
		potential_energy = sum(interactions%energy)
		call calculate_temperature(temperature,kinetic_energy,atoms,groups(moving_atoms_group_num))
		total_energy = potential_energy+kinetic_energy
		if (integrator_name=='nvt') call calculate_nose_hoover_chain_energy(nhc)
		conserved_energy = total_energy+nhc%e
		
		if (mod(md_step,integrators(integrator_index)%period_log)==0) then
			write(log_id,'(A6,i9,5f24.6)',advance='no') trim(integrator_name),md_step,dt%simulation_time,&
			total_energy,potential_energy,kinetic_energy,temperature
			write(log_id,interactions_energies_format) interactions%energy
		endif
		
		if (mod(md_step,out_period)==0) then
			call calculate_force_sum(fs,atoms,groups(all_atoms_group_num))
			write(out_id,'(f12.2,A,A6,i12,f27.9,e14.4)') omp_get_wtime()-starttime,' ',&
			trim(integrator_name),md_step,conserved_energy,sqrt(sum(fs**2))
			call nlists_load(out_id,interactions)
		endif
		
	endif

	if (mod(md_step,integrators(integrator_index)%period_snapshot)==0 .and. md_step/=0) then
		write(filename,'(A,A,i6.6,A)') trim(output_prefix),'snapshot_',md_step,'.xyz'
		call write_particle_group(filename,atoms,groups(all_atoms_group_num),cell)
	endif

enddo
mdtime = omp_get_wtime()-mdtime 
if (integrator_name=='nvms') write(out_id,*) 'potential_energy_difference:', potential_energy-prev_potential_energy
write(out_id,*) 'steps number:', md_step
write(out_id,*) 'MD FINISH TIME: ',omp_get_wtime()-starttime,' SECONDS'
write(out_id,*) 'TIME STEPS PER HOUR: ',md_step/mdtime*3600

write(all_out_id,'(A32,i9,5f20.9)',advance='no') trim(output_prefix),md_step-1,dt%simulation_time,&
total_energy,potential_energy,kinetic_energy,temperature
write(all_out_id,interactions_energies_format,advance='no') interactions%energy

write(filename,'(A,A,A)') trim(output_prefix),'final_',trim(init_xyz_filename)
if (md_step/=0) call write_particle_group(filename,atoms,groups(all_atoms_group_num),cell)
close(log_id)

end subroutine md
end module md_simulation