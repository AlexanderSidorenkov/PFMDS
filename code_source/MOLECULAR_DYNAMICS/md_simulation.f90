module md_simulation
use perfomance_settings
use md_general
use md_integrators
use md_read_write
use md_interactions
implicit none 
contains

!> Молекулярная динамика.
!> \param[in] out_id Идентификатор основного потока вывода
!> \param[in] all_out_id Идентификатор потока вывода для конечных значений энергий, температуры и прочего
!> \param[in] input_path Путь к папке с входными файлами
!> \param[in] settings_filename Имя файла настроек
!> \param[in] output_prefix	Префикс выходных файлов
!> \param[in] out_period Период основного вывода в шагах МД
!> \param[in] num_of_omp_treads Количество OpenMP потоков
!> \param[in] rand_seed Число для инициализации ГСЧ
subroutine md(out_id,all_out_id,input_path,settings_filename,output_prefix,out_period,num_of_omp_treads,rand_seed)

type(simulation_cell)					:: cell
type(time_steps)						:: dt
type(particles)							:: atoms
type(particle_group),allocatable		:: groups(:)
type(interaction),allocatable			:: interactions(:)
type(integrator_params),allocatable		:: integrators(:)
type(nose_hoover_chain),allocatable		:: nhc_thermostats(:)
real									::	exe_t,exe_time_start,exe_time_md,exe_time_pos_vel,&
											exe_time_nlists,exe_time_nlsearch,exe_time_nldistance,exe_time_forces,exe_time_energy,&
											conserved_energy,nose_hoover_energy,total_energy,kinetic_energy,potential_energy,prev_potential_energy,&
											ms_de,fs(3),mcv(3),mc(3),mav_vel,initial_temperature,nhc_temperature,temperature,nhc_q1
integer									::	i,rand_seed,num_of_omp_treads,md_step,md_step_limit,integrators_num,integrator_index,&
											file_id,out_id,log_id,out_period,all_out_id,&
											all_atoms_group_num,nhc_atoms_group_num,nhc_M,nhc_num,&
											all_moving_atoms_group_num,xyz_moving_atoms_group_num,z_moving_atoms_group_num,&
											traj_group_num,period_traj,change_group_num,&
											zero_momentum_period
integer,allocatable						::	group_change_from(:),group_change_to(:),change_ts1(:),change_ts2(:),change_frec(:)
character(len=128)						::	str,filename,logfilename,init_xyz_filename,input_path,settings_filename,output_prefix
character(len=32)						::	integrator_name,interactions_energies_format,nose_hoover_energies_format
logical 								::	new_velocities,invert_z_vel
real									::	omp_get_wtime

exe_time_start = omp_get_wtime()
log_id = 108
file_id = 2017

open(file_id,file=trim(input_path)//settings_filename);	write(out_id,'(A,A)') 'settings_filename: ',trim(settings_filename)
read(file_id,*) str,md_step_limit;						write(out_id,'(A32,i12)') str,md_step_limit
read(file_id,*) str,logfilename;						write(out_id,'(A32,A,A)') str,'	',trim(logfilename)
read(file_id,*) str,init_xyz_filename;					write(out_id,'(A32,A,A)') str,'	',trim(init_xyz_filename)
call read_box_size(cell,trim(input_path)//init_xyz_filename);	write(out_id,'(A,3f16.6)') 	'box_size: ',cell%box_size
call read_particles(atoms,trim(input_path)//init_xyz_filename);	write(out_id,'(A,i12)') 	'particles_num: ',atoms%N
read(file_id,*) str,new_velocities;						write(out_id,'(A32,l8)') str,new_velocities
read(file_id,*) str,zero_momentum_period;				write(out_id,'(A32,i12)') str,zero_momentum_period
call create_groups(groups,file_id,out_id,atoms)
read(file_id,*) str,all_moving_atoms_group_num;			write(out_id,'(A32,i12)') str,all_moving_atoms_group_num
read(file_id,*) str,xyz_moving_atoms_group_num;			write(out_id,'(A32,i12)') str,xyz_moving_atoms_group_num
read(file_id,*) str,z_moving_atoms_group_num;			write(out_id,'(A32,i12)') str,z_moving_atoms_group_num
read(file_id,*) str,all_atoms_group_num;				write(out_id,'(A32,i12)') str,all_atoms_group_num
read(file_id,*) str,traj_group_num;						write(out_id,'(A32,i12)') str,traj_group_num
read(file_id,*) str,period_traj;						write(out_id,'(A32,i12)') str,period_traj
read(file_id,*) str,change_group_num;					write(out_id,'(A32,i12)') str,change_group_num
allocate(group_change_from(change_group_num),group_change_to(change_group_num),&
change_ts1(change_group_num),change_ts2(change_group_num),change_frec(change_group_num))
do i=1,change_group_num
	read(file_id,*) str,group_change_from(i),group_change_to(i);
	write(out_id,'(A32,2i12)') str,group_change_from(i),group_change_to(i)
	read(file_id,*) str,change_ts1(i),change_ts2(i),change_frec(i)
	write(out_id,'(A32,3i12)') str,change_ts1(i),change_ts2(i),change_frec(i)
enddo
read(file_id,*) str,invert_z_vel;						write(out_id,'(A32,l8)') str,invert_z_vel
read(file_id,*) str,integrators_num;					write(out_id,'(A32,i12)') str,integrators_num
allocate(integrators(0:integrators_num)); integrators(0)%int_name='none'; integrators(0)%dt=0.; integrators(0)%l=0
read(file_id,'(A)') str;								write(out_id,'(A)') str
do i=1,integrators_num
	call read_integrator_params(integrators(i),file_id)
	write(out_id,'(A,A6,f10.5,4i9)') '  ',integrators(i)%int_name,&
	integrators(i)%dt,integrators(i)%l,integrators(i)%period_snapshot,integrators(i)%period_log
enddo
read(file_id,*) str,ms_de;								;write(out_id,'(A32,es16.6)') str,ms_de
read(file_id,*) str,nhc_num								;write(out_id,'(A32,i12)') str,nhc_num
allocate(nhc_thermostats(nhc_num))
do i=1,nhc_num
	read(file_id,*) nhc_atoms_group_num,nhc_temperature,nhc_M,nhc_q1
	write(out_id,'(i6,f16.6,i6,f16.6)') nhc_atoms_group_num,nhc_temperature,nhc_M,nhc_q1
	call create_nose_hoover_chain(nhc_thermostats(i),nhc_M)
	call set_nose_hoover_chain(nhc_thermostats(i),nhc_temperature,nhc_q1,nhc_atoms_group_num,groups(nhc_atoms_group_num)%N)
enddo
read(file_id,*) str,initial_temperature;				write(out_id,'(A32,f16.6)') str,initial_temperature
if (initial_temperature<0.) initial_temperature=0.
if (new_velocities) call set_new_temperature(atoms,groups(all_moving_atoms_group_num),initial_temperature,rand_seed)
call create_interactions(interactions,groups,file_id,out_id,input_path)
close(file_id)

write(interactions_energies_format,'("(",i0,"f20.9)")') size(interactions)
write(nose_hoover_energies_format,'("(",i0,"f20.9)")') size(nhc_thermostats)
write(str,'(A,A)') trim(output_prefix),trim(logfilename)
open(log_id,file=str)
write(out_id,'(A24,1f10.2,A)') 'PREPARATIONS TIME: ',omp_get_wtime()-exe_time_start,' S '
write(out_id,'(A24,i6,A)') 'RUNNING ON ',num_of_omp_treads,' OPENMP THREADS'
write(out_id,*)
call set_openmp_perfomance(num_of_omp_treads,atoms%N)
exe_time_md = omp_get_wtime()
exe_time_pos_vel = 0.
exe_time_nlists = 0.
exe_time_nlsearch = 0.
exe_time_nldistance = 0.
exe_time_forces = 0.
exe_time_energy = 0.
integrator_index = 0
dt%simulation_time = 0.

do md_step=0,md_step_limit

	do i=1,change_group_num
		call change_particle_group_N(groups(group_change_to(i)),md_step,change_ts1(i),change_ts2(i),&
		change_frec(i),groups(group_change_from(i)))
	enddo
	
	if (md_step-1==sum(integrators(0:integrator_index)%l) .or. integrator_index==0) then
		integrator_index = integrator_index+1
		do i=integrator_index,integrators_num
			if (integrators(i)%l>0) then
				integrator_name = integrators(i)%int_name
				call init_time_steps(dt,integrators(i)%dt)
				!if (integrator_name=='nvt') call set_nose_hoover_chain(nhc,nhc_temperature,nhc_q1,groups(nhc_atoms_group_num)%N)
				exit
			endif
		enddo
		if (i>integrators_num .or. (i==integrators_num .and. integrators(integrators_num)%l<1)) then
			write(out_id,*) 'no more integrators'
			exit
		endif
		integrator_index = i
	endif

	exe_t = omp_get_wtime()
	call check_positions(out_id,atoms,cell)
	if (invert_z_vel) call invert_z_velocities(atoms,0.8*cell%box_size(3),0.9*cell%box_size(3))
	if (md_step/=0) then
		if (integrator_name=='nvms') then
			if (abs(potential_energy-prev_potential_energy)<=ms_de) then
				write(out_id,*) 'potential energy diffrence is small enough'
				exit
			endif
		endif
		if (integrator_name=='nvt') then
			do i=1,nhc_num
				call integrate_nose_hoover_chain(nhc_thermostats(i),atoms,groups(nhc_thermostats(i)%group_num),dt)
			enddo
		endif
		call integrate_verlet_xyz_velocities(atoms,groups(xyz_moving_atoms_group_num),dt)
		call integrate_verlet_z_velocities(atoms,groups(z_moving_atoms_group_num),dt)
		call integrate_verlet_xyz_positions(atoms,groups(xyz_moving_atoms_group_num),dt,cell)
		call integrate_verlet_z_positions(atoms,groups(z_moving_atoms_group_num),dt,cell)
	endif
	prev_potential_energy = potential_energy
	exe_time_pos_vel = exe_time_pos_vel+omp_get_wtime()-exe_t

	exe_t = omp_get_wtime()
	call update_interactions_neighbour_lists(md_step,interactions,atoms,groups,cell,exe_time_nlsearch,exe_time_nldistance)
	exe_time_nlists = exe_time_nlists+omp_get_wtime()-exe_t
	exe_t = omp_get_wtime()
	if (mod(md_step,zero_momentum_period)==0) call zero_momentum(atoms,groups(all_atoms_group_num))
	call zero_forces(atoms,groups(all_atoms_group_num))
	call calculate_forces(atoms,interactions)
	call calculate_forces_numerically(atoms,interactions)
	exe_time_forces = exe_time_forces+omp_get_wtime()-exe_t

	exe_t = omp_get_wtime()
	if (md_step/=0) then
		call integrate_verlet_xyz_velocities(atoms,groups(xyz_moving_atoms_group_num),dt)
		call integrate_verlet_z_velocities(atoms,groups(z_moving_atoms_group_num),dt)
		if (integrator_name=='nvt') then
			do i=1,nhc_num
				call integrate_nose_hoover_chain(nhc_thermostats(i),atoms,groups(nhc_thermostats(i)%group_num),dt)
			enddo
		endif
		if (integrator_name=='nvms') then
			call molecular_static_xyz_velocities(atoms,groups(xyz_moving_atoms_group_num))
			call molecular_static_1D_velocities(atoms,groups(z_moving_atoms_group_num))
		endif
		dt%simulation_time = dt%simulation_time+dt%ts(1)
	endif
	exe_time_pos_vel = exe_time_pos_vel+omp_get_wtime()-exe_t
	
	if (mod(md_step,integrators(integrator_index)%period_log)==0 .or. mod(md_step,out_period)==0 .or. integrator_name=='nvms') then

		exe_t = omp_get_wtime()
		call calculate_potential_energies(interactions)
		potential_energy = sum(interactions%energy)
		call calculate_temperature(temperature,kinetic_energy,atoms,groups(all_moving_atoms_group_num))
		if (integrator_name=='nvt') then
			do i=1,nhc_num
				call calculate_nose_hoover_chain_energy(nhc_thermostats(i))
			enddo
		endif
		nose_hoover_energy = sum(nhc_thermostats%e)
		total_energy = potential_energy+kinetic_energy
		conserved_energy = total_energy+nose_hoover_energy
		exe_time_energy = exe_time_energy+omp_get_wtime()-exe_t
		
		if (mod(md_step,integrators(integrator_index)%period_log)==0) then
			write(log_id,'(A6,i9,7f24.6)',advance='no') trim(integrator_name),md_step,dt%simulation_time,&
			conserved_energy,nose_hoover_energy,total_energy,potential_energy,kinetic_energy,temperature
			write(log_id,interactions_energies_format,advance='no') interactions%energy
			write(log_id,nose_hoover_energies_format) nhc_thermostats%e
		endif
		
		if (mod(md_step,out_period)==0) then
			call calculate_force_sum(fs,atoms,groups(all_atoms_group_num))
			call calculate_mass_center(mc,atoms,groups(all_atoms_group_num))
			call calculate_mass_center_velocity(mcv,atoms,groups(all_atoms_group_num))
			write(out_id,'(A,f10.2)') 'exe time (s) = ',omp_get_wtime()-exe_time_start
			write(out_id,'(A,i12,A,A6)') 'step = ',md_step,' integrator = ',trim(integrator_name)
			write(out_id,'(A,es21.9)') 'conserved energy (eV) = ',conserved_energy
			write(out_id,'(A,es16.6)') 'forces sum (eV/A) = ',sqrt(sum(fs**2))
			write(out_id,'(A,es16.6)') 'c.o.m. velocity (A/fs) = ',sqrt(sum(mcv**2))
			write(out_id,'(A,3f12.6)') 'c.o.m. position (A) = ',mc
			call find_max_velocity(mav_vel,atoms)
			write(out_id,'(A,f12.6)') 'maximum velocity (A/fs) = ',mav_vel
			call nlists_load(out_id,interactions)
			if (integrator_name/='nvms' .and. sqrt(sum(fs**2))>10.**(-14)) &
			write(out_id,*) 'WARNING: forces sum is too big'
			if (integrator_name=='nvt' .and. xyz_moving_atoms_group_num==all_atoms_group_num .and. sqrt(sum(mcv**2))>10.**(-14)) &
			write(out_id,*) 'WARNING: center of mass velocity is too big'
			write(out_id,*)
		endif
		
	endif

	if (mod(md_step,integrators(integrator_index)%period_snapshot)==0 .and. md_step/=0) then
		write(filename,'(A,A,i6.6,A)') trim(output_prefix),'snapshot_',md_step,'.xyz'
		call write_particle_group(filename,atoms,groups(all_atoms_group_num),cell)
	endif

	if (mod(md_step,period_traj)==0) then
		write(filename,'(A,A,i2.2,A)') trim(output_prefix),'traj_',traj_group_num,'.xyz'
		call write_particle_group_append(filename,atoms,groups(traj_group_num),cell,md_step)
	endif
	
enddo

exe_time_md = omp_get_wtime()-exe_time_md 
if (integrator_name=='nvms') write(out_id,*) 'potential energy difference: ',potential_energy-prev_potential_energy
write(out_id,*) 'steps number:', md_step-1

write(out_id,*)
write(out_id,*) 'PERFOMANCE:'
write(out_id,'(A24,f10.2,A,f10.2,A)') 'MD:',exe_time_md,' S ',exe_time_md/exe_time_md*100,'%'
write(out_id,'(A24,f10.2,A,f10.2,A)') 'POSITION AND VELOCITY:',exe_time_pos_vel,' S ',exe_time_pos_vel/exe_time_md*100,'%'
write(out_id,'(A24,f10.2,A,f10.2,A)') 'NEIGHBOURS:',exe_time_nlists,' S ',exe_time_nlists/exe_time_md*100,'%'
write(out_id,'(A24,f10.2,A,f10.2,A)') 'NEIGHBOURS SEARCH:',exe_time_nlsearch,' S ',exe_time_nlsearch/exe_time_md*100,'%'
write(out_id,'(A24,f10.2,A,f10.2,A)') 'NEIGHBOURS DISTANCE:',exe_time_nldistance,' S ',exe_time_nldistance/exe_time_md*100,'%'
write(out_id,'(A24,f10.2,A,f10.2,A)') 'FORCES:',exe_time_forces,' S ',exe_time_forces/exe_time_md*100,'%'
write(out_id,'(A24,f10.2,A,f10.2,A)') 'ENERGY:',exe_time_energy,' S ',exe_time_energy/exe_time_md*100,'%'
write(out_id,'(A24,f10.2,A,f10.2,A)') 'REST:',exe_time_md-exe_time_pos_vel-exe_time_nlists-exe_time_forces-exe_time_energy,&
' S ',(exe_time_md-exe_time_pos_vel-exe_time_nlists-exe_time_forces-exe_time_energy)/exe_time_md*100,'%'
write(out_id,'(A24,f16.2)') 'TIME STEPS PER HOUR:',md_step/exe_time_md*3600
write(out_id,*)

write(all_out_id,'(A32,i9,5f20.9)',advance='no') trim(output_prefix),md_step-1,dt%simulation_time,&
total_energy,potential_energy,kinetic_energy,temperature
write(all_out_id,interactions_energies_format,advance='no') interactions%energy
write(log_id,'(A32,i9,5f20.9)',advance='no') trim(output_prefix),md_step-1,dt%simulation_time,&
total_energy,potential_energy,kinetic_energy,temperature
write(log_id,interactions_energies_format,advance='no') interactions%energy

write(filename,'(A,A,A)') trim(output_prefix),'final_',trim(init_xyz_filename)
if (md_step/=0) call write_particle_group(filename,atoms,groups(all_atoms_group_num),cell)
close(log_id)

end subroutine md

end module md_simulation