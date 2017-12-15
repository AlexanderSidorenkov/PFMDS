module fit_gr_moire
use md_simulation
use graphene_on_surface_analysis
implicit none

integer					:: sim_num,out_period,num_of_omp_treads,out_id,final_out_id,oid
integer					:: ar_c_num(2)
character(len=256)		:: interaction_name,ar_settings_filename(2),output_prefix,input_path,out_path,&
						   ar_final_file(2),param_file,ar_start_xyz_file(2),ar_xyz_file(2)						
real 					:: z,ar_zero_energy_level(2),be0,ar_grd0(2),Rcut(2)
logical					:: simplified
character(len=80)		:: line

contains

subroutine calc_error(error,from_init_xyz,params)
character(len=256)		::  settings_filename,start_xyz_file,xyz_file,final_file,str,tempstring,op,prevop
real 					::	params(4),error,arr1(3),arr2(3),be,bd,grd
logical					::	from_init_xyz
integer					::	i,rnm_err

	error = 0
	sim_num = sim_num+1
	write(out_id,*)
	write(out_id,*) trim(input_path)//trim(param_file)

	if(interaction_name=='ljc') then
		write(out_id,*) sim_num,params(3),params(1),params(2)
		write(oid,'(i6,3f21.6)',advance='no') sim_num,params(3),params(1),params(2)
		open(1234,file=trim(input_path)//trim(param_file))
		write(1234,*) params(3),params(1),params(2)
		write(1234,*) Rcut(1),Rcut(2)
		write(1234,*) simplified
		close(1234)
	endif
	
	if(interaction_name=='morsec') then
		write(out_id,*) sim_num,params(3),params(1),params(4),params(2)
		write(oid,'(i6,4f21.6)',advance='no') sim_num,params(3),params(1),params(4),params(2)
		open(1234,file=trim(input_path)//trim(param_file))
		write(1234,*) params(3),params(1),params(4),params(2)
		write(1234,*) Rcut(1),Rcut(2)
		write(1234,*) simplified
		close(1234)
	endif
	
	do i=1,2
	
		settings_filename = trim(ar_settings_filename(i))!;	write(out_id,*) settings_filename
		start_xyz_file = trim(ar_start_xyz_file(i))!;		write(out_id,*) start_xyz_file
		xyz_file = trim(ar_xyz_file(i))!;					write(out_id,*) xyz_file
		final_file = trim(ar_final_file(i))!;				write(out_id,*) final_file
		
		if(from_init_xyz) then
			rnm_err = rename(trim(input_path)//trim(start_xyz_file),trim(input_path)//trim(xyz_file))
		else
			write(str,'(i6.6)') sim_num-1
			write(prevop,'(A)') trim(out_path)//trim(output_prefix)//trim(str)//'_'
			rnm_err = rename(trim(prevop)//'final_'//trim(xyz_file),trim(input_path)//trim(xyz_file))
		endif
		
		write(str,'(i6.6)') sim_num
		write(op,'(A)') trim(out_path)//trim(output_prefix)//trim(str)//'_'
		open(final_out_id,file=trim(op)//trim(final_file))
		write(out_id,'(A)') line 
		call md(out_id,final_out_id,input_path,settings_filename,op,out_period,num_of_omp_treads)
		write(out_id,'(A)') line
		write(tempstring,'(A)') trim(op)//'final_'//trim(xyz_file)
		call gr_on_cu_analysis(arr1,arr2,tempstring,z)
		write(out_id,*) 'gr_on_cu_analysis:',arr1,arr2
		write(final_out_id,'(3f16.6)') arr1-arr2(1)
		close(final_out_id)

		open(final_out_id,file=trim(op)//trim(final_file))
		read(final_out_id,'(A61,f20.9,A)') tempstring,be,tempstring
		close(final_out_id)

		if(from_init_xyz) then
			rnm_err = rename(trim(input_path)//trim(xyz_file),trim(input_path)//trim(start_xyz_file))
		else
			rnm_err = rename(trim(input_path)//trim(xyz_file),trim(prevop)//'final_'//trim(xyz_file))
		endif
		
		bd = arr1(1)-arr2(1)
		grd = arr1(3)-arr1(2)
		be = (be-ar_zero_energy_level(i))/ar_c_num(i)
		
		write(oid,'(3f21.6)',advance='no') be,bd,grd
		
		if(i==1) then
			error = (be/be0-1.)**2+(grd/ar_grd0(i)-1.)**2
			write(oid,'(2f21.6)',advance='no') (be/be0-1.)**2,(grd/ar_grd0(i)-1.)**2
		endif
		if(i==2) then
			error = error+(grd/ar_grd0(i)-1.)**2
			write(oid,'(f21.6)',advance='no') (grd/ar_grd0(i)-1.)**2
		endif
		
	enddo
	
	write(oid,'(3f21.6)') error
	rnm_err = rename(trim(input_path)//trim(param_file),trim(input_path)//trim(str)//trim(param_file))	
	
end subroutine calc_error

subroutine set_fitting_parameters(fitting_parameters_file_name,init_min_params,init_max_params)
real 					:: init_min_params(4),init_max_params(4)
character(len=256)		:: fitting_parameters_file_name,min_param_file,max_param_file,str
integer					:: i

	final_out_id = 1211
	oid = 111
	do i=1,80
		line(i:i) = '_'
	enddo

	open(9,file=fitting_parameters_file_name)
	read(9,'(A16,A)') str,input_path;		write(out_id,*) trim(str),'/t',trim(input_path)
	read(9,'(A16,A)') str,out_path;			write(out_id,*) trim(str),'/t',trim(out_path)
	read(9,*) str,interaction_name;			write(out_id,*) trim(str),'/t',trim(interaction_name)
	read(9,*) str,output_prefix;			write(out_id,*) trim(str),'/t',trim(output_prefix)
	read(9,*) str,ar_settings_filename(1);	write(out_id,*) trim(str),'/t',trim(ar_settings_filename(1))
	read(9,*) str,ar_start_xyz_file(1);		write(out_id,*) trim(str),'/t',trim(ar_start_xyz_file(1))
	read(9,*) str,ar_c_num(1);				write(out_id,*) trim(str),'/t',ar_c_num(1)
	read(9,*) str,ar_zero_energy_level(1);	write(out_id,*) trim(str),'/t',ar_zero_energy_level(1)
	read(9,*) str,ar_final_file(1);			write(out_id,*) trim(str),'/t',trim(ar_final_file(1))
	read(9,*) str,ar_settings_filename(2);	write(out_id,*) trim(str),'/t',trim(ar_settings_filename(2))
	read(9,*) str,ar_start_xyz_file(2);		write(out_id,*) trim(str),'/t',trim(ar_start_xyz_file(2))
	read(9,*) str,ar_c_num(2);				write(out_id,*) trim(str),'/t',ar_c_num(2)
	read(9,*) str,ar_zero_energy_level(2);	write(out_id,*) trim(str),'/t',ar_zero_energy_level(2)
	read(9,*) str,ar_final_file(2);			write(out_id,*) trim(str),'/t',trim(ar_final_file(2))
	read(9,*) str,min_param_file;			write(out_id,*) trim(str),'/t',trim(min_param_file)
	read(9,*) str,max_param_file;			write(out_id,*) trim(str),'/t',trim(max_param_file)
	read(9,*) str,z;						write(out_id,*) trim(str),'/t',z
	read(9,*) str,be0;						write(out_id,*) trim(str),'/t',be0
	read(9,*) str,ar_grd0(1);				write(out_id,*) trim(str),'/t',ar_grd0(1)
	read(9,*) str,ar_grd0(2);				write(out_id,*) trim(str),'/t',ar_grd0(2)
	close(9)
	

	simplified = .false.
	if(interaction_name=='ljc') then
		open(1234,file=trim(input_path)//trim(min_param_file))
		read(1234,*) init_min_params(3),init_min_params(1),init_min_params(2)
		read(1234,*) Rcut(1),Rcut(2)
		close(1234)
		open(1234,file=trim(input_path)//trim(max_param_file))
		read(1234,*) init_max_params(3),init_max_params(1),init_max_params(2)
		close(1234)
		init_min_params(4) = 0.
		init_max_params(4) = 0.
	endif
	if(interaction_name=='morsec') then
		open(1234,file=trim(input_path)//trim(min_param_file))
		read(1234,*) init_min_params(3),init_min_params(1),init_min_params(4),init_min_params(2)
		read(1234,*) Rcut(1),Rcut(2)
		close(1234)
		open(1234,file=trim(input_path)//trim(max_param_file))
		read(1234,*) init_max_params(3),init_max_params(1),init_min_params(4),init_max_params(2)
		close(1234)
	endif
	
	sim_num = 0
	
	open(1234,file=trim(input_path)//trim(ar_settings_filename(1)))
	read(1234,*);read(1234,*);read(1234,*) str,ar_xyz_file(1)
	do while(.true.)
		read(1234,'(A128)') str
		if(str(1:3)==trim(interaction_name)) then
			read(str(4:),*) param_file
			exit
		endif
		if(str(1:6)==trim(interaction_name)) then
			read(str(7:),*) param_file
			exit
		endif
	enddo
	close(1234)
	open(1234,file=trim(input_path)//trim(ar_settings_filename(2)))
	read(1234,*);read(1234,*);read(1234,*) str,ar_xyz_file(2)
	close(1234)

end subroutine set_fitting_parameters

end module fit_gr_moire