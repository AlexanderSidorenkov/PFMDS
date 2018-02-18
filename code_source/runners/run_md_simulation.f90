program md_run
use md_simulation
implicit none

integer					::	i,out_period,num_of_omp_treads,out_id,all_out_id,set_num
character(len=128)		::	settings_filename,settings_files_list,all_out_file,output_prefix,input_path,out_path,str
character(len=32)		::	arg
character(len=80)		::	line
integer					::	omp_get_max_threads

do i=1,80
	line(i:i) = '_'
enddo
out_id = 6
all_out_id = 66
out_period = 1
settings_files_list = ''
settings_filename = 'md_run_settings.txt'
all_out_file = 'all_out.txt'
output_prefix = ''
input_path = ''
out_path = ''
num_of_omp_treads = omp_get_max_threads()
write(out_id,'(A)') line
i = 0
do while (i<=command_argument_count())
	i = i+1; call get_command_argument(i,arg)
	select case (arg)
	case('-i','--input')
		i=i+1;call get_command_argument(i,settings_filename);
	case('-p','--prefix')
		i=i+1;call get_command_argument(i,output_prefix)
		write(out_id,'(A,A)')	'output_prefix:  ',trim(output_prefix)
	case('-op','--out_period')
		i=i+1;call get_command_argument(i,str)
		write(out_id,'(A,A)')			'out_period:         ',trim(str)
		read(str,*) out_period
	case('-omp_n','--openmp_threads_num')
		i=i+1;call get_command_argument(i,str)
		write(out_id,'(A,A)')			'openmp_threads_num: ',trim(str)
		read(str,*) num_of_omp_treads
	case('-ipath','--input_path') 
		i=i+1;call get_command_argument(i,input_path)
		write(out_id,'(A,A)')			'input_path:         ',trim(input_path)
	case('-ilist','--input_list') 
		i=i+1;call get_command_argument(i,settings_files_list)
		write(out_id,'(A,A)')			'input_list:         ',trim(settings_files_list)
	case('-opath','--out_path') 
		i=i+1;call get_command_argument(i,out_path)
		write(out_id,'(A,A)')			'out_path:           ',trim(out_path)
	case('-ofile','--all_out_file') 
		i=i+1;call get_command_argument(i,all_out_file)
		write(out_id,'(A,A)')			'all_out_file:   ',trim(all_out_file)
	end select
enddo
if (settings_files_list/='') then
	open(1234,file=trim(input_path)//trim(settings_files_list))
	open(all_out_id,file=trim(out_path)//trim(all_out_file))
	read(1234,*) set_num
	do i=1,set_num
		read(1234,*) settings_filename,str
		write(output_prefix,'(A,A)') trim(out_path),trim(str)
		write(out_id,'(A)') line
		write(out_id,'(i6,A,A32,A32)') i,'	',settings_filename,output_prefix
		call md(out_id,all_out_id,input_path,settings_filename,output_prefix,out_period,num_of_omp_treads)
		write(all_out_id,*)
		write(out_id,'(A)') line
	enddo
	close(1234)
	close(all_out_id)
else
	write(out_id,'(A)') line
	call md(out_id,out_id,input_path,settings_filename,output_prefix,out_period,num_of_omp_treads)
	write(out_id,*)
	write(out_id,'(A)') line
endif
write(out_id,'(A)') line
end program md_run