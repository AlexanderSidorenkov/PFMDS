program md_run
use md_simulation
implicit none
include 'mpif.h' 

integer					::	i,out_period,num_of_omp_treads,out_id,all_out_id,settings_files_list_id,set_num
character(len=128)		::	settings_filename,settings_files_list,all_out_file,output_prefix,list_p,input_path,out_path,str
character(len=32)		::	arg,str_node_id
character(len=80)		::	line
integer					::	omp_get_max_threads
integer					::	node_id, number_of_nodes
integer					::	ierr

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,node_id,ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,number_of_nodes,ierr)
node_id = node_id+1
write(str_node_id,'(i4.4,A)') node_id,'-'

do i=1,80
	line(i:i) = '_'
enddo
out_id = 2018
all_out_id = 66
settings_files_list_id = 1234
out_period = 1
settings_files_list = ''
settings_filename = 'md_run_settings.txt'
all_out_file = 'all_out.txt'
output_prefix = ''
input_path = ''
out_path = ''
num_of_omp_treads = omp_get_max_threads()

i = 0
do while (i<=command_argument_count())
	i = i+1; call get_command_argument(i,arg)
	select case (arg)
	case('-i','--input')
		i=i+1;call get_command_argument(i,settings_filename);
	case('-p','--prefix')
		i=i+1;call get_command_argument(i,output_prefix)	
	case('-op','--out_period')
		i=i+1;call get_command_argument(i,str);read(str,*) out_period
	case('-omp_n','--openmp_threads_num')
		i=i+1;call get_command_argument(i,str);read(str,*) num_of_omp_treads
	case('-ipath','--input_path') 
		i=i+1;call get_command_argument(i,input_path)
	case('-ilist','--input_list') 
		i=i+1;call get_command_argument(i,settings_files_list)
	case('-opath','--out_path') 
		i=i+1;call get_command_argument(i,out_path)
	case('-ofile','--all_out_file') 
		i=i+1;call get_command_argument(i,all_out_file)
	end select
enddo

open(out_id,file=trim(output_prefix)//trim(str_node_id)//'out.txt')

write(out_id,'(A,A)')  'output_prefix:      ',trim(output_prefix)
write(out_id,'(A,i9)') 'out_period:         ',out_period
write(out_id,'(A,i9)') 'openmp_threads_num: ',num_of_omp_treads
write(out_id,'(A,A)')  'input_path:         ',trim(input_path)	
write(out_id,'(A,A)')  'input_list:         ',trim(settings_files_list)
write(out_id,'(A,A)')  'out_path:           ',trim(out_path)
write(out_id,'(A,A)')  'all_out_file:       ',trim(all_out_file)

write(out_id,'(A)') line
if (settings_files_list/='') then
	open(all_out_id,file=trim(out_path)//trim(output_prefix)//trim(str_node_id)//trim(all_out_file))
	
	open(settings_files_list_id,file=trim(input_path)//trim(settings_files_list))
	read(settings_files_list_id,*) set_num
	if (number_of_nodes<=set_num) then
		do i=1,set_num
			read(settings_files_list_id,*) settings_filename,list_p
			if (mod(i-1,number_of_nodes)==node_id-1) then
				write(out_id,'(A)') line
				write(out_id,'(A,i6,A,i6,A)') 'RUNNING ON NODE ',node_id,' OUT OF', number_of_nodes,' NODES'
				str = trim(out_path)//trim(output_prefix)//trim(str_node_id)//trim(list_p)
				write(out_id,'(i6,A,A,A,A)') i,'	',trim(settings_filename),'	',trim(str)
				call md(out_id,all_out_id,input_path,settings_filename,str,out_period,num_of_omp_treads)
				write(all_out_id,*)
				write(out_id,'(A)') line
			endif
		enddo
	else
		write(out_id,'(A,i6,A,i6,A,A)') 'error: too many mpi nodes (',number_of_nodes,') for this list (',set_num,')',settings_files_list
	endif
	close(settings_files_list_id)
	
	close(all_out_id)
else
	write(out_id,'(A)') line
	write(out_id,'(A,i6,A,i6,A)') 'RUNNING ON NODE ',node_id,' OUT OF', number_of_nodes,' NODES. EACH NODE RUNS THE SAME SIMULATION.'
	str = trim(out_path)//trim(output_prefix)//trim(str_node_id)//trim(list_p)
	write(out_id,'(i6,A,A,A,A)') i,'	',trim(settings_filename),'	',trim(str)
	call md(out_id,out_id,input_path,settings_filename,str,out_period,num_of_omp_treads)
	write(out_id,*)
	write(out_id,'(A)') line
endif
write(out_id,'(A)') line

close(out_id)

call MPI_FINALIZE(ierr)
end program md_run