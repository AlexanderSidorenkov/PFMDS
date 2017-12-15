program run_gr_moire_fitting
use fit_gr_moire
implicit none

real,parameter			:: gold = (sqrt(5.)-1.)/(sqrt(5.)+1.)
real 					:: init_min_params(4),init_max_params(4),params(4),min_params(4),max_params(4)
real 					:: error,de,prev_error_fit,delta_error_gold,delta_error_fit,error_array(4),delta_error
integer					:: i,k,gold_i,omp_get_max_threads
character(len=256)		:: arg,str,fitting_parameters_file_name

	out_id = 6
	!out_id = 66!
	!open(out_id,file='out.log')!
	out_period = 1
	num_of_omp_treads = omp_get_max_threads()
	i = 0
	do while (i<=command_argument_count())
		i = i+1; call get_command_argument(i,arg)
		select case (arg)
		case('-op','--out_period')
			i=i+1;call get_command_argument(i,str);	read(str,*) out_period
			write(out_id,'(A,A)')			'out_period: ',trim(str)
		case('-omp_n','--openmp_threads_num')
			i=i+1;call get_command_argument(i,str);	read(str,*) num_of_omp_treads
			write(out_id,'(A,A)')			'openmp_threads_num: ',trim(str)
		case('-fpfn')
			i=i+1;call get_command_argument(i,fitting_parameters_file_name)
			write(out_id,'(A,A)')			'fitting_parameters_file_name: ',trim(fitting_parameters_file_name)
		case('-delta_error_gold')
			i=i+1;call get_command_argument(i,str);	read(str,*) delta_error_gold
			write(out_id,'(A,A)')			'delta_error_gold: ',trim(str)
		case('-delta_error_fit')
			i=i+1;call get_command_argument(i,str);	read(str,*) delta_error_fit
			write(out_id,'(A,A)')			'delta_error_fit: ',trim(str)
		end select
	enddo

	call set_fitting_parameters(fitting_parameters_file_name,init_min_params,init_max_params)
	
	open(oid,file=trim(out_path)//trim(output_prefix)//'fit_out.txt')
	write(out_id,'(A)') line
	
	min_params = init_min_params
	max_params = init_max_params
	params = (init_min_params+init_max_params)/2
	
	do while(abs(prev_error_fit-error)>delta_error_fit .or. sim_num==0)
		do k=1,3		
			prev_error_fit = error
			if(sim_num==0) prev_error_fit = 1000000.
			min_params(k) = init_min_params(k)
			max_params(k) = init_max_params(k)
			
			params(k) = min_params(k)
			call calc_error(error_array(1),.true.,params)
			params(k) = max_params(k)
			call calc_error(error_array(4),.true.,params)	
			params(k) = min_params(k)+(max_params(k)-min_params(k))*gold
			call calc_error(error_array(2),.true.,params)
			params(k) = max_params(k)-(max_params(k)-min_params(k))*gold
			call calc_error(error_array(3),.true.,params)
			de = delta_error(error_array)
			write(out_id,*) ' delta_error: ',de
			gold_i = 0
			do while(de>delta_error_gold .and. gold_i<=20)
				gold_i = gold_i+1
				!if(error_array(2)<error_array(1) .and. error_array(3)<error_array(4)) then
					if(error_array(2)<error_array(3)) then
						error_array(4) = error_array(3)
						error_array(3) = error_array(2)
						max_params(k) = max_params(k)-(max_params(k)-min_params(k))*gold
						params(k) = min_params(k)+(max_params(k)-min_params(k))*gold
						call calc_error(error_array(2),.false.,params)
						!call calc_error(error_array(2),.true.,params)
					else
						error_array(1) = error_array(2)
						error_array(2) = error_array(3)
						min_params(k) = min_params(k)+(max_params(k)-min_params(k))*gold
						params(k) = max_params(k)-(max_params(k)-min_params(k))*gold
						call calc_error(error_array(3),.false.,params)
						!call calc_error(error_array(3),.true.,params)
					endif
				!else
				!	write(out_id,*) ' no monotonicity ',error_array
					!stop
				!	exit
				!endif
				de = delta_error(error_array)
				write(out_id,*) ' delta_error: ',de
			enddo
			error = minval(error_array)
			select case (minloc(error_array,1))
			case(2);	params(k) = min_params(k)+(max_params(k)-min_params(k))*gold
			case(3);	params(k) = max_params(k)-(max_params(k)-min_params(k))*gold
			case(1);	params(k) = min_params(k)
			case(4);	params(k) = max_params(k)
			end select
			
			write(out_id,*) ' parameters: ',params
		enddo
	enddo

	close(oid)
	write(out_id,'(A)') line
end program run_gr_moire_fitting

function delta_error(error_array)
	real:: delta_error,error_array(4)
	delta_error = maxval(error_array)-minval(error_array)
	return
end function delta_error