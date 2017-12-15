program run_gr_moire_fitting
use fit_ljc_gr_moire
implicit none

real,parameter			:: gold = (sqrt(5.)-1.)/(sqrt(5.)+1.)
real 					:: init_min_ljp(3),init_max_ljp(3),ljp(3),min_ljp(3),max_ljp(3)
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

	call set_fitting_parameters(fitting_parameters_file_name,init_min_ljp,init_max_ljp)
	
	open(oid,file=trim(out_path)//trim(output_prefix)//'fit_out.txt')
	write(out_id,'(A)') line
	
	min_ljp = init_min_ljp
	max_ljp = init_max_ljp
	ljp = (init_min_ljp+init_max_ljp)/2
	
	do while(abs(prev_error_fit-error)>delta_error_fit .or. sim_num==0)
		do k=1,3		
			prev_error_fit = error
			if(sim_num==0) prev_error_fit = 1000000.
			min_ljp(k) = init_min_ljp(k)
			max_ljp(k) = init_max_ljp(k)
			
			ljp(k) = min_ljp(k)
			call calc_error(error_array(1),.true.,ljp)
			ljp(k) = max_ljp(k)
			call calc_error(error_array(4),.true.,ljp)	
			ljp(k) = min_ljp(k)+(max_ljp(k)-min_ljp(k))*gold
			call calc_error(error_array(2),.true.,ljp)
			ljp(k) = max_ljp(k)-(max_ljp(k)-min_ljp(k))*gold
			call calc_error(error_array(3),.true.,ljp)
			de = delta_error(error_array)
			write(out_id,*) ' delta_error: ',de
			gold_i = 0
			do while(de>delta_error_gold .and. gold_i<=20)
				gold_i = gold_i+1
				!if(error_array(2)<error_array(1) .and. error_array(3)<error_array(4)) then
					if(error_array(2)<error_array(3)) then
						error_array(4) = error_array(3)
						error_array(3) = error_array(2)
						max_ljp(k) = max_ljp(k)-(max_ljp(k)-min_ljp(k))*gold
						ljp(k) = min_ljp(k)+(max_ljp(k)-min_ljp(k))*gold
						call calc_error(error_array(2),.false.,ljp)
						!call calc_error(error_array(2),.true.,ljp)
					else
						error_array(1) = error_array(2)
						error_array(2) = error_array(3)
						min_ljp(k) = min_ljp(k)+(max_ljp(k)-min_ljp(k))*gold
						ljp(k) = max_ljp(k)-(max_ljp(k)-min_ljp(k))*gold
						call calc_error(error_array(3),.false.,ljp)
						!call calc_error(error_array(3),.true.,ljp)
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
			case(2);	ljp(k) = min_ljp(k)+(max_ljp(k)-min_ljp(k))*gold
			case(3);	ljp(k) = max_ljp(k)-(max_ljp(k)-min_ljp(k))*gold
			case(1);	ljp(k) = min_ljp(k)
			case(4);	ljp(k) = max_ljp(k)
			end select
			
			write(out_id,*) ' ljp: ',ljp
		enddo
	enddo

	close(oid)
	write(out_id,'(A)') line
end program run_gr_moire_fitting

function delta_error(error_array)
	real:: delta_error,me(2),error_array(4)
	!me = 1000000000.
	!do i=1,size(error_array)
	!	if(error_array(i)<me(2)) then
	!		me(1) = me(2)
	!		me(2) = error_array(i)
	!	elseif(error_array(i)<me(1)) then
	!		me(1) = error_array(i)
	!	endif
	!enddo
	!delta_error = abs(me(1)-me(2))
	delta_error = maxval(error_array)-minval(error_array)
	return
end function delta_error