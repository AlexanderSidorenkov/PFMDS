program run_gr_analysis
use graphene_on_surface_analysis
implicit none

integer									:: i,n,inid,outid
real									:: z,arr1(3),arr2(3)
character(len=256)						:: filelist,filename,path,outfilename
character(len=32)						:: arg,str

inid = 111
outid = 222
path = ''
filelist = 'filelist.txt'
outfilename = 'outfilename.txt'
i = 0
do while (i<=command_argument_count())
	i = i+1; call get_command_argument(i,arg)
	select case (arg)
	case('-path')
		i=i+1;call get_command_argument(i,path);
	case('-fl')
		i=i+1;call get_command_argument(i,filelist);
	case('-o')
		i=i+1;call get_command_argument(i,outfilename);
	case('-z')
		i=i+1;call get_command_argument(i,str); read(str,*) z
	end select
enddo

print*,trim(trim(path)//filelist)
print*,trim(trim(path)//outfilename)

open(inid,file=trim(path)//filelist)
open(outid,file=trim(path)//outfilename)
read(inid,*) n
print*,n
do i=1,n
	read(inid,*) filename
	call gr_on_cu_analysis(arr1,arr2,trim(path)//filename,z)
	write(*,'(i10)',advance='no') i
	write(outid,*) trim(filename),'    ',arr1-arr2(1)
enddo
close(inid)
close(outid)

end program run_gr_analysis
