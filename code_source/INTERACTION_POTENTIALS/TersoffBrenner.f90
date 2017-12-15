module TersoffBrenner
use cut_off_function
use md_general
implicit none

type TersoffBrenner_parameters
	real d,s,b,r0,delt,a0,c0,d0,R1,R2
	real c02,d02
end type TersoffBrenner_parameters

contains

subroutine read_TB_parameters(TBp,filename)
	type(TersoffBrenner_parameters):: TBp
	character(*):: 	filename
	
	open(1,file=filename)
	read(1,*) TBp%d,TBp%s,TBp%b,TBp%r0,TBp%delt,TBp%a0,TBp%c0,TBp%d0
	read(1,*) TBp%R1,TBp%R2
	close(1)
	TBp%c02 = TBp%c0**2
	TBp%d02 = TBp%d0**2
	
end subroutine read_TB_parameters

subroutine TB_energy(energy,nl,TBp)
	type(neibour_list):: nl
	type(TersoffBrenner_parameters):: TBp
	integer:: i,p,j,q
	real:: energy,energy_priv,B(nl%neib_num_max,nl%N),a

	energy = 0.
	energy_priv = 0.
	!$OMP PARALLEL firstprivate(energy_priv,i,p,j,q,a)
	!$OMP DO
	do i=1,nl%N
		B(:,i) = 0. !look rebosc 
		do p=1,nl%nnum(i)
			!B(p,i) = 0.
			if (nl%moddr(p,i)<TBp%R2) then
				do q=1,nl%nnum(i)
					if(p/=q .AND. nl%moddr(q,i)<TBp%R2) then
						B(p,i) = B(p,i)+f_cut(nl%moddr(q,i),TBp%R1,TBp%R2)*&
						(1.+TBp%c02/TBp%d02-TBp%c02/(TBp%d02+(1.+sum(nl%dr(:,p,i)*nl%dr(:,q,i))/(nl%moddr(p,i)*nl%moddr(q,i)))**2))	
					endif
				enddo
				B(p,i) = (1.+TBp%a0*B(p,i))**(-TBp%delt)
			endif
		enddo
	enddo
	!$OMP END DO	
	!$OMP DO
	do i=1,nl%N
		do p=1,nl%nnum(i)
			if (nl%moddr(p,i)<TBp%R2) then
				j = nl%nlist(p,i)
				if (j>i) then
					do q=1,nl%nnum(j)
						if(nl%nlist(q,j)==i) exit
					enddo
					a = -sqrt(2.*TBp%s)*TBp%b*(nl%moddr(p,i)-TBp%r0)
					energy_priv = energy_priv+f_cut(nl%moddr(p,i),TBp%R1,TBp%R2)*TBp%d/(TBp%s-1.)*(exp(a)-(B(p,i)+B(q,j))/2*TBp%s*exp(a/TBp%s))
				endif
			endif
		enddo
	enddo
	!$OMP END DO
	!$OMP ATOMIC
		energy = energy+energy_priv
	!$OMP END PARALLEL
	
end subroutine TB_energy

subroutine TB_forces(atoms,nl,TBp)
	type(particles)::	atoms
	type(neibour_list):: nl
	type(TersoffBrenner_parameters):: TBp
	integer:: i,p,j,q,l
	real:: B(nl%neib_num_max,nl%N),dB(3),a,dff,rr,cosi,f_c,df_c
	
	!$OMP PARALLEL firstprivate(i,p,j,q,l,a,dff,rr,cosi,f_c,df_c,dB)
	!$OMP DO
	do i=1,nl%N
		do p=1,nl%nnum(i)
			B(p,i) = 0.
			if (nl%moddr(p,i)<TBp%R2) then
				do q=1,nl%nnum(i)
					if(p/=q .AND. nl%moddr(q,i)<TBp%R2) then
						B(p,i) = B(p,i)+f_cut(nl%moddr(q,i),TBp%R1,TBp%R2)*&
						(1.+TBp%c02/TBp%d02-TBp%c02/(TBp%d02+(1.+sum(nl%dr(:,p,i)*nl%dr(:,q,i))/(nl%moddr(p,i)*nl%moddr(q,i)))**2))	
					endif
				enddo
				B(p,i) = (1.+TBp%a0*B(p,i))**(-TBp%delt)
			endif
		enddo
	enddo
	!$OMP END DO	
	!$OMP DO
	do i=1,nl%N	
		do p=1,nl%nnum(i)
			dB=0.d0
			do q=1,nl%nnum(i)
				if (p/=q) then
					rr = 1./nl%moddr(p,i)/nl%moddr(q,i)
					cosi = sum(nl%dr(:,p,i)*nl%dr(:,q,i))*rr
					dB = dB+ f_cut(nl%moddr(q,i),TBp%R1,TBp%R2)*2.d0*TBp%a0*TBp%c02*(1.+cosi)/(TBp%d02+(1.+cosi)**2)**2*&
						( (nl%dr(:,p,i)+nl%dr(:,q,i))*rr-cosi*(nl%dr(:,p,i)/nl%moddr(p,i)**2+nl%dr(:,q,i)/nl%moddr(q,i)**2) )&
						+ df_cut(nl%moddr(q,i),TBp%R1,TBp%R2)*TBp%a0*(1.+TBp%c02/TBp%d02-TBp%c02/(TBp%d02+(1.+cosi)**2))*nl%dr(:,q,i)
				endif
			enddo
			dB=dB*B(p,i)**(1./TBp%delt+1.)
			j=nl%nlist(p,i)
			do l=1,nl%nnum(j)
				if(nl%nlist(l,j)==i) exit
			enddo
			do q=1,nl%nnum(j)
				if (q/=l) then
					rr = 1./nl%moddr(l,j)/nl%moddr(q,j)
					cosi = sum(nl%dr(:,l,j)*nl%dr(:,q,j))*rr 
					dB = dB+(B(l,j)**(1./TBp%delt+1.)*f_cut(nl%moddr(q,j),TBp%R1,TBp%R2)*&
						2.*TBp%a0*TBp%c02*(1.+cosi)/(TBp%d02+(1.+cosi)**2)**2)*(-nl%dr(:,q,j)*rr + cosi*nl%dr(:,l,j)/nl%moddr(l,j)**2)
				endif
			enddo
			dB = dB*(-TBp%delt)/2.
			if (nl%moddr(p,i)<TBp%R2) then
				f_c = f_cut(nl%moddr(p,i),TBp%R1,TBp%R2)
				dff = df_cut(nl%moddr(p,i),TBp%R1,TBp%R2)/f_c
				a = -sqrt(2.*TBp%s)*TBp%b*(nl%moddr(p,i)-TBp%r0)
				atoms%forces(:,nl%particle_index(i)) = atoms%forces(:,nl%particle_index(i))+&
				f_c*TBp%d/(TBp%s-1.)*((nl%dr(:,p,i)*dff-sqrt(2.*TBp%s)*TBp%b/nl%moddr(p,i)*nl%dr(:,p,i))*exp(a)&
				-( db+(B(p,i)+B(l,j))/2*(nl%dr(:,p,i)*dff-sqrt(2./TBp%s)*TBp%b/nl%moddr(p,i)*nl%dr(:,p,i)) )*TBp%s*exp(a/TBp%s))
			endif
			do q=1,nl%nnum(j)
				if (q/=l) then	
					f_c = f_cut(nl%moddr(l,j),TBp%R1,TBp%R2)
					df_c = df_cut(nl%moddr(l,j),TBp%R1,TBp%R2)
					rr = 1./nl%moddr(l,j)/nl%moddr(q,j)
					cosi = sum(nl%dr(:,l,j)*nl%dr(:,q,j))*rr				
					atoms%forces(:,nl%particle_index(i)) = atoms%forces(:,nl%particle_index(i))+TBp%delt/2*B(q,j)**(1./TBp%delt+1.)*&
					( f_c*2.*TBp%a0*TBp%c02*(1.+cosi)/(TBp%d02+(1.+cosi)**2)**2 *( -nl%dr(:,q,j)*rr + cosi*nl%dr(:,l,j)/nl%moddr(l,j)**2 )+&
					nl%dr(:,p,i)*df_c*TBp%a0*(1.+TBp%c02/TBp%d02-TBp%c02/(TBp%d02+(1.+cosi)**2)) )*&
					f_cut(nl%moddr(q,j),TBp%R1,TBp%R2)*TBp%d/(TBp%s-1.)*TBp%s*exp(-sqrt(2.*TBp%s)*TBp%b*(nl%moddr(q,j)-TBp%r0)/TBp%s)
				endif
			enddo		
		enddo
	enddo
	!$OMP END DO
	!$OMP END PARALLEL
	
end subroutine TB_forces

end module TersoffBrenner