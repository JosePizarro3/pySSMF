!
!
! In this program, we are going to calculate, in a self-consistenly way, the quasiparticle
! weight Z_ms for a given filling and interactions U,U' and J_H, where we will consider 
! the rotational invariance:
!		U' = U - 2*JH
! in the single-site approximation of a U(1) slave-spin formalism.
!
!
! We will use the dispersion relation corresponding with the system in which we are; in this
! case, we will use the tight-binding model from Graser et al. NJP, 11 (2009). 
!
!
! Within the Slave-Spin formalism, we obtain 2 coupled hamiltonians: H_f (corresponding
! to the auxiliary fermions f and described by the spinons f_kmm') and H_ps (corresponding
! to the pseudospin and described by the matrices Sz_ms and xi_ms). We are going to self-
! -consistenly solve this problem by passing through H_f and H_PS 1 time each iteration step,
! Once we fulfill convergence criterias (for the constraint_ms, Z_ms and lambda_ms), we 
! will write Z_ms vs U in a file (as well as other interesting quantities).
!
!
! If we have N orbitals (in a paramagnetic case), the hamiltonians H_f and H_ps will have 
! the following dimensions:
!		H_f		-->		N x N
!		H_ps		-->		2^2N x 2^2N
! so that, it is convinient to perform an special construction (because of dim(H_PS)) when
! we have more than 2 orbitals in our problem is a hard issue (LAURA FANFARILLO Bit-method).
!
!
! Also, in the xi_ms matrices, we have to build in the following way:
!	xi_ms_conj = <P+_m> S+_ms <P-_m>
! where:
!	P+_m = 1.0/sqrt(0.5 + 1.0e-8 + Sz_ms)
!	P-_m = 1.0/sqrt(0.5 + 1.0e-8 - Sz_ms)
! So we will calculate the initial <P+_m><P-_m> using some results from Z2 slave-spin:
! 	<P+_m><P-_m> = c_gauge_ms + 1
!
!
! Finally, we have to vary the Lagrange multiplier lambda to fulfill the local constraint:
!		<nf_ms> = <Sz_ms> + 1/2		-->		constraint_ms = <nf_ms> - <Sz_ms> - 1/2
! in the following way:
!		lambda_ms = (1.0 - alpha*constraint_ms)*lambda_initial_ms
! Or in this other way:
!		lambda_ms = lambda_initial_ms - alpha*constraint_ms
! where alpha will be close to 0.1-1.0.
!
!
!
!
! ** Notice that all the averages are taken with respect to the hamiltonian where the magnitude
! is defined.
!
!





program slavespin_Norbital
	implicit none



! Some parameters and dimensions that we will use in the subroutines
	real*8,parameter :: pi = 4.0*atan(1.0)
	integer,parameter :: grid_u = 400
	integer,parameter :: grid_jh = 20
	integer,parameter :: grid = 500
	integer,parameter :: pointsk = grid + 1
	integer,parameter :: nk = pointsk**2
	integer,parameter :: sites = 1
	integer,parameter :: orbitals = 3
	integer,parameter :: spin = 2
	integer,parameter :: N_orb = orbitals*sites*spin
	integer,parameter :: N_bit = 2*N_orb/spin
	integer,parameter :: hd = 2**N_bit



! Interaction parameters.
	real*8 uint,uprime,jh				
	real*8 factorjh								! jh = factorjh*uint
	real*8 delta_uint,delta_jh							! deltas defined as variations of the corresponding quantity
	

! Slave-spin variables. 
	real*8 average_O_ms(N_orb),average_Sz_ms(N_orb)
	real*8 average_Sp_ms(N_orb),average_Pp_ms(N_orb),average_Pm_ms(N_orb)
	real*8 Z_ms(N_orb),lambda_ms(N_orb),correction_ms(N_orb),c_gauge_ms(N_orb),h_ms(N_orb)
	real*8 Z_ms_old(N_orb),lambda_ms_old(N_orb),correction_ms_old(N_orb)
	real*8 Z_ms_ps(N_orb),lambda_ms_ps(N_orb),correction_ms_ps(N_orb)
	real*8 var_lambda_ms(N_orb),var_chempot_ms(N_orb),tot_magn,magnetization_m(orbitals)

	real*8 lambda_ms_initial(N_orb)
	real*8 constraint_ms(N_orb),constraint_ms_initial(N_orb)

	real*8 filling_ms0(N_orb),lambda_ms0(N_orb),Z_ms0(N_orb),correction_ms0(N_orb)
	real*8 average_Pp_ms0(N_orb),average_Pm_ms0(N_orb),c_gauge_ms0(N_orb)
	real*8 constraint_ms0(N_orb)
	

! Control definitions in the convergence of our quantities.
	real*8 norm_Z(N_orb),norm_l(N_orb),norm_corr(N_orb)						! norm_X_ms = abs(X_ms - X_ms_old)
	real*8 norm_const(N_orb)									! norm_const_ms = abs(constraint_ms)
		

! Mixing quantities.
	real*8,parameter :: mix_Z = 0.1 		! X = (1.0 - mix)*X_old + mix*X_new
	real*8,parameter :: mix_l = 0.1
	real*8,parameter :: mix_corr = 0.1

	real*8 alpha(N_orb)							! lambda_ms = lambda_initial_ms - alpha*constraint_ms
	integer signoconst(N_orb),signoconst_initial(N_orb)
	
	
! Convergence criteria, iteration steps and controllers.			
	real*8 delta_auxf_Z,delta_auxf_l,delta_auxf_corr,delta_auxf
	real*8 delta_ps,delta_ps_init,delta_ps_min
	
	real*8 iter_auxf,iter_ps
	real*8,parameter :: iter_auxf_max = 1.0e4
	real*8,parameter :: iter_ps_max = 2.0e4
	integer control_auxf,control_ps,control_ps_min,control_exit
		
	
! Filling_ms and c_gauge_ms
	real*8 electrons,filling
	real*8 delta_fill
	real*8 filling_ms(N_orb)
	real*8 onsite_ms(N_orb)

	
! Tight-binding and fermion auxiliary eigenenergies, as well as the chemical potential 
! and the occupation number.
	complex*8 tight_binding(N_orb,N_orb,nk)
	real*8 eigenvalues_auxf(N_orb,nk)
	complex*8 weight_auxf(N_orb,N_orb,N_orb,nk)
	real*8 eigenvalues_auxf_tot(N_orb*nk),energies_auxf(N_orb*nk)
	real*8 enF,chem_pot
	real*8 sum_bandf(N_orb)
	complex*8 sum_orbf(N_orb,N_orb)		! sum_k {<->}, in the band (diagonal Hf) and in the orbital (non-diagonal Hf) basis


! We define paramters to use LAPACK diagonalization (soubrutine DSYEV) as follows.
	real*8 h_ps(hd,hd),h_ps_old(hd,hd),ham_sz_ms(hd,N_orb)
	integer l_dim,inf,gs_pos
	real*8 work(hd*(3+hd/2))
	real*8 eigenvector_ps(hd,hd)
	real*8 eigenvalues_ps(hd)
	real*8 groundstate_ps(hd),gsenergy_ps
	

! Correlation functions
	real*8 matS(hd,hd),matS2(hd,hd),matnT(hd,hd),matnT2(hd,hd)
	real*8 nintra(hd,hd),ninter(hd,hd),nintra2(hd,hd),ninter2(hd,hd)
	real*8 corr_spin,corr_totcharge,corr_intracharge,corr_intercharge,corr_totcharge_deg
	real*8 corr_spin0,corr_totcharge0,corr_intracharge0,corr_intercharge0,corr_totcharge_deg0
	real*8 av_S2,av_S,av_nT2,av_nT,av_nintra2,av_ninter2,av_nintra,av_ninter


! Declaration of integers	
	integer i_orb,j_orb,i_uint,i_jh,i_fill
	integer i,j,l,r,m,q,k,ltot
	
	integer lorden(N_orb*nk),l2,i_fermi


! Chars to name the files
	character(40) namei,nameorb,namee,namejh	
	
	
	
! Times involved (in seconds).	
	real start_global,finish_global
	real timemin_global,timemin_global_auxf,timemin_global_ps



!	We will consider onsite_ms = 0.0:
!	1-spinmayoritary -> 1
!	1-spinminoritary -> 2
	do i_orb=1,N_orb
		onsite_ms(i_orb) = 0.0
	end do


	
	
! We define the variatonal change as:
!		delta_X = X_max/grid
	delta_uint = 8.0*4.0/grid_u
	delta_jh = 1.0/(3.0*grid_jh)

	
	
! And the initial and maximun convergence criteria
	delta_auxf = 1.0e-3
	delta_auxf_Z = delta_auxf 
	delta_auxf_l = delta_auxf 
	delta_auxf_corr = delta_auxf 
	
	delta_ps_init = 1.0e-4
	delta_ps_min = 5.0e-6
	
	delta_ps = delta_ps_init



!\\\\\\\\\\\\\\//////////////!	
!\\\\\\\\\\\\\\//////////////!
!\\\\\ STARTING PROGRAM /////!
!\\\\\\\\\\\\\\//////////////!
!\\\\\\\\\\\\\\//////////////!




			
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! We define the total filling (per unit spin) as follows
	electrons = 2.00
	filling = electrons/N_orb		! filling = [0.0,1.0], where half-fill = 0.5 



! In order to calculate the effective fields h_ms(i,j) we will need a subroutine
! which calculates the dispersion relations without interactions. Then, we have 
! to obtain the tight-binding energies.
	call tightbinding(pi,grid,nk,N_orb,tight_binding)		
								
						
! WE USE THE DATA FROM TIGHT-BINDING RESULTS.	
! And also, we will calculate the reference filling when we are out of the slave-spin;
! This is because the technique usually creates some unphysical lambda_ms when U=0.0,
! so that we have to calculate in a iterative way the correction_ms and lambda_ms which
! are better compensated.
	call formula_correction(pi,grid,pointsk,nk,N_orb,filling,filling_ms0,lambda_ms0,onsite_ms)	


! And we will initially consider Z_ms = 1.0 and correction_ms = lambda_ms	
	do i_orb=1,N_orb
		Z_ms0(i_orb) = 1.0

		correction_ms0(i_orb) = lambda_ms0(i_orb)
		
		average_Pp_ms0(i_orb) = sqrt(1.0/sqrt(filling_ms0(i_orb)*(1.0 - filling_ms0(i_orb)) + 1.0e-8))
		average_Pm_ms0(i_orb) = sqrt(1.0/sqrt(filling_ms0(i_orb)*(1.0 - filling_ms0(i_orb)) + 1.0e-8))

		c_gauge_ms0(i_orb) = - 1.0 + average_Pm_ms0(i_orb)*average_Pp_ms0(i_orb)
		
		constraint_ms0(i_orb) = 1.0
				
! And we initialice the value of alpha_ms
		alpha(i_orb) = 0.5
	end do


			
	write(nameorb,1) N_orb
	write(namee,3) electrons			
! "Z..." saves the Z_ms values in terms of uint and "noconv..." saves the same, but when iter_ps and iter_global reach his maximun value.
	open(110, file="data_"//trim(adjustl(nameorb))//"orb_ne="//trim(adjustl(namee))//".dat")
	open(120, file="noconv_"//trim(adjustl(nameorb))//"orb_ne="//trim(adjustl(namee))//".dat")

	open(130, file="corr_data_"//trim(adjustl(nameorb))//"orb_ne="//trim(adjustl(namee))//".dat")
	open(140, file="corr_noconv_"//trim(adjustl(nameorb))//"orb_ne="//trim(adjustl(namee))//".dat")

	do i_orb=1,N_orb
		write(namei,2) i_orb
		


! "testL...", "testZ..." and "testcorr..." files save the convergence data during all the process	
		open(i_orb, file="testAUXF"//trim(adjustl(namei))//"_"//trim(adjustl(nameorb))//"orb_ne="//trim(adjustl(namee))//".dat")
!		open(i_orb+10, file="testPS"//trim(adjustl(namei))//"_"//trim(adjustl(nameorb))//"orb_ne="//trim(adjustl(namee))//".dat")
	end do



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!	<	<	<	<	<	<	<	<	<	SECTION 1	>	>	>	>	>	>	>	>   >	>	!
!	<	<	<	<	<	<	<	<	Initial declarations	>	>	>	>	>	>	>	>	!
	
	
	jh_loop: do i_jh=0,0
		factorjh = 0.25


		do i_orb=1,N_orb
			Z_ms(i_orb) = Z_ms0(i_orb)

			lambda_ms(i_orb) = lambda_ms0(i_orb)
			correction_ms(i_orb) = correction_ms0(i_orb)
		
			average_Pp_ms(i_orb) = average_Pp_ms0(i_orb)
			average_Pm_ms(i_orb) = average_Pm_ms0(i_orb)

			filling_ms(i_orb) = filling_ms0(i_orb)
			c_gauge_ms(i_orb) = c_gauge_ms0(i_orb)
		
			constraint_ms(i_orb) = constraint_ms0(i_orb)
		end do


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		control_exit = 1	
		uint_loop: do i_uint=0,grid_u			! we will take U=[0.0,5.0]

			call cpu_time(start_global)	
	
				
			uint = i_uint*delta_uint	
			jh = factorjh*uint
			uprime = uint - 2.0*jh

						

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!	<	<	<	<	<	<	<	<	<	SECTION 1	>	>	>	>	>	>	>	>   >	>	!
!	<	<	<	<	<	<	<	<	Diagonalization Hf	>	>	>	>	>	>	>	>   >	!


			control_auxf = 1
			control_ps_min = 1
			do i_orb=1,N_orb
				if(Z_ms(i_orb) <= 1.0e-2) then
					control_auxf = control_auxf + 3
					control_ps_min = control_ps_min + 1
				end if	
			end do
			if(control_auxf <= 3*N_orb) then
				control_auxf = 1
			end if
			if(control_ps_min <= N_orb) then
				control_ps_min = 1
			end if


			iter_auxf = 0.0
			iter_ps = 0.0
			loop_auxf: do while(iter_auxf < iter_auxf_max .and. (control_auxf <= 3*N_orb .or. control_ps_min <= N_orb))
											
				
				control_auxf = 1
				control_ps_min = 1
				iter_auxf = iter_auxf + 1.0	


! Now, we diagonalize the Hf hamiltonian
				call diago_auxf(pi,grid,nk,N_orb,Z_ms,onsite_ms,correction_ms,lambda_ms,eigenvalues_auxf,weight_auxf)
! With this subroutine, we have obtained eigenvalues_auxf(i_band,l) and eigenvectors_auxf(i_band,j_band,l)
! where l is the possible values of k(l<=nk).

! And we rename Z_ms_old and lambda_ms_old.
				do i_orb=1,N_orb
					average_O_ms(i_orb) = sqrt(Z_ms(i_orb))

					Z_ms_old(i_orb) = Z_ms(i_orb)
					lambda_ms_old(i_orb) = lambda_ms(i_orb)
					correction_ms_old(i_orb) = correction_ms(i_orb)
				end do				
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				
				
				
				
				
				
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!	<	<	<	<	<	<	<	<	<	SECTION 2	>	>	>	>	>	>	>	>   >	>	!
!	<	<	<	<	<	<	<	<	Fermi level in H_f	>	>	>	>	>	>	>	>   >	!

! Now, we have to save the eigenvalues in the same array "eigenvalues_auxf_tot(l)" in order
! to obtain the chem_pot.
				ltot=1
				do i=1,N_orb
					do l=1,nk
						eigenvalues_auxf_tot(ltot) = eigenvalues_auxf(i,l)
						ltot = ltot + 1
					end do
				end do
				
! Obtaining the Fermi level.		
				do l=1,N_orb*nk
					lorden(l) = l
				end do
				call indexx(N_orb*nk,eigenvalues_auxf_tot,lorden)
				 
				do l=1,N_orb*nk
					l2 = lorden(l)
					energies_auxf(l) = eigenvalues_auxf_tot(l2)
				end do
				i_fermi = int(filling*N_orb*nk)
	
				enF = energies_auxf(i_fermi)
				chem_pot = enF
					
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!	<	<	<	<	<	<	<	<	<	SECTION 3	>	>	>	>	>	>	>	>   >	>	>	!
!	<	<	<	<	<	<	<   Calculation h_eff and c_gauge	>	>	>	>	>	>	>	>	!


! Then, we obtain all the effective fields h_ms(i,j) using the following change of basis. 
! In the following subroutine is explained the method:					
				call calculation_eff_field(nk,N_orb,chem_pot,tight_binding,eigenvalues_auxf,weight_auxf,sum_orbf,average_O_ms,h_ms)
	
! And then, we can calculate the filling_ms = sum_orbf(i,i)/nk and c_gauge_ms (when U,JH = 0.0)
				do i_orb=1,N_orb
					filling_ms(i_orb) = real(sum_orbf(i_orb,i_orb))/nk
					
					c_gauge_ms(i_orb) = - 1.0 + 1.0/(sqrt(filling_ms(i_orb)*(1.0-filling_ms(i_orb)))+1.0e-8)
				end do
						
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! We first construct the Hps matrix....					
				call construct_h_ps(N_orb,N_bit,hd,h_ms,c_gauge_ms,h_ps_old,uint,uprime,jh,matS,matS2,matnT,matnT2,nintra,ninter,nintra2,ninter2)
				call construct_sz_ms(N_orb,N_bit,hd,ham_sz_ms)
			
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!	<	<	<	<	<	<	<	<	<	SECTION 4	>	>	>	>	>	>	>	>   >	>	!
!	<	<	<	<	<	<	<	<	Diagonalization Hps	>	>	>	>	>	>	>	>   >	!
				delta_ps = delta_ps_init/(iter_auxf**2)
				if(delta_ps <= delta_ps_min) then
					delta_ps = delta_ps_min
				end if	


				control_ps = 1
				do i_orb=1,N_orb
					if(Z_ms(i_orb) <= 1.0e-2) then
						control_ps = control_ps + 1
					end if	
				end do
				if(control_ps <= N_orb) then
					control_ps = 1
				end if

	
				iter_ps = 0.0
				loop_ps: do while(iter_ps < iter_ps_max .and. control_ps <= N_orb)! .and. control_exit<=N_orb/sites)
				
					control_ps = 1
					iter_ps = iter_ps + 1.0	

					do i_orb=1,N_orb
						lambda_ms_initial(i_orb) = lambda_ms(i_orb)

						constraint_ms_initial(i_orb) = constraint_ms(i_orb)
					end do


					do j=1,hd
						do k=1,hd	
							h_ps(j,k) = h_ps_old(j,k)
						end do
						do i_orb=1,N_orb
							h_ps(j,j) = h_ps(j,j) + ham_sz_ms(j,i_orb)*lambda_ms_initial(i_orb)
						end do
					end do
				

! ...then, we diagonalize the hamiltonian matrix...
					eigenvector_ps(:,:) = h_ps(:,:)
					l_dim = hd*(3+hd/2)
	
					call dsyev('V','U',hd,eigenvector_ps,hd,eigenvalues_ps,work,l_dim,inf)
! ... and finally, we asign to the ground state the eigenvector with the lowest eigenvalue.
					gs_pos = sum(MINLOC(eigenvalues_ps))
					do i=1,hd
						groundstate_ps(i) = eigenvector_ps(i,gs_pos)
					end do
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!<	<	<	<	<	<	<	<	<	SECTION 5	>	>	>	>	>	>	>	>   >	>	!
!<	<	<	<	<	<	<	Calculation of <O> and <Sz>	>	>	>	>	>	>	>	>	!
	
! Now, let's recalculate <O> and calculate <Sz> (over the groundstate_ps).
					call calculation_average_ps(hd,N_orb,N_bit,groundstate_ps,average_Sz_ms,average_Sp_ms)
					do i_orb=1,N_orb	
						average_Pp_ms(i_orb) = 1.0/sqrt(0.5 + 1.0e-8 + average_Sz_ms(i_orb))	
						average_Pm_ms(i_orb) = 1.0/sqrt(0.5 + 1.0e-8 - average_Sz_ms(i_orb))
				
						c_gauge_ms(i_orb) = - 1.0 + average_Pp_ms(i_orb)*average_Pm_ms(i_orb)

						average_O_ms(i_orb) = average_Pp_ms(i_orb)*average_Pm_ms(i_orb)*average_Sp_ms(i_orb)
						Z_ms_ps(i_orb) = average_O_ms(i_orb)**2

						correction_ms_ps(i_orb) = - 2.0*abs(average_O_ms(i_orb))*abs(h_ms(i_orb))*(filling_ms(i_orb) - 0.5)*((1.0 + c_gauge_ms(i_orb))**2)
					end do	
									
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





	 
	 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!<	<	<	<	<	<	<	<	<	SECTION 6	>	>	>	>	>	>	>	>   >	>  !
!<	<	<	<	<	<	<	<	Self-Cons. check	>	>	>	>	>	>	>	>   >  !	

! Now, we calculate constraint_ms
! We will obtain the relative values to convergence the problem as norm_X
					do i_orb=1,N_orb
						if(Z_ms_ps(i_orb) > 1.0e-2) then
							constraint_ms(i_orb) = average_Sz_ms(i_orb) + 0.5 - filling_ms(i_orb)
						else
							constraint_ms(i_orb) = 0.0
						
							Z_ms_ps(i_orb) = 0.0
							correction_ms_ps(i_orb) = 0.0
						end if
						norm_const(i_orb) = abs(constraint_ms(i_orb))

! And also, the sign of alpha is decided with respect to the change on
! the constraint_ms.
						if(abs(constraint_ms(i_orb)) > abs(constraint_ms_initial(i_orb))) then
							alpha(i_orb) = - alpha(i_orb)
						else
							alpha(i_orb) = alpha(i_orb)
						end if




! CONSTRAINT CONVERGENCE:			
						if(norm_const(i_orb) <= delta_ps .or. i_uint==0) then 
							control_ps = control_ps + 1
						end if
						if(abs(lambda_ms_initial(i_orb))>1.0e-1) then
							lambda_ms(i_orb) = (1.0 - alpha(i_orb)*norm_const(i_orb))*lambda_ms_initial(i_orb)
						else
							lambda_ms(i_orb) = lambda_ms_initial(i_orb) - alpha(i_orb)*norm_const(i_orb)
						end if					
						lambda_ms_ps(i_orb) = lambda_ms(i_orb)	


!						if(iter_ps>=iter_ps_max/10.0) then
!							write(i_orb+10,11) uint/8.0,iter_auxf,iter_ps,delta_ps,norm_const(i_orb),average_Sz_ms(i_orb),filling_ms(i_orb),lambda_ms_old(i_orb),lambda_ms_initial(i_orb),lambda_ms_ps(i_orb),Z_ms_ps(i_orb),correction_ms_ps(i_orb),c_gauge_ms(i_orb),h_ms(i_orb),average_O_ms(i_orb),average_Pp_ms(i_orb),average_Pm_ms(i_orb),average_Sp_ms(i_orb)
!						end if
					end do	
				end do loop_ps






! Z and LAMBDA CONVERGENCE:	
				do i_orb=1,N_orb
					if(norm_const(i_orb) <= delta_ps_min .or. i_uint==0) then 
						control_ps_min = control_ps_min + 1
					end if


					if(Z_ms(i_orb) > 1.0e-2) then
!						norm_Z(i_orb) = abs((Z_ms_ps(i_orb) - Z_ms_old(i_orb))/(Z_ms_old(i_orb) + 1.0e-8))
!						norm_l(i_orb) = abs((lambda_ms_ps(i_orb) - lambda_ms_old(i_orb))/(lambda_ms_old(i_orb) + 1.0e-8))
!						norm_corr(i_orb) = abs((correction_ms_ps(i_orb) - correction_ms_old(i_orb))/(correction_ms_old(i_orb) + 1.0e-8))
						norm_Z(i_orb) = abs(Z_ms_ps(i_orb) - Z_ms_old(i_orb))
						norm_l(i_orb) = abs(lambda_ms_ps(i_orb) - lambda_ms_old(i_orb))
						norm_corr(i_orb) = abs(correction_ms_ps(i_orb) - correction_ms_old(i_orb))
					else
						norm_Z(i_orb) = 0.0
						norm_l(i_orb) = 0.0
						norm_corr(i_orb) = 0.0
					end if

					Z_ms(i_orb) = (1.0 - mix_Z)*Z_ms_old(i_orb) + mix_Z*Z_ms_ps(i_orb)
					if(norm_Z(i_orb) <= delta_auxf_Z) then
						control_auxf = control_auxf + 1
					end if
							
					lambda_ms(i_orb) = (1.0 - mix_l)*lambda_ms_old(i_orb) + mix_l*lambda_ms_ps(i_orb)
					if(norm_l(i_orb) <= delta_auxf_l) then
						control_auxf = control_auxf + 1
					end if

					correction_ms(i_orb) = (1.0 - mix_corr)*correction_ms_old(i_orb) + mix_corr*correction_ms_ps(i_orb)
					if(norm_corr(i_orb) <= delta_auxf_corr) then
						control_auxf = control_auxf + 1
					end if

							
! And we save the convergence process data
					if(iter_auxf>=iter_auxf_max/10.0) then
						write(i_orb,12) uint/8.0,iter_auxf,iter_ps,norm_Z(i_orb),norm_l(i_orb),norm_corr(i_orb),Z_ms_old(i_orb),Z_ms_ps(i_orb),Z_ms(i_orb),lambda_ms_old(i_orb),lambda_ms_ps(i_orb),lambda_ms(i_orb),correction_ms_old(i_orb),correction_ms_ps(i_orb),correction_ms(i_orb),chem_pot,c_gauge_ms(i_orb),h_ms(i_orb),filling_ms(i_orb),average_Sz_ms(i_orb),average_O_ms(i_orb),av_S,tot_magn
					end if	

					

					if(Z_ms(i_orb) <= 1.0e-2) then
						Z_ms(i_orb) = 0.0
						correction_ms(i_orb) = 0.0
					end if	

					var_lambda_ms(i_orb) = -lambda_ms(i_orb) + correction_ms(i_orb)	
					var_chempot_ms(i_orb) = var_lambda_ms(i_orb) - chem_pot	
				end do

				tot_magn = 0.0
				do i_orb=1,orbitals
					magnetization_m(i_orb) = filling_ms(2*i_orb-1) - filling_ms(2*i_orb)	

					tot_magn = tot_magn + magnetization_m(i_orb)
				end do	
					



				av_S2 = 0.0
				av_S = 0.0
				av_nT2 = 0.0
				av_nT = 0.0
				av_nintra2 = 0.0
				av_nintra = 0.0
				av_ninter2 = 0.0
				av_ninter = 0.0
				do i=1,hd	
					av_S2 = av_S2 + matS2(i,i)*(groundstate_ps(i)**2)
					av_S = av_S + matS(i,i)*(groundstate_ps(i)**2)

					av_nT2 = av_nT2 + matnT2(i,i)*(groundstate_ps(i)**2)
					av_nT = av_nT + matnT(i,i)*(groundstate_ps(i)**2)

					av_nintra2 = av_nintra2 + nintra2(i,i)*(groundstate_ps(i)**2)
					av_nintra = av_nintra + nintra(i,i)*(groundstate_ps(i)**2)

					av_ninter2 = av_ninter2 + ninter2(i,i)*(groundstate_ps(i)**2)
					av_ninter = av_ninter + ninter(i,i)*(groundstate_ps(i)**2)
				end do

				corr_spin = av_S2 - av_S**2

				corr_totcharge = av_nT2 - av_nT**2

				corr_intracharge = av_nintra2 - av_nintra**2
				corr_intercharge = av_ninter2 - av_nintra*av_ninter

				corr_totcharge_deg = N_orb*(corr_intracharge + (N_orb - 1)*corr_intercharge)



!				if(av_S<1.0e-2) then
!					do i_orb=1,orbitals
!						Z_ms(2*i_orb-1) = (Z_ms(2*i_orb)+Z_ms(2*i_orb-1))/2
!						Z_ms(2*i_orb) = Z_ms(2*i_orb-1)
!
!						lambda_ms(2*i_orb-1) = (lambda_ms(2*i_orb)+lambda_ms(2*i_orb-1))/2
!						lambda_ms(2*i_orb) = lambda_ms(2*i_orb-1)
!
!						correction_ms(2*i_orb-1) = (correction_ms(2*i_orb)+correction_ms(2*i_orb-1))/2
!						correction_ms(2*i_orb) = correction_ms(2*i_orb-1)
!
!						filling_ms(2*i_orb-1) = (filling_ms(2*i_orb)+filling_ms(2*i_orb-1))/2
!						filling_ms(2*i_orb) = filling_ms(2*i_orb-1)
!
!						c_gauge_ms(2*i_orb-1) = (c_gauge_ms(2*i_orb)+c_gauge_ms(2*i_orb-1))/2
!						c_gauge_ms(2*i_orb) = c_gauge_ms(2*i_orb-1)
!
!						h_ms(2*i_orb-1) = (h_ms(2*i_orb)+h_ms(2*i_orb-1))/2
!						h_ms(2*i_orb) = h_ms(2*i_orb-1)
!
!						average_Pp_ms(2*i_orb-1) = (average_Pp_ms(2*i_orb)+average_Pp_ms(2*i_orb-1))/2
!						average_Pp_ms(2*i_orb) = average_Pp_ms(2*i_orb-1)
!
!						average_Pm_ms(2*i_orb-1) = (average_Pm_ms(2*i_orb)+average_Pm_ms(2*i_orb-1))/2
!						average_Pm_ms(2*i_orb) = average_Pm_ms(2*i_orb-1)
!					end do
!				end if	

!				if(i_uint==0) then
!					corr_spin0 = corr_spin
!					corr_totcharge0 = corr_totcharge
!					corr_intracharge0 = corr_intracharge
!					corr_intercharge0 = corr_intercharge
!					corr_totcharge_deg0 = corr_totcharge_deg 
!				end if

!				if(iter_auxf < iter_auxf_max) then
!					write(130,13) uint/8.0,corr_spin/corr_spin0,corr_totcharge/corr_totcharge0,corr_totcharge_deg/corr_totcharge_deg0,corr_intracharge,corr_intercharge,corr_spin0,corr_totcharge0,corr_totcharge_deg0,corr_intracharge0,corr_intercharge0,av_S2,av_S,av_nT2,av_nT,av_nintra2,av_ninter2,av_nintra,av_ninter
!				else
!					write(140,13) uint/8.0,corr_spin/corr_spin0,corr_totcharge/corr_totcharge0,corr_totcharge_deg/c!corr_totcharge_deg0,corr_intracharge,corr_intercharge,corr_spin0,corr_totcharge0,corr_totcharge_deg0,corr_intracharge0,corr_intercharge0,av_S2,av_S,av_nT2,av_nT,av_nintra2,av_ninter2,av_nintra,av_ninter
!				end if	
			end do loop_auxf
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




	

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!<	<	<	<	<	<	<	<	<	SECTION 7	>	>	>	>	>	>	>	>  	>	>  !
!<	<	<	<	<	<	<	<	Data save section & timers	>	>	>	>	>	>	>	>   	>  !



! We save the time that it tooks each self-consistenly problem and the number of iterations needed			
			call cpu_time(finish_global)
			timemin_global = (finish_global - start_global)/60.0	

			if(iter_auxf < iter_auxf_max) then
				write(130,13) uint/8.0,corr_spin,corr_totcharge,corr_totcharge_deg,corr_intracharge,corr_intercharge,av_S2,av_S,tot_magn,tot_magn/av_S,av_nT2,av_nT,av_nintra2,av_ninter2,av_nintra,av_ninter
			else
				write(140,13) uint/8.0,corr_spin,corr_totcharge,corr_totcharge_deg,corr_intracharge,corr_intercharge,av_S2,av_S,tot_magn,tot_magn/av_S,av_nT2,av_nT,av_nintra2,av_ninter2,av_nintra,av_ninter
			end if

! If all the values in the previously section 7 are OK, then we can keep the values of Z for that 
! given interaction ratio uint/bandwidth.	
			if(iter_auxf < iter_auxf_max) then
				write(110,14) uint/8.0,factorjh,Z_ms,lambda_ms,correction_ms,var_lambda_ms,var_chempot_ms,2.0*filling_ms,magnetization_m,abs(magnetization_m),1.0+c_gauge_ms,h_ms,average_Pp_ms,average_Pm_ms,chem_pot,iter_auxf,iter_ps,timemin_global
			else
				write(120,15) uint/8.0,factorjh,Z_ms,lambda_ms,correction_ms,var_lambda_ms,var_chempot_ms,2.0*filling_ms,magnetization_m,abs(magnetization_m),1.0+c_gauge_ms,h_ms,average_Pp_ms,average_Pm_ms,chem_pot,iter_auxf,iter_ps,timemin_global,constraint_ms,norm_Z,norm_l,norm_corr
			end if
		end do uint_loop
	end do jh_loop

			
		
		
	1 format(i2)
	2 format(i1)
	3 format(f4.2)


	11 format(f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8)
	12 format(f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8)
	13 format(f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8)
	14 format(f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8)
	15 format(f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8,f20.8)

	
	do i_orb=1,N_orb
		close(i_orb)
!		close(i_orb+10)
	end do
!	close(100)
	close(110)
	close(120)
	close(130)
	close(140)
	

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end program slavespin_Norbital













!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!										FUNCTIONS								!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! We define now the function "bit-spin" which changes a logical variable in a PS state.  
function bit_spin(b)
	implicit none 

	logical b
	integer bit_spin

! We want to describe:
!	b true --> bit_spin = +1 (PS up)
!	b false --> bit_spin = -1 (PS down)

	if(b) then
		bit_spin =  1
	else 
		bit_spin = - 1 
	end if
end function bit_spin





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now, we need an special logical step in order to differentiate from diagonal terms (which
! will build the Sz_ms terms) and the off-diagonal terms (which will build the O_ms terms).
! We have to use an XOR (with function IEOR and a compilation term "-fpscomp logicals")
! between the bits "n" and "m"; this XOR compares the bits of n,m and counts the times that
! the bits flip.
! TRANSFER is the same as CASTING one type of variable to another in C-language -> it makes the
! sum function to value between 0 and imax.
! When fliping=0 we will have diagonal terms, and when fliping=1 we will have off-diagonal
! terms (check 'construc_h_ps' subroutine).
function fliping(n,m,imax)
	implicit none 
	
	integer n,m,imax,i,sum
	integer fliping
	
	sum = 0
	do i=0,imax
	  	  sum = sum + TRANSFER(BTEST(IEOR(n,m),i),0)
	end do
	fliping = sum
end function fliping















	




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!										SUBROUTINES								!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! We will write down a subroutine to define the hoppings and crystal-field parameters.
! We will define:
!	A -> 134	Fe1
!	B -> 134	Fe1
!	C -> 25		Fe1
!	D -> 25		Fe1
!	E -> 134	Fe1
!	F -> 134	Fe2
!	G -> 134	Fe2
!	H -> 25		Fe2
!	I -> 25		Fe2
!	J -> 134	Fe2
subroutine construct_ham(N_orb,kx,ky,ham)
	implicit none
	
	
	integer N_orb
	real*8 kx,ky
	complex*16 ham(N_orb,N_orb)
	
	
	integer i,j,l,m,n
	
	complex*16,parameter :: zi = (0.d0,1.d0)
	
! Initialize hamiltonian to zero (complex)
	do i=1,N_orb
		do j=1,N_orb
			ham(i,j)=(0.d0,0.d0)
		end do
	end do	


	do m=1,N_orb
		ham(m,m) = -2*(cos(kx) + cos(ky))	
	end do

end subroutine construct_ham








!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! We will calculate the tight-binding solution
subroutine tightbinding(pi,grid,nk,N_orb,tight_binding)
	implicit none
	
	
	real*8 pi
	integer grid
	integer nk
	integer N_orb
	complex*8 tight_binding(N_orb,N_orb,nk)
	
	
	complex*16 ham(N_orb,N_orb)
	real*8 kx,ky
	real*8 deltak
	
	integer i,j,l,i_orb,j_orb
	
	deltak = 2.0*pi/real(grid+1)
	
	
	
! We map our BZ onto a grid in X and Y directions as follows; we will
! consider kx,ky = [-pi,pi], and then, we will keep the values of the 
! hamiltonian (see construct_ham) in the array tight_binding(l), where
! l will give us the k-dependence.
	l = 1
	do i=0,grid
		kx = - pi + i*deltak
		do j=0,grid
			ky = - pi + j*deltak
		
			call construct_ham(N_orb,kx,ky,ham)
			
			do i_orb=1,N_orb
				do j_orb=1,N_orb
					tight_binding(i_orb,j_orb,l) = ham(i_orb,j_orb)
				end do
			end do
			l = l + 1
		end do
	end do
	
!	write(*,*) "Tight-binding calculated!!"
end subroutine tightbinding









!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! We will calculate the tight-binding solution
subroutine formula_correction(pi,grid,pointsk,nk,N_orb,filling,filling_ms,lambda_ms,onsite_ms)
	implicit none



! Some parameters and dimensions that we will use in the subroutines
	real*8 pi 
	integer grid
	integer pointsk 
	integer nk 
	integer N_orb
	real*8 filling
	real*8 filling_ms(N_orb)
	real*8 lambda_ms(N_orb)
	real*16 onsite_ms(N_orb)
	
	
	
! Tight-binding and fermion auxiliary eigenenergies, as well as the chemical potential 
! and the occupation number	
	real*8 kx,ky
	real*8 deltak
	complex*16 tightbinding(N_orb,N_orb,nk)
	real*8 eigenvalues_tb(N_orb,nk)
	complex*8 weight_tb(N_orb,N_orb,N_orb,nk)
	
	real*8 eigenvalues_tb_tot(N_orb*nk)
	real*8 energies_tb(N_orb*nk)
	real*8 enF,chem_pot
	
	complex*8 sum_orbf(N_orb,N_orb), sum_f(N_orb,N_orb)					! sum_k {<->}, in the band (diagonal Hf) and in the orbital (non-diagonal Hf) basis
	complex*8 nkorb_f(N_orb,N_orb,nk)
	real*8 nkband_f(N_orb,nk)
	real*8 eff_field(N_orb)
	real*8 c_gauge_ms(N_orb)								! We will define it as beta = 1.0/kB*T    where kB = 8.6173324e-5 eV/K
	real*8 average_O_ms(N_orb)

	real*8 h_ms(N_orb),c_ms(N_orb),a_ms(N_orb),n_ms(N_orb)


! Declaration of integers	
	integer i_orb,j_orb,i_band
	integer i,j,l,r,m,q,k,ltot
	
	integer lorden(N_orb*nk),l2,i_fermi
	
	real start_diago,finish_diago,start_chempot,finish_chempot
	real timemin_diago,timemin_chempot

! We define paramters to use LAPACK diagonalization (soubrutine zheev) as follows
	complex*16 eigenvectors(N_orb,N_orb), h(N_orb,N_orb)
	integer lwork
	real*8 eigenenergies(N_orb)
	real*8 rwork(3*N_orb - 2)
	complex*16 work(2*N_orb - 1)
	character jobz,uplo
	integer lda,info	
	jobz = 'V'
	uplo = 'U'
	lda = N_orb
	lwork = 2*N_orb - 1	
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now we diagonalize the tight-binding hamiltonian, obtaining the eigenvalues and eigenvectors.		
	deltak = 2.0*pi/real(grid+1)
	

	call cpu_time(start_diago)

	l = 1
	do i=0,grid
		kx = - pi + i*deltak

		do j=0,grid
			ky = - pi + j*deltak

	

! We 1st construct the tight-binding hamiltonian for this lattice.	
			call construct_ham(N_orb,kx,ky,h)
			
! Then, we write in the aux_fermion basis.
			do i_orb=1,N_orb
				do j_orb=1,N_orb
					tightbinding(i_orb,j_orb,l) = h(i_orb,j_orb)
				end do
				h(i_orb,i_orb) = h(i_orb,i_orb) + onsite_ms(i_orb)
				do j_orb=1,N_orb
					eigenvectors(i_orb,j_orb) = h(i_orb,j_orb)
				end do
			end do
			
			
	
! Finally, we diagonalize using LAPACK subroutine zheev.	
			call zheev(jobz,uplo,lda,eigenvectors,lda,eigenenergies,work,lwork,rwork,info)
			
! And we save the eigenvalues and eigenvectors in their corresponding array.				
			do i_band=1,N_orb
				eigenvalues_tb(i_band,l) = eigenenergies(i_band)
				
				do i_orb=1,N_orb
					do j_orb=1,N_orb
						weight_tb(i_orb,j_orb,i_band,l) = conjg(eigenvectors(i_orb,i_band))*eigenvectors(j_orb,i_band)
					end do
				end do
			end do
			l = l + 1
		end do
	end do
	




	call cpu_time(finish_diago)
	timemin_diago = (finish_diago - start_diago)/60.0
	print '("I have spent = ",f10.2," min (DIAGONALIZING H0)")',timemin_diago


	call cpu_time(start_chempot)
	do i_orb=1,N_orb
		average_O_ms(i_orb) = 1.0
	end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
				
! Now, we have to save the eigenvalues in the same array "eigenvalues_auxf_tot(l)" in order
! to obtain the chem_pot.
	ltot=1
	do i=1,N_orb
		do l=1,nk
			eigenvalues_tb_tot(ltot) = eigenvalues_tb(i,l)
			
			ltot = ltot + 1
		end do
	end do
				
		
	do l=1,N_orb*nk
		lorden(l) = l
	end do
	call indexx(N_orb*nk,eigenvalues_tb_tot,lorden)

	do l=1,N_orb*nk
		l2 = lorden(l)
		energies_tb(l) = eigenvalues_tb_tot(l2)
	end do
	i_fermi = int(filling*N_orb*nk)
	
	enF = energies_tb(i_fermi)
	chem_pot = enF


	call cpu_time(finish_chempot)
	timemin_chempot = (finish_chempot - start_chempot)/60.0
	print '("I have spent = ",f10.2," min (CALCULATING CHEM_POT)")',timemin_chempot
	write(*,*) "Chemical potential =",chem_pot


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



! 1st, we will fill the bands as follows:
	do i_band=1,N_orb
		do l=1,nk
			if(eigenvalues_tb(i_band,l) < chem_pot) then
				nkband_f(i_band,l) = 1.0
			else if(eigenvalues_tb(i_band,l) > chem_pot) then
				nkband_f(i_band,l) = 0.0
			else
				nkband_f(i_band,l) = 0.5
			end if
		end do
	end do


	
! And then, we obtain the sum_orbf_mns following the next expressions (see personal notes):
	do i_orb=1,N_orb
		eff_field(i_orb) = 0.0
		do j_orb=1,N_orb
			sum_orbf(i_orb,j_orb) = 0.0
			do l=1,nk
				nkorb_f(i_orb,j_orb,l) = 0.0
				do i_band=1,N_orb
					nkorb_f(i_orb,j_orb,l) = nkorb_f(i_orb,j_orb,l) + weight_tb(i_orb,j_orb,i_band,l)*nkband_f(i_band,l)
				end do
				sum_orbf(i_orb,j_orb) = sum_orbf(i_orb,j_orb) + nkorb_f(i_orb,j_orb,l)
				sum_f(i_orb,j_orb) = sum_f(i_orb,j_orb) + tightbinding(i_orb,j_orb,l)*nkorb_f(i_orb,j_orb,l)
			end do
			eff_field(i_orb) = eff_field(i_orb) + abs(average_O_ms(j_orb))*real(sum_f(i_orb,j_orb))/nk
		end do
		filling_ms(i_orb) = real(sum_orbf(i_orb,i_orb))/nk
		c_ms(i_orb) = 1.0/(filling_ms(i_orb)*(1.0 - filling_ms(i_orb)) + 1.0e-8)

		lambda_ms(i_orb) = - 2.0*abs(eff_field(i_orb))*(filling_ms(i_orb) - 0.5)*c_ms(i_orb)
!		write(*,3000) i_orb,onsite_ms(i_orb),2.0*filling_ms(i_orb),lambda_ms(i_orb),c_ms(i_orb),eff_field(i_orb)
	end do
	

!	3000 format(i3,f20.8,f20.8,f20.8,f20.8,f20.8)
end subroutine formula_correction

















!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now we calculate the energies of the auxiliary fermions as follows.
subroutine diago_auxf(pi,grid,nk,N_orb,Z,onsite,correction,lambda,eigenvalues_auxf,weight_auxf)
	implicit none
	
	
	real*8 pi
	integer grid
	integer nk
	integer N_orb
	real*8 Z(N_orb)
	real*16 onsite(N_orb)
	real*8 correction(N_orb)
	real*8 lambda(N_orb)
	real*8 eigenvalues_auxf(N_orb,nk)
	complex*8 weight_auxf(N_orb,N_orb,N_orb,nk)
	
	real*8 kx,ky
	real*8 deltak
	
	integer i,j,l,i_orb,j_orb,i_band


! We define paramters to use LAPACK diagonalization (soubrutine zheev) as follows
	complex*16 eigenvectors(N_orb,N_orb), h(N_orb,N_orb)
	integer lwork
	real*8 eigenenergies(N_orb)
	real*8 rwork(3*N_orb - 2)
	complex*16 work(2*N_orb - 1)
	character jobz,uplo
	integer lda,info	
	jobz = 'V'
	uplo = 'U'
	lda = N_orb
	lwork = 2*N_orb - 1
	
	
	
	
	deltak = 2.0*pi/real(grid+1)
	


	l = 1
	do i=0,grid
		kx = - pi + i*deltak

		do j=0,grid
			ky = - pi + j*deltak


			call construct_ham(N_orb,kx,ky,h)
			
			do i_orb=1,N_orb
				do j_orb=1,N_orb
					h(i_orb,j_orb) = sqrt(Z(i_orb)*Z(j_orb))*h(i_orb,j_orb)
				end do
! Now, we add the onsite energies and lambda
				h(i_orb,i_orb) = h(i_orb,i_orb) + onsite(i_orb) - lambda(i_orb) + correction(i_orb)
			end do	


			do i_orb=1,N_orb
				do j_orb=1,N_orb
					eigenvectors(i_orb,j_orb) = h(i_orb,j_orb)
				end do
			end do
			
			call zheev(jobz,uplo,lda,eigenvectors,lda,eigenenergies,work,lwork,rwork,info)
						
			do i_band=1,N_orb
				eigenvalues_auxf(i_band,l) = eigenenergies(i_band)
					
				do i_orb=1,N_orb
					do j_orb=1,N_orb
						weight_auxf(i_orb,j_orb,i_band,l) = conjg(eigenvectors(i_orb,i_band))*eigenvectors(j_orb,i_band)
					end do
				end do
			end do
			l = l + 1
		end do
	end do



	
!	write(*,*) "Hf diagonalization on!!"
end subroutine diago_auxf

















!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now, we define the subroutine to order the eigenvalues in order to obtain the chemical potential.
subroutine indexx(n,arr,indx)
	implicit real*8 (a-h,o-z)
    dimension  arr(n),indx(n)
	PARAMETER (M=7,NSTACK=50)
	dimension istack(NSTACK)
    do 11 j=1,n
		indx(j)=j
11	continue
    jstack=0
    l=1
    ir=n
1	if(ir-l.lt.M) then
	do 13 j=l+1,ir
		indxt=indx(j)
		a=arr(indxt)
        do 12 i=j-1,1,-1
			if(arr(indx(i)).le.a) goto 2
			indx(i+1)=indx(i)
12			continue
			i=0
2			indx(i+1)=indxt
13		continue
		if(jstack.eq.0)return
		ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
end subroutine indexx




















!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! We will calculate sum_orbf_mns = sum_k (<f_kms†  f_kns>) in order to obtain the eff fields.
subroutine calculation_eff_field(nk,N_orb,chem_pot,tight_binding,eigenvalues_auxf,weight_auxf,sum_orbf,average_O,eff_field)
	implicit none


	integer nk
	integer N_orb,sites
	real*8 chem_pot
	complex*8 tight_binding(N_orb,N_orb,nk)
	real*8 eigenvalues_auxf(N_orb,nk)
	complex*8 weight_auxf(N_orb,N_orb,N_orb,nk)
	complex*8 sum_orbf(N_orb,N_orb)
	real*8 average_O(N_orb)
	real*8 eff_field(N_orb)
	

	complex*8 nkorb_f(N_orb,N_orb,nk)
	real*8 nkband_f(N_orb,nk)					! nk... = <...>
	complex*8 sum_f(N_orb,N_orb)								! sum_f = sum_k [tight_binding(k) * <f_kms† f_kns>]


	integer i_orb,j_orb,i_band,l




! 1st, we will fill the bands as follows:
	do i_band=1,N_orb
		do l=1,nk
			if(eigenvalues_auxf(i_band,l) < chem_pot) then
				nkband_f(i_band,l) = 1.0
			else if(eigenvalues_auxf(i_band,l) > chem_pot) then
				nkband_f(i_band,l) = 0.0
			else
				nkband_f(i_band,l) = 0.5
			end if
		end do
	end do

	
! And then, we obtain the sum_orbf_mns following the next expressions (see personal notes):
	do i_orb=1,N_orb
		eff_field(i_orb) = 0.0
		do j_orb=1,N_orb
			sum_orbf(i_orb,j_orb) = 0.0

			sum_f(i_orb,j_orb) = 0.0
			do l=1,nk
				nkorb_f(i_orb,j_orb,l) = 0.0
				do i_band=1,N_orb
					nkorb_f(i_orb,j_orb,l) = nkorb_f(i_orb,j_orb,l) +  weight_auxf(i_orb,j_orb,i_band,l)*nkband_f(i_band,l) !+ conjg(eigenvectors_auxf(i_orb,i_band,l))*eigenvectors_auxf(j_orb,i_band,l)*nkband_f(i_band,l)
				end do
				sum_orbf(i_orb,j_orb) = sum_orbf(i_orb,j_orb) + nkorb_f(i_orb,j_orb,l)
				
				sum_f(i_orb,j_orb) = sum_f(i_orb,j_orb) + tight_binding(i_orb,j_orb,l)*nkorb_f(i_orb,j_orb,l)
			end do
!			eff_field(i_orb) = eff_field(i_orb) + abs(average_O(j_orb))*real(sum_f(i_orb,j_orb))/nk
			eff_field(i_orb) = eff_field(i_orb) + average_O(j_orb)*real(sum_f(i_orb,j_orb))/nk
		end do
	end do		

end subroutine calculation_eff_field

















!LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
!	METHOD: The BIT MAPPING (***extracted from Laura Fanfarillo notes***)
!	The Hilbert space of 2N_orb pseudospin admits 2**{2N_orbitals} states 
!	corresponding to the possible arrangements of all the pseudospins. 
!	Thus every states can be written as a vector of 2N_orbitals components: 
!	|Sz_{m,sigma}> = |Sz_{1-},Sz_{1+},...Sz_{N_orb,-},Sz_{N_orb,+}> 
!	Using the mapping: 
!	Sz up <-> TRUE		Sz down <-> FALSE
!	We can associate every |Sz_{m,sigma}> vector to a 2N_orb-bit string i.e.
!	every state to an integer number I 
!	|I>=|bit_{i}> 
!	with the convention:  
!	I	|	m      sigma
!	0	|	1	-
!	1	|	1	+
!	2	|	2	-
!	3	|	2	+
!	...	|	...	...
!	2N_orb-1|	N_orb	+
!
!	Example1: 
!	Given the state |I>, the value of Sz_{M,S} is given by the value of the 
!	bit in position pos((M-1)*2) if S= -, pos((M-1)*2 + 1) if S= + .
!	Example2: 
!	Given the state |I>, the value of the bit in position pos(i) gives the 
!	value of Sz of the pseudospin in the orbital m=(i/2)+1 and with spin 
!	sigma = -/+ if i is even/odd respectively. 
!
!	Following this procedure the states are already ordered simply counting 
!	from |0> up to |2N_orbital-1> (-1 since we start from 0). Moreover the 
!	integers can be use as labels for the Hamiltonian elements always having 
!	in mind that the first line(column) j=1 is given by the state <0|(|0>)  
!	i.e. j=I+1. 
!	NB: This mapping allows us to treat a problem of N-spin with N <= 31 
!	since a standard integer has 31 bits available (the 32-th is useless 
!	since is the sign of I). In our multiorbila case where N= 2N_orbitals 
!	this means N_orbitals <=15. 
!	NB: In what follows we use these intrinsic bit manipulation functions: 
!	BTEST(I, POS)	.TRUE. if the position number POS of I is 1
!	IEOR(I, J)		gives an INTEGER whose bits are the logical exclusive OR 
!					between I and J (i.e. the bits series of this new integer
!					has .FALSE. in each position but the ones in which I and J 
!					has a different bit value)
!	IBCHNG(I,POS)	reverses the value of the bit POS in the integer I, 
!					gives back a new INTEGER (this function needs the intel 
!					compiler)
!LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
!FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! We now define the H_ps matrix with the following subroutine (taken from L.Fanfarillo notes).
subroutine construct_h_ps(N_orb,N_bit,hd,eff_field,c_gauge,h_ps,uint,uprime,jh,matS,matS2,matnT,matnT2,nintra,ninter,nintra2,ninter2)
	implicit none
	
	integer N_orb
	integer N_bit
	integer hd
	real*8 eff_field(N_orb)
	real*8 c_gauge(N_orb)
	real*8 h_ps(hd,hd)
	real*8 matS(hd,hd),matnT(hd,hd)
	real*8 matS2(hd,hd),matnT2(hd,hd)
	real*8 ninter(hd,hd),nintra(hd,hd)
	real*8 ninter2(hd,hd),nintra2(hd,hd)
	real*8 uint,uprime,jh
	
! We will consider the following arrays in order to use them with the bit-logic.
	real*8 eff_field_ms(N_bit),c_gauge_ms(N_bit)


	integer i,j,k,l,f,i_orb,i_bit
	integer m,s
! i_max is the maximun number of bits that we need starting from the bit 0-th
	integer i_max	

! Definition of pseudospin Arrays
	integer Sz_ms(N_bit)
	integer Sz_m(N_orb)
	integer Sz_s(2)

	real*8 sum_tot,sumu1,sumu2,sumjh
	real*8 ma,na,na1,na2
	real*8 sum_s,sum_n,sum_m
	integer flag
	
	integer bit_spin,fliping


	complex*8 sumtrans
	

! Now we write the values of lambda_ms, h_ms and c_ms for every orbital and spin
! degree of freedom as follows
	do i=1,N_orb
		c_gauge_ms(i) = c_gauge(i)
		

		eff_field_ms(i) = eff_field(i)
	end do
	


	i_max = N_bit - 1
	do j=1,hd
		do k=1,hd
			sum_tot = 0.0

! We will obtain the sums described in the hamiltonian Hps (see notes):
!	sumszdiag	---->	\sum_ms [Sz_ms]
!	sumu1 		---->	0.5*uprime*[\sum_ms Sz_ms]**2
!	sumu2		---->	0.5*(uint-uprime)* \sum_m [\sum_s Sz_ms]**2
!	sumjh 		---->	0.5*jh* \sum_s [\sum_m Sz_ms]**2
			sumu2 = 0.0
			sumjh = 0.0 
			sumu1 = 0.0
			
			sum_n = 0.0
			sum_m = 0.0
			sum_s = 0.0

			ma = 0.0
			na = 0.0

			na1 = 0.0
			na2 = 0.0

! "flips" account for the diagonal elements (if flips=0), off-diagonal
! elements (if flips=1) and integer>1 elsewhere.
			f = fliping(j-1,k-1,i_max)
			
			
!!!!!!!!!!!!!!!!!!!!!!!	DIAGONAL
			if (f==0) then 
! We first write the Sz_ms matrices in the basis described in the introduction (note
! that is defined without the "1/2" term of the matrix).
				do i=0,i_max
					Sz_ms(i+1) = bit_spin(BTEST(j-1,i))
				end do
				
				na1 = na1 + 0.5*(Sz_ms(1) + Sz_ms(2)) + 1.0
				na2 = na2 + 0.5*(Sz_ms(3) + Sz_ms(4)) + 1.0
				do l=1,N_bit/2
					ma = ma + 0.5*(Sz_ms(2*l) - Sz_ms(2*l-1)) 	
					na = na + 0.5*(Sz_ms(2*l) + Sz_ms(2*l-1)) + 1.0
				end do 

				do l=1,N_bit
					sumu1 = sumu1 + 0.5*Sz_ms(l) 
				end do 
				sumu1 = sumu1**2

! Now, in order to make it easier, we are going to define Sz_m and Sz_s as:
!	Sz_m	---->	\sum_s Sz_ms
!	Sz_s	---->	\sum_m Sz_ms
				do m=1,N_orb/2
					Sz_m(m) = Sz_ms(2*m-1) + Sz_ms(2*m)
				
					sum_n = sum_n + Sz_ms(2*m-1)
					sum_m = sum_m + Sz_ms(2*m)
				end do
				Sz_s(1) = sum_n
				Sz_s(2) = sum_m

				do m=1,N_orb/2 
					sumu2 = sumu2 + (0.5*Sz_m(m))**2
				end do
				do s=1,2 
					sumjh = sumjh + (0.5*Sz_s(s))**2
				end do

! Then, we finally obtain the sum_tot
				sum_tot = 0.5*uint*sumu2 + 0.5*uprime*(sumu1 - sumu2) - 0.5*jh*sumjh	
	
!!!!!!!!!!!!!!!!!!!!!!!	OFF-DIAGONAL
			else if(f==1) then 
! This loop (with IEOR) locate the position of the only one TRUE bit of the integer.
				do i=0,i_max
					if(BTEST(IEOR(j-1,k-1),i)) then
						flag = i
					end if 
				end do
				

				sum_tot = (1.0 + c_gauge_ms(flag+1))*eff_field_ms(flag+1)

				ma = 0.0
				na = 0.0

				na1 = 0.0
				na2 = 0.0
			else
				sum_tot = 0.0

				ma = 0.0
				na = 0.0

				na1 = 0.0
				na2 = 0.0
			end if
! Then, every term in the hamiltonian Hps will be this sum_tot that we calculated.
			h_ps(j,k) = sum_tot

			matS(j,k) = ma
			matS2(j,k) = ma**2
			matnT(j,k) = na
			matnT2(j,k) = na**2

			nintra(j,k) = na1
			nintra2(j,k) = na1**2
			ninter(j,k) = na2
			ninter2(j,k) = na1*na2
		end do 
	end do

end subroutine construct_h_ps














!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! We now define the H_ps matrix with the following subroutine (taken from L.Fanfarillo notes).
subroutine construct_sz_ms(N_orb,N_bit,hd,ham_sz_ms)
	implicit none
	
	integer N_orb
	integer N_bit
	integer hd
	real*8 ham_sz_ms(hd,N_orb)
	

	integer i,j,k,l,f,i_orb,i_bit
	integer m,s
! i_max is the maximun number of bits that we need starting from the bit 0-th
	integer i_max	

! Definition of pseudospin Arrays
	integer Sz_ms(N_bit)

	real*8 suml(N_orb)

	integer bit_spin



	i_max = N_bit - 1
	do j=1,hd
! We will obtain the sums described in the hamiltonian Hps (see notes):
!	suml 	---->	\sum_ms [Sz_ms + 1/2]
		do i=1,N_orb
			suml(i) = 0.0
		end do
			
!!!!!!!!!!!!!!!!!!!!!!!	DIAGONAL
! We first write the Sz_ms matrices in the basis described in the introduction (note
! that is defined without the "1/2" term of the matrix).
		do i=0,i_max
			Sz_ms(i+1) = bit_spin(BTEST(j-1,i))
		end do
		
		do i=1,N_orb
			suml(i) = suml(i) + 0.5*(Sz_ms(i) + 1.0)
		end do
		

		do i=1,N_orb
			ham_sz_ms(j,i) = suml(i)
		end do
	end do

end subroutine construct_sz_ms


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! We will calculate the averages <O> and <Sz> over the Hps in order to obtain the parameters involved 
! in the self-consistence process (check notes to see how to calculate <->ps).
subroutine calculation_average_ps(hd,N_orb,N_bit,groundstate_ps,average_Sz,average_Sp)
	implicit none

	integer hd
	integer N_orb
	integer N_bit
	real*8 groundstate_ps(hd)
	real*8 average_Sz(N_orb)
	real*8 average_Sp(N_orb)

	integer bit_spin,b
	integer i,j,k,l
	integer i_max

 	real*8 sum_Sz,sum_Sp
	real*8 average_Sz_ms(N_bit),average_Sp_ms(N_bit)




! We will calculate the sums (sum_Sz and sum_Sp) over <k-1| and |j-1>
	i_max = N_bit - 1
	do i=0,i_max
		average_Sz_ms(i+1) = 0.0
		average_Sp_ms(i+1) = 0.0

		sum_Sz = 0.0
		sum_Sp = 0.0
! Our ground state can be written as |GS> = sum_j groundstate_ps(j) |j>; Thus we have to compute: 
!		<GS|Sz_ms|GS> = sum_jk groundstate_ps(j)*groundstate_ps(k) <k-1|Sz_ms|j-1> = sum_j groundstate_ps(j)**2 Sz_ms|_(j-1)
! where Sz_ms|_(j-1) = +/- and:
!		<GS|S+-_ms|GS> = sum_jk groundstate_ps(k)*groundstate_ps(j) <k-1|S+-_ms|j-1> = sum_j groundstate_ps(j)groundstate_ps(k*)
! where the state |k*> is simply |j+1> with the spin Sz_ms flipped. The presence of c_gauge_ms depends if we flip a spin + or a spin -  
		do j=1,hd
			b = BTEST(j-1,i) 
			sum_Sz = sum_Sz + 0.5*bit_spin(b)*groundstate_ps(j)**2

! Now, O mix states with j and j+1,j-1, so we need to define the following k.		
			k = IBCHNG(j-1,i) + 1 
			
			if(b) then
				sum_Sp = sum_Sp + groundstate_ps(j)*groundstate_ps(k)
			end if 
		end do
		
! We save the <Sz> and <O> for every orbital and spin degree of freedom (ms).
		average_Sz_ms(i+1) = average_Sz_ms(i+1) + sum_Sz 
		average_Sp_ms(i+1) = average_Sp_ms(i+1) + sum_Sp
	end do

! We will take the medium value of <Sz_m> and <O_m> for each orbital m (we have an spin-degenerate system).
	do i=1,N_orb
		average_Sz(i) = average_Sz_ms(i)
		average_Sp(i) = average_Sp_ms(i)
	end do
end subroutine calculation_average_ps

