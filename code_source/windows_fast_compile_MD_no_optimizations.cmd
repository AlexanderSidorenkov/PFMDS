@echo off
SET md=MOLECULAR_DYNAMICS\
SET int=INTERACTION_POTENTIALS\
SET rns=runners\
SET anl=graphene_on_surface_analysis\
SET lmf=ljc_and_morsec_moire_graphene_fitting\
SET mod=mods\my_pc_gfortran\
SET exe=..\executables\

( date /t & time /t & gfortran -v 
) 1>compilation.txt 2>&1
SET optall=-fopt-info-optall=compilation_optall.txt
SET missed=-fopt-info-missed=compilation_missed.txt

SET compiler=gfortran -fdefault-real-8 -fopenmp -Wsurprising -Wunused -J %mod% -I %mod%
SET modules=

%compiler% -c %md%perfomance_settings.f90 -o %mod%perfomance_settings.o
SET modules=%modules% %mod%perfomance_settings.o
%compiler% -c %md%IFPORT_illusion.f90 -o %mod%IFPORT_illusion.o
SET modules=%modules% %mod%IFPORT_illusion.o
%compiler% -c %md%md_general.f90 -o %mod%md_general.o
SET modules=%modules% %mod%md_general.o
%compiler% -c %md%md_integrators.f90 -o %mod%md_integrators.o
SET modules=%modules% %mod%md_integrators.o
%compiler% -c %md%md_neighbours.f90 -o %mod%md_neighbours.o
SET modules=%modules% %mod%md_neighbours.o
%compiler% -c %md%md_read_write.f90 -o %mod%md_read_write.o
SET modules=%modules% %mod%md_read_write.o
%compiler% -c %int%cut_off_function.f90 -o %mod%cut_off_function.o
SET modules=%modules% %mod%cut_off_function.o
%compiler% -c %int%cut_off_poly.f90 -o %mod%cut_off_poly.o
SET modules=%modules% %mod%cut_off_poly.o
%compiler% -c %int%LennardJones.f90 -o %mod%LennardJones.o
SET modules=%modules% %mod%LennardJones.o
%compiler% -c %int%LennardJones_1g.f90 -o %mod%LennardJones_1g.o
SET modules=%modules% %mod%LennardJones_1g.o
%compiler% -c %int%graphenenorm.f90 -o %mod%graphenenorm.o
SET modules=%modules% %mod%graphenenorm.o
%compiler% -c %int%LennardJonesCosine.f90 -o %mod%LennardJonesCosine.o
SET modules=%modules% %mod%LennardJonesCosine.o
%compiler% -c %int%MorseCosine.f90 -o %mod%MorseCosine.o
SET modules=%modules% %mod%MorseCosine.o
%compiler% -c %int%TersoffBrenner.f90 -o %mod%TersoffBrenner.o
SET modules=%modules% %mod%TersoffBrenner.o
%compiler% -c %int%REBOsolidcarbon.f90 -o %mod%REBOsolidcarbon.o
SET modules=%modules% %mod%REBOsolidcarbon.o
%compiler% -c %int%RosatoGuillopeLegrand.f90 -o %mod%RosatoGuillopeLegrand.o
SET modules=%modules% %mod%RosatoGuillopeLegrand.o

%compiler% -c %md%md_interactions.f90 -o %mod%md_interactions.o
SET modules=%modules% %mod%md_interactions.o
%compiler% -c %md%md_simulation.f90 -o %mod%md_simulation.o
SET modules=%modules% %mod%md_simulation.o

%compiler% -o %exe%run_md_simulation.exe %rns%run_md_simulation.f90 %modules%

pause