SET md=MOLECULAR_DYNAMICS\
SET int=INTERACTION_POTENTIALS\
SET rns=runners\
SET anl=graphene_on_surface_analysis\
SET lmf=ljc_and_morsec_moire_graphene_fitting\
SET mod=mods\my_pc_gfortran\
SET exe=..\executables\

SET compiler=gfortran -fdefault-real-8 -fopenmp -O3 -Wsurprising -Wunused -J %mod% -I %mod%

%compiler% -c %md%md_general_gfortran.f90 -o %mod%md_general.o

%compiler% -c %int%cut_off_function.f90 -o %mod%cut_off_function.o
%compiler% -c %int%LennardJones.f90 -o %mod%LennardJones.o
%compiler% -c %int%graphenenorm.f90 -o %mod%graphenenorm.o
%compiler% -c %int%LennardJonesCosine.f90 -o %mod%LennardJonesCosine.o
%compiler% -c %int%MorseCosine.f90 -o %mod%MorseCosine.o
%compiler% -c %int%TersoffBrenner.f90 -o %mod%TersoffBrenner.o
%compiler% -c %int%REBOsolidcarbon.f90 -o %mod%REBOsolidcarbon.o
%compiler% -c %int%RosatoGuillopeLegrand.f90 -o %mod%RosatoGuillopeLegrand.o

%compiler% -c %md%md_interactions.f90 -o %mod%md_interactions.o
%compiler% -c %md%md_simulation.f90 -o %mod%md_simulation.o

%compiler% -o %exe%run_md_simulation.exe %rns%run_md_simulation.f90 %mod%md_general.o %mod%md_simulation.o %mod%md_interactions.o %mod%graphenenorm.o %mod%cut_off_function.o %mod%LennardJones.o %mod%LennardJonesCosine.o %mod%MorseCosine.o %mod%RosatoGuillopeLegrand.o %mod%TersoffBrenner.o %mod%REBOsolidcarbon.o

%compiler% -c %anl%graphene_on_surface_analysis.f90 -o %mod%graphene_on_surface_analysis.o
%compiler% -o %exe%run_gr_analysis.exe %rns%run_gr_analysis.f90 %mod%graphene_on_surface_analysis.o %mod%md_general.o

%compiler% -c %lmf%fit_gr_moire.f90 -o %mod%fit_gr_moire.o
%compiler% -o %exe%run_gr_moire_fitting.exe %rns%run_gr_moire_fitting.f90 %mod%fit_gr_moire.o %mod%md_general.o %mod%md_simulation.o %mod%md_interactions.o %mod%graphenenorm.o %mod%cut_off_function.o %mod%LennardJones.o %mod%LennardJonesCosine.o %mod%MorseCosine.o %mod%RosatoGuillopeLegrand.o %mod%TersoffBrenner.o %mod%REBOsolidcarbon.o %mod%graphene_on_surface_analysis.o

pause