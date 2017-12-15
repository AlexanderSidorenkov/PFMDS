#!/bin/bash
md='MOLECULAR_DYNAMICS/'
int='INTERACTION_POTENTIALS/'
rns='runners/'
anl='graphene_on_surface_analysis/'
lmf='ljc_and_morsec_moire_graphene_fitting/'
mod='mods/'
exe='../executables/'

#compiler='ifort'
compiler='gfortran'

if [ $compiler = 'ifort' ]
then
module add intel
mod="$mod"lomonosov_ifort/
p=" -fpconstant -r8 -openmp -O3 -module $mod -I $mod "
$compiler -c "$md"md_general_ifort.f90 -o "$mod"md_general.o
elif [ $compiler = 'gfortran' ]
then
module add gcc
mod="$mod"lomonosov_gfortran/
p=" -fdefault-real-8 -fopenmp -O3 -Wsurprising -Wunused -J $mod -I $mod "
$compiler $p -c "$md"md_general_gfortran.f90 -o "$mod"md_general.o
fi

$compiler $p -c "$int"cut_off_function.f90 -o "$mod"cut_off_function.o
$compiler $p -c "$int"LennardJones.f90 -o "$mod"LennardJones.o
$compiler $p -c "$int"graphenenorm.f90 -o "$mod"graphenenorm.o
$compiler $p -c "$int"LennardJonesCosine.f90 -o "$mod"LennardJonesCosine.o
$compiler $p -c "$int"MorseCosine.f90 -o "$mod"MorseCosine.o
$compiler $p -c "$int"TersoffBrenner.f90 -o "$mod"TersoffBrenner.o
$compiler $p -c "$int"REBOsolidcarbon.f90 -o "$mod"REBOsolidcarbon.o
$compiler $p -c "$int"RosatoGuillopeLegrand.f90 -o "$mod"RosatoGuillopeLegrand.o

$compiler $p -c "$md"md_interactions.f90 -o "$mod"md_interactions.o
$compiler $p -c "$md"md_simulation.f90 -o "$mod"md_simulation.o

$compiler $p -o "$exe"run_md_simulation_"$compiler" "$rns"run_md_simulation.f90 "$mod"md_general.o "$mod"md_simulation.o "$mod"md_interactions.o "$mod"graphenenorm.o "$mod"cut_off_function.o "$mod"LennardJones.o "$mod"LennardJonesCosine.o "$mod"MorseCosine.o "$mod"RosatoGuillopeLegrand.o "$mod"TersoffBrenner.o "$mod"REBOsolidcarbon.o

$compiler $p -c "$anl"graphene_on_surface_analysis.f90 -o "$mod"graphene_on_surface_analysis.o
$compiler $p -o "$exe"run_gr_analysis_"$compiler" "$rns"run_gr_analysis.f90 "$mod"graphene_on_surface_analysis.o "$mod"md_general.o

$compiler $p -c "$lmf"fit_gr_moire.f90 -o "$mod"fit_gr_moire.o
$compiler $p -o "$exe"run_gr_moire_fitting_"$compiler" "$rns"run_gr_moire_fitting.f90 "$mod"fit_gr_moire.o "$mod"md_general.o "$mod"md_simulation.o "$mod"md_interactions.o "$mod"graphenenorm.o "$mod"cut_off_function.o "$mod"LennardJones.o "$mod"LennardJonesCosine.o "$mod"MorseCosine.o "$mod"RosatoGuillopeLegrand.o "$mod"TersoffBrenner.o "$mod"REBOsolidcarbon.o "$mod"graphene_on_surface_analysis.o

module rm intel
module rm gcc


