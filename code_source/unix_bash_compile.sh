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
modules=''

if [ $compiler = 'ifort' ]
then
module add intel
mod="$mod"lomonosov_ifort/
p=" -fpconstant -r8 -openmp -fast -unroll-aggressive -module $mod -I $mod "
elif [ $compiler = 'gfortran' ]
then
module add gcc
gfortran -v
mod="$mod"lomonosov_gfortran/
#p=" -fdefault-real-8 -fopenmp -Ofast -funroll-all-loops -ffast-math -funsafe-math-optimizations -flto -march=native -Wsurprising -Wunused -J $mod -I $mod "
p=" -fdefault-real-8 -fopenmp -O3 -funroll-all-loops -flto -march=native -Wsurprising -Wunused -J $mod -I $mod "
$compiler $p -c "$md"IFPORT_illusion.f90 -o "$mod"IFPORT_illusion.o
modules="$modules $mod"'IFPORT_illusion.o '
fi

$compiler $p -c "$md"perfomance_settings.f90 -o "$mod"perfomance_settings.o
modules="$modules $mod"'perfomance_settings.o '
$compiler $p -c "$md"md_general.f90 -o "$mod"md_general.o
modules="$modules $mod"'md_general.o '
$compiler $p -c "$md"md_integrators.f90 -o "$mod"md_integrators.o
modules="$modules $mod"'md_integrators.o '
$compiler $p -c "$md"md_neighbours.f90 -o "$mod"md_neighbours.o
modules="$modules $mod"'md_neighbours.o '
$compiler $p -c "$md"md_read_write.f90 -o "$mod"md_read_write.o
modules="$modules $mod"'md_read_write.o '
$compiler $p -c "$int"cut_off_function.f90 -o "$mod"cut_off_function.o
modules="$modules $mod"'cut_off_function.o '
$compiler $p -c "$int"cut_off_poly.f90 -o "$mod"cut_off_poly.o
modules="$modules $mod"'cut_off_poly.o '
$compiler $p -c "$int"LennardJones.f90 -o "$mod"LennardJones.o
modules="$modules $mod"'LennardJones.o '
$compiler $p -c "$int"LennardJones_1g.f90 -o "$mod"LennardJones_1g.o
modules="$modules $mod"'LennardJones_1g.o '
$compiler $p -c "$int"graphenenorm.f90 -o "$mod"graphenenorm.o
modules="$modules $mod"'graphenenorm.o '
$compiler $p -c "$int"LennardJonesCosine.f90 -o "$mod"LennardJonesCosine.o
modules="$modules $mod"'LennardJonesCosine.o '
$compiler $p -c "$int"MorseCosine.f90 -o "$mod"MorseCosine.o
modules="$modules $mod"'MorseCosine.o '
$compiler $p -c "$int"TersoffBrenner.f90 -o "$mod"TersoffBrenner.o
modules="$modules $mod"'TersoffBrenner.o '
$compiler $p -c "$int"REBOsolidcarbon.f90 -o "$mod"REBOsolidcarbon.o
modules="$modules $mod"'REBOsolidcarbon.o '
$compiler $p -c "$int"RosatoGuillopeLegrand.f90 -o "$mod"RosatoGuillopeLegrand.o
modules="$modules $mod"'RosatoGuillopeLegrand.o '

$compiler $p -c "$md"md_interactions.f90 -o "$mod"md_interactions.o
modules="$modules $mod"'md_interactions.o '
$compiler $p -c "$md"md_simulation.f90 -o "$mod"md_simulation.o
modules="$modules $mod"'md_simulation.o '

$compiler $p -o "$exe"run_md_simulation_"$compiler" "$rns"run_md_simulation.f90 $modules

$compiler $p -c "$anl"graphene_on_surface_analysis.f90 -o "$mod"graphene_on_surface_analysis.o
modules="$modules $mod"'graphene_on_surface_analysis.o '
$compiler $p -o "$exe"run_gr_analysis_"$compiler" "$rns"run_gr_analysis.f90 "$mod"graphene_on_surface_analysis.o "$mod"md_general.o "$mod"md_read_write.o

$compiler $p -c "$lmf"fit_gr_moire.f90 -o "$mod"fit_gr_moire.o
modules="$modules $mod"'fit_gr_moire.o'
$compiler $p -o "$exe"run_gr_moire_fitting_"$compiler" "$rns"run_gr_moire_fitting.f90 $modules

module rm intel
module rm gcc