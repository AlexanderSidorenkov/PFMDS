# PFMDS

## Parallel Fortran Molecular Dynamics Simulation

"Real Programmers use FORTRAN."

"If you can’t do it in FORTRAN, do it in assembly language. If you can’t do it in assembly language, it isn’t worth doing."

(Ed Post, Datamation Volume 29 number 7 July 1983)

## Что это?

Это набор программ для молекулярной динамики (МД), подгона потенциалов и анализа результатов.

Поддерживается OpenMP параллелизм. Несколько независимых прогонов МД могут быть запущены параллельно с помощью MPI.

## Как скачать?

Исходники и скомпилированные .exe для Windows можно скачать по ссылке [стабильная версия](https://github.com/AlexanderSidorenkov/PFMDS/archive/stable.zip). Или здесь: https://github.com/AlexanderSidorenkov/PFMDS/tree/stable , зеленая кнопка Clone or Download. 

## Документация

Скачиваем описание [отсюда](https://github.com/AlexanderSidorenkov/PFMDS_manual/archive/master.zip). Разархивируем и жмем на documentation. Просмотр с помощью браузера.

## Список программ:

### run_md_simulation и run_md_simulation_mpi

Моделирование методом МД.

Потенциалы:

* Леннарда-Джонса (между двумя группами атомов и внутри одной) - lj, lj1g
* модифицированный Леннарда-Джонса c косинусом (графен - металл) - ljc
* модифицированный Морзе c косинусом (графен - металл) - morsec
* Розато-Жиллопа-Легранда для металлов (один тип атомов) - rjl
* Терсоффа-Бреннера (углерод) - tb
* REBO (углерод, только подсчет энергии, сил нет) - rebosc

Интеграторы:

* Верле в скоростной форме - nve
* Верле в скоростной форме с цепочкой термостатов Нозе-Гувера (пока возможно применение только одной цепочки) - nvt
* Молекулярная статика - nvms

### run_gr_moire_fitting

Подгон параметров модифицированных потенциалов Леннарда-Джонса или Морзе.

### run_gr_analysis

Анализ расположения графена на поверхности меди.

## Дополнительный софт

### Компилятор

* ifort

Работает без проблем (проверено на кластере Ломоносов-1), но с gfortran программа получается быстрее, это вопрос оптимизации.

* gfortran 

В gfortran версии 4.4 есть баг при аллокации omp private массивов в параллельной части программы. Рекомендуется использовать более поздние версии (проверено с gcc version 6.3.0).

### Просмотр .xyz файлов

Для информации о атомах в вычислительной ячейки испльзован формат .xyz (http://libatoms.github.io/QUIP/io.html#module-ase.io.extxyz). Просматривать можно, например, с помощью [OVITO](https://ovito.org).

## Как компилировать?

На кластере Ломоносов запустить скрипт unix_bash_compile.sh в папке code_source.
На Windows запустить скрипт windows_compile.cmd в папке code_source.

Изменить испльзуемый компилятор, путь к нему, опции компиляции и т. д. можно внутри этих скриптов.

Исполняемые файлы появятся в папке executables.

## Как запуситить моделирование?

Исполняемые файлы в папке executables - результат компиляции программ в папке code_source/runnes, в которых вызываются необходимые подпрограммы. Для запуска моделирования достаточно использовать программу run_md_simulation или run_md_simulation_mpi.

### Запуск run_md_simulation:

Запускает программу МД на одном процессоре (вычислительном узле). Если указан список файлов настроек, то будет по очереди запущено несколько прогонов.

Опции run_md_simulation. В скобках указаны значения по умолчанию.

-i,--input (md_run_settings.txt) Имя файла настроек.

-p,--prefix () Префикс выходных файлов.

-op,--out_period (1) Период вывода в шагах МД.

-omp_n,--openmp_threads_num (максимальное количество потоков) Количество OpenMP потоков.

-ipath,--input_path () Путь к папке с входными файлами

-ilist,--input_list () Имя файла со списком файлов настроек.

-opath,--out_path () Путь к папке с выходными файлами.

-ofile,--all_out_file (all_out.txt) Имя файла для общего вывода нескольких прогонов.

### Запуск run_md_simulation_mpi:

Запускает программу МД на нескольких процессорах (вычислительных узлах). Если указан список файлов настроек, то несколько прогонов будут распределены между доступными узлами (количество прогонов должно быть равно или больше количества вычислительных узлов). Если список не указан, то один прогон будет независимо запущен на всех узлах.

run_md_simulation_mpi имеет те же опции, что и run_md_simulation

### Входные файлы

#### Файл настроек

Имеет фиксированный порядок строк и переменных в строках, который определяется порядком чтения в md().

Названия параметров в файле настроек не проверяются при чтении. Единственное что важно, это порядок! 

Пример:

Двухкомпонентный Л-Дж газ из атомов A и B.

Сначала проводится моделирование при постоянной температуре, потом при постоянной энергии, потом все охаждается до 0 К молекулярной статикой.

```
md_step_limit:				1000000
logfilename:				md_run.log
init_xyz_filename:			AB.xyz
new_velocities:				T
zero_momentum_period:		1000000
particle_types_num:			2
groups_num:					4
	1		A	B	
	2		A	#	
	3		B	#	
	4		#	#	
all_moving_atoms_group_num:		1
xyz_moving_atoms_group_num:		1
z_moving_atoms_group_num:		4
termo_atoms_group_num:		1
all_atoms_group_num:		1
traj_group_num:				3
period_traj:				200
change_group_num:			0
invert_z_vel:				F
integrators_num:			3
 name    dt         len     snap     log
 nvt     0.50000    10000 	500000    1000
 nve     0.50000    10000 	500000    1000
 nvms    0.50000    10000 	500000    1000
ms_de:						1.e-8
nhc_params:					100.	3	10000.
initial_temerature:			100.
interactions_num:			3
lj  	parameters_LJ_A-B.txt
2	3	8		7.5		20
3	2	64		7.5		20
lj1g	parameters_LJ_A-A.txt
2	2	64		7.5		20
lj1g	parameters_LJ_B-B.txt
3	3	8		7.5		20
```

md_step_limit - лимит шагов МД.

logfilename - имя файла для регулярного вывода энергий и температуры.

init_xyz_filename - имя .xyz файла с начальными расположениями и скоростями атомов.

new_velocities - T или F, присваивать ли атомам новые скорости. Если T, то атомы из группы all_moving_atoms_group_num будут иметь скорости, распределенные по Максвеллу для температуры initial_temerature.

zero_momentum_period - период обнуления скорости центра масс группы атомов all_atoms_group_num. При долгих временах моделирования могут накопиться ошибки. Угловые скорости не обнуляются.

particle_types_num - количество рассматриваемых типов атомов.

groups_num - количество групп атомов. Один атом может входить в несколько групп или не входить ни в одну. Скорее всего понадобится пустая группа.

После groups_num идет таблица с groups_num строками и particle_types_num+1 столбцами. Первый столбец - вспомогательный, рекомендуется писать в него номера строк, они же номера групп. Остальные столбцы - имена атомов.

all_moving_atoms_group_num - номер группы со всеми двигающимися атомами.

xyz_moving_atoms_group_num - номер группы с атомами, двигающимися во всех направлениях

z_moving_atoms_group_num - номер группы с атомами, двигающимися только по оси Z. Группы xyz_moving_atoms_group_num и z_moving_atoms_group_num не должны пересекаться! Термостат не учитывает уменьшение количества степеней свободы для этих атомов!

termo_atoms_group_num - термостатируемая группа атомов.

all_atoms_group_num - группа со всеми моделируемыми атомами.

traj_group_num - группа атомов для вывода положений атомов в один файл каждые period_traj шагов МД.

period_traj - период вывода положений атомов из группы traj_group_num.

change_group_num - количество групп с изменяемым количеством атомов (для добавления атомов в ячейку). См. ???????.

invert_z_vel - T или F, создать или нет абсолютно упругую стенку в верхней части ячейки (для предотвращения прохода атомов через Z границу ячейки)

integrators_num - количество интеграторов.

Далее идет вспомогательная строка.

Далее integrators_num строк с названием интегратора, шагом интегрирования, количеством шагов, периодом вывода положений и скоростей атомов в отдельный файл и периодом вывода в logfilename.

ms_de - разница энергии для остановки молекулярной статики.

nhc_params - параметры цепочки термостатов Нозе-Гувера: температура, количество звеньев и масса первого звена.

initial_temerature - начальная температура группы all_moving_atoms_group_num если в new_velocities указано T.

interactions_num - количество взаимодействий между атомами.

Далее идут interactions_num блоков с параметрами для каждого взаимодействия: название взаимодействия, название файла с параметрами потенциала; в селующих строках - параметры списков соседей. 

Список соседей имеет следующие параметры (по порядку): номер первой группы атомов, второй группы атомов, максимальное количество соседей из второй группы для атомов из первой, максимальный радиус включения в список соседей (должен быть меньше радиуса обрезания потенциала и как минимум в два раза меньше минимального размера ячейки), период обновления списка в шагах.

#### Файл .xyz с начальными положениями и скоростями атомов.

Его имя указывается в файле настроек (init_xyz_filename). Также содержит размеры вычислительной ячейки, массы атомов и их названия.

Пример:

4 атома A и 4 атома B с нулевыми скоростями в ячейке 16*16*16
```
          8
Lattice=" 16.0 0.0 0.0 0.0 16.0 0.0 0.0 0.0 16.0 " Properties=pos:R:3:vel:R:3:mass:R:1:species:S:1
        6.0000000        6.0000000        6.0000000        0.0000000        0.0000000        0.0000000       10.0000000    B
        6.0000000        6.0000000       10.0000000        0.0000000        0.0000000        0.0000000       10.0000000    B
        6.0000000       10.0000000        6.0000000        0.0000000        0.0000000        0.0000000       10.0000000    B
       10.0000000        6.0000000        6.0000000        0.0000000        0.0000000        0.0000000       10.0000000    B
        0.0000000        0.0000000        0.0000000        0.0000000        0.0000000        0.0000000        1.0000000    A
        4.0000000        0.0000000        0.0000000        0.0000000        0.0000000        0.0000000        1.0000000    A
        8.0000000        0.0000000        0.0000000        0.0000000        0.0000000        0.0000000        1.0000000    A
       12.0000000        0.0000000        0.0000000        0.0000000        0.0000000        0.0000000        1.0000000    A
```

Первая строка - количество атомов.

Вторая строка - размеры ячейки и названия столбцов.

Ячейка задается тремя векторами (в общем случае ячейка триклинная). Пока поддерживается только прямоугольные ячейки. Размер по X - первое число, по Y - пятое, по Z - девятое. Остальные числа - нули.

Типы столбцов используются следующие - Properties=pos:R:3:vel:R:3:mass:R:1:species:S:1. То есть три координаты, три проекции скорости, масса и имя атома (в таком порядке)

С третьей строки идет информация об атомах.

Для генерации ячеек написаны отдельные программы (https://github.com/AlexanderSidorenkov/xyz_create).

#### Файлы с параметрами потенциалов

См. описания потенциалов взаимодействия.

### Запуск

На Windows запускать программы из командной строки или подготовить .cmd файл. См. в примерах.

На кластере удобнее всего запускать скриптом оболочки. Например:
```
#!/bin/bash
EXEPATH='/mnt/msu/users/your_path_to/PFMDS/executables/'
EXE='run_md_simulation_mpi_gfortran'
module add gcc
module add openmpi/1.8.4-gcc
export I_MPI_PIN=off

prefix='test_'
sbatch -o "$prefix"out.txt -e "$prefix"err.txt -p regular4 -t 3-0 -N 48 ompi --bind-to none "$EXEPATH""$EXE" -p "$prefix" -op 50000 -ipath input_files/ -opath results/ -ilist settings_list.txt &>"$prefix"job_id.txt

module rm gcc
module rm openmpi/1.8.4-gcc
```
Запускать из папки с папками input_files (там входные файлы) и results (создать пустую для рзультатов). Приведенный скрипт поставит в очередь задачу run_md_simulation_mpi_gfortran на 48 узлов, выходные файлы будут иметь префикс "test_".

### Выходная информация

Информация в основном потоке вывода.

Лог файл (logfilename).

Конечные значения энергий в общий файл вывода.

Конечные положения и скорости атомов.

Промежуточные положения и скорости атомов, выводимые с периодом указанным в параметрах интегратора.

Промежуточные положения и скорости атомов определенной группы, выводимые в один файл.

## Примеры входных и выходных файлов

Примеры находятся в папке examples в репозитории [документации](https://github.com/AlexanderSidorenkov/PFMDS_manual).

## Контакты

e-mail av.sidorenkov@physics.msu.ru

github https://github.com/AlexanderSidorenkov
