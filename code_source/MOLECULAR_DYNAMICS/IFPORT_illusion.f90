!> \brief Заглушка для удобства компиляции.
!> \detailed В IFPORT находится функция rand() при компиляции с ifort. При компиляции с gfortran такого модуля нет. Этот пустой модуль нужен чтобы не убирать use IFPORT в md_general при компиляции с gfortran.
module IFPORT
implicit none

contains

!nothing
!compile this for gfortran to leave "use IFPORT" in md_general uncommented
!do not compile it for ifort

end module IFPORT