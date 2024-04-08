In order to perform the bands calculation, make sure the INCAR, and KPOINTS setting are correct.

For SOC calculation, I will turn off the `ISMEAR = 0`, `ICHARG=11` to read `CHGCAR` from SCF calculation. Linked them by `ln -s ../CHGCAR ./`.

Specific setting for SOC bands calculation. `MAGMOM = 3*IONS*0`, 3 for $x$, $y$, and $z$ direction. `LSORBIT =.TRUE.`. Executable file in batchrun has to be `vasp_ncl`.
