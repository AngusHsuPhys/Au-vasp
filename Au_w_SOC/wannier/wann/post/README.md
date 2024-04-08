First we need to take the OUTCAR from bandstructure calculation. Find the KPOINTS in reciprocal lattice.
Extract the information with the following command:
```
awk 'NR < 1234 { next } { print } NR == 1449 { exit }' OUTCAR > kpointsweneed
```
where the ```1234``` and ```1449``` need to be replace by the first and last position of your KPOINTS.
Convert the 

Take the wannier90.chk, and wannier90.eig, from calculation in wann. Run the post Wannier90 to get wannier90_geninterp.dat.
