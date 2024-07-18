## Usage

pre-requisites: [`ASE`](https://wiki.fysik.dtu.dk/ase/)

### chgcent.py
Calculate charge electron and ion charge centers.
- Read:
  - CHGCAR
  - OUTCAR
- Output:
  - charge center
  - dipole moment [elect]
  - dipole moment [ionic]
  - total dipole moment
- Usage:
```
chgcent.py CHGCAR OUTCAR
```

### chgcent_cube.py
Calculate charge electron and ion charge centers. Suitable for cube files.
- Read:
  - filename.cube
  - ZVAL [list of integers, manual input]
  - ions_per_type [list of integers, manual input]
- Output:
  - charge center
  - dipole moment [elect]
  - dipole moment [ionic]
  - total dipole moment
- Usage:
```
chgcent.py filename.cube "1 2" "3 4"
```

### chgdiff.py
Calculate differential charge density. `CHGDIFF = CHG1 - CHG2`
- Read:
  - CHGCAR
- Output:
  - CHGDIFF.vasp
- Usage:
```
chgdiff.py CHGCAR1 CHGCAR2
```

### impose_sym.py
Find, print and impose symmetry to POSCAR.
- Read:
  - POSCAR
- Output:
  - POSCAR_0.1.vasp
- Usage:
```
impose_sym.py -c POSCAR -s 0.1
```

### p4vmod.py
Convert different naming conventions from v541 to v544 for `vasprun.xml`.
- Read:
  - vasprun.xml
- Output:
  - vasprun.541.xml
- Usage:
```
p4vmod.py vasprun.xml
```

### spincar.py
Calculate spin densities.
- Read:
  - CHGCAR
- Output:
  - SPINCAR*.vasp
- Usage:
```
spincar.py CHGCAR
```

### view_atoms.py
Quickly plot geometry with GUI from command line, light weight.
- Read:
  - POSCAR
- Output:
  - None
- Usage:
```
view_atoms.py POSCAR
```

### vtotav.py
Calculate planar average and macroscopic average for density files.
- Read:
  - LOCPOT
- Output:
  - LOCPOT_*
- Option:
  - `-macro --len DISTANCE` for macroscopic average
- Optional output:
  - `LOCPOT_*` interpolated planar average
  - `LOCPOT_*_macro` macroscopic average
- Usage:
```
vtotav.py -c LOCPOT -dir z -macro --len 5.0
```


### plotwfc
*use vaspwfc from [QijingZheng/VaspBandUnfolding](https://github.com/QijingZheng/VaspBandUnfolding/blob/master/vaspwfc.py).*
Plot Bloch wavefunctions.
- Read:
  - WAVECAR
  - POSCAR
- Output:
  - `WAV.*.vasp`
- Option:
  - `-soc` spinor wavefunction.
  - `-chg` output charge density, doesn't work with `-soc`
- Usage:
```
plotwfc -b 12 -k 1 -s 1 --soc
```

### cell_tool.py
Make supercells, measure bond length on the fly.
- Read:
  - POSCAR/CONTCAR
- Output:
  - POSCAR.super.vasp
  - distance between atoms.
- Usage:
```
cell_tool.py -d --dim '1 2 1' -v -di --dis '0 1'
cell_tool.py -d --dim '1 2 0 0 1 0 0 0 1' -v
```

### get_pp.py
Generate pseudopotential(POTCAR) file by combining different POTCAR files.
- Read:
  - POSCAR
- Output:
  - POTCAR
- Usage:
```
get_pp.py -c "POSCAR" -o -s "manual"
get_pp.py -i "Au Fe" -o -s "recommended" -xc "LDA"
```

### kp_gen.py
generate KPOINTS file.
- Read:
  - POSCAR
- Output:
  - KPOINTS_k_path
- Usage:
```
kp_gen.py -c POSCAR -r 0.1
```

### vasp_clean.py
clean current dir and leave only vasp input files.
- Read:
  - Current dir
- Output:
  - None
- Usage:
```
vasp_clean.py
vasp_clean.py -f
vasp_clean.py -a REPORT XDATCAR
vasp_clean.py -f -a REPORT XDATCAR
```

