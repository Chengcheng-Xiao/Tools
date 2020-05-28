## Usage
### chgcent.py
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

### chgdiff.py
- Read:
  - CHGCAR
- Output:
  - CHGDIFF.vasp
- Usage:
```
chgdiff.py CHGCAR1 CHGCAR2
```

### impose_sym.py
- Read:
  - POSCAR
- Output:
  - POSCAR_0.1.vasp
- Usage:
```
impose_sym.py -c POSCAR -s 0.1
```

### p4vmod.py
- Read:
  - vasprun.xml
- Output:
  - vasprun.541.xml
- Usage:
```
p4vmod.py vasprun.xml
```

### spincar.py
- Read:
  - CHGCAR
- Output:
  - SPINCAR*.vasp
- Usage:
```
spincar.py CHGCAR
```

### view_atoms.py
- Read:
  - POSCAR
- Output:
  - None
- Usage:
```
view_atoms.py POSCAR
```

### vtotav.py
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
- Read:
  - POSCAR
- Output:
  - POTCAR
- Usage:
```
get_pp.py -c "POSCAR" -o -s "manual"

get_pp.py -i "Au Fe" -o -s "recommended" -xc "LDA"

```
