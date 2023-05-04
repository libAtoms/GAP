# Simple test job to test compilation

8-atom SiC cell, running 50 steps of MD.

This is intended to check that your installation is correct and the needed executables are in your path. There is no
scientific meaning of this test, the basis set and all settings are dumbed down a lot.

However, this even runs on a laptop in serial.

*n.b. If you are on Castep v22, then use the alternative param file `sic_md.param-v22` instead, which is the same in
essence but has the relevant keywords.*

## how to execute

You should have the following in your path:

- `gap_fit`: compiled serial or with OpenMP (no MPI linked to it)
- `hybrid-md`: python package installs this executable, you need to activate the correct Python env to have it

```bash
mpirun -n 2 castep.mpi sic_md
```

## Success criteria

Pretty much that no error is seen with the above, plus a number of files showing up:

- `GAP_model_step000.../` directories showing up - these are the model versions retired at given steps of the run
- `GAP.xml` + sparse files: this is the latest model
- `sic_md.castep` castep output file
- `sic_md.hybrid-md.xyz` ab-initio observations gathered during the calculation

