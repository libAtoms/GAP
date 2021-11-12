#  Hybrid MD decision making package
#
#  Copyright (c) Tamas K. Stenczel 2021.

#rm -v *.stderr *.stdout sic_md.hybrid-md.xyz std*.txt *.xyz.idx debug_output.txt *.usp sic_md.md sic_md.*err sic_md.castep train.xyz sic_md.hybrid-md-state.yaml
mpirun -n 16 castep.mpi sic_md