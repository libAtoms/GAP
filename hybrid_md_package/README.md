# Hybrid MD 

This is a python pakage for running "Hybrid MD": accelerate MD in Quantum Mechanics (QM) codes with
Force Fields (FF) and allow on the fly refinement as well.

So far this has been used with GAP, particularly in the CASTEP code.

The QM codes need the FF interfaced for calculation of energies and forces, which 
is independent of this package, and need to allow switching between the QM and FF forces for 
the MD steps. Preferably, this switching should not be costly an no expensive QM calculations 
should be done in cases when the FF is used.

At specified points in the MD loop, the QM code needs to make decisions, for which this package 
was written. Given a simple interfce, this package allows rapid development of methods and 
tayloring those for specific questions or problems, all without touching the source code of the 
QM package used.

The calls to this package can be done through a system call to the installed executable of 
`hybrid-md` which has subcommand for :
- `initialise`: before the start of MD loop, decides to start with QM or FF
- `pre-step`: at the start of the MD step, decides if QM calculation is required and if any 
  comparison will be done later.
- `post-step`: after QM calculation, but before first Velocity-Verlet step. This performs the
  comparison, decides if there is a need for refitting and tells the QM code if the FF
  parameters need updating.
  
The subcommands need the file name prefix used in for the calculation's inputs and outputs,
called seed here, and the MD iteration number (apart from at initialisation). For full details see
the `--help` pages of the executable after installation or read the code.

The QM code should write it's results into an XYZ file for the executable to read if needed,
and the settings are stored in a yaml file `<seed>.hybrid-md-input.yaml`.

For further details, see the code and feel free to tweak it to your needs.
  
## Installation

Download the source code and run

```bash
pip install .
```