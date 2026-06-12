----- Dipolar-Depletion Colloids -----

Simulation and Analysis scripts for simulating Dipolar-Depletion Colloids.

----- Description of Current Working Stuff -----

Nothing Yet ...

----- Description of Things to Add and Change -----

Add Lammps Scripts

Add gr, sq, and msd scripts
    fix sq script to work in the same way as gr and msd

Update Lammps script to:

    print thermal data along with dump files, logarithmically

    set up to read and output restart files

Add Bash Scripts to:

    Automatically create grid of statepoints (phix/epsy). For a given mu value
        Start at eps0 ---> run sim ---> generate restart file to read in by next directory input
        In each statepoint, there will be an input, and output subdir, in output theres configs and thermal data

    An analysis script that computes msd, gr, sq (when ready) for the statepoints

Add general Plotting scheme
    A python file with functions to make a plot (from matplotlib.pyplot)