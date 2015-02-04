Info
====

A simulator of diffusion of Ornstein-Uhlenbeck-like diffusion.
Particles are elastically bound to the origin. The simulation is based on Smoluchowski dynamics (i.e. noise is fed into the momentum, not displacement as in Einstein view).

    dd x = dy = -k*x*dt + (2*D*dt)^(1/2)* dW
    dx = y * dt


Technicalities
=============
1. The simulation is done in C++. The ploting is in Python.

2. The implementation in the `main.cpp` looks following:

    particles[tt][nn] = particles[tt-1][nn] + p->dt * velocity_particles[nn] ; 

    velocity_particles[nn] += p->dt * p->drift_term_velocity(particles[tt-1][nn], velocity_particles[nn]) \
                                                            + p->diffusion_term( dW );

The functions (`drift_term_velocity`, `diffusion_term`) are defined within the `Params.cpp` object file.

2. Simulation results are store in an SQLite database. Each run should be assigned a label (under the `-l` flag) that is used as a table name. Within each run table, 
each entry corresponds to a time point.
Data particle positions are mapped on a grid, and counts within grid cells are stored as long integers within 'blob'-type fields in an SQLite database. Three raw particle coordinates are saved for visualization purposes.


