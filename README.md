Info
====

A simulator of diffusion of Ornstein-Uhlenbeck-like diffusion.
Particles are elastically bound to the origin. The simulation is based on Smoluchowski dynamics (i.e. noise is fed into the momentum, not displacement as in Einstein view).

    dd x = dy = -k*x*dt + (2*D*dt)^(1/2)* dW
    dx = y * dt

The simulation is done in C++.
The ploting is in Python.

Technicalities
=============
The implementation in the `main.cpp` looks following:

    particles[tt][nn] = particles[tt-1][nn] + p->dt * velocity_particles[nn] ; 

    velocity_particles[nn] += p->dt * p->drift_term_velocity(particles[tt-1][nn], velocity_particles[nn]) \
                                                            + p->diffusion_term( dW );

The functions (`drift_term_velocity`, `diffusion_term`) are defined within the `Params.cpp` object file.
