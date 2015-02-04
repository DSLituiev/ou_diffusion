a simulator of diffusion of Ornstein-Uhlenbeck-like diffusion.
Particles are elastically bound to the origin. The simulation is based on Smoluchowski dynamics (i.e. noise is fed into the momentum, not displacement as in Einstein view).

    dd x = dy = -k*x*dt + (2*D*dt)^(1/2)* dW
    dx = y * dt

The simulation is done in C++.
The ploting is in Python.
