## TODO List

Work needs to be done
- reference tracking using V(x) = (x-x_ref)'Q(x-x_ref) + u'Ru instead of B_ref

Issues to be considered/addressed:
- Alternative control scenarios (e.g. constant velocity pendulum 
  equilibrium),
- Disturbance rejection of e.g. sinusoidal input
- Insight into physically realistic actuator limitations (roll, pitch and 
  yaw rate limits etc based on estimated rotational and translational 
  inertia of the quadrotor, rotor thrust curves, blade dynamics etc).
- A Lyapunov analysis of the region of attraction/stability to get an idea
  of what the actuation limits should be (since large actuation might cause
  the system to move to states where the nonlinearities are non-neglible).

The experimental parameters of the quadrotor as used by the team at ETH 
Zurich can be found in literature directory of this repository 
