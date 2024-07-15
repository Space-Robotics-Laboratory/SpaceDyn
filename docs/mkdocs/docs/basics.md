# Basics for SpaceDyn


## Equation of motion (EOM) for floating-based system

To be added...

## Generalized Jacobian [2]

To be added...

!!! Note
    If you use Genralized Jacobian in an academic context, please put the following citation.
    
[2] Y. Umetani and K. Yoshida, Resolved motion rate control of space manipulators with generalized Jacobian matrix, _IEEE Transactions on Robotics and Automation_, vol. 5, no. 3, pp. 303--314, 1989.

    @ARTICLE{generalizedJacobian,
      author={Umetani, Y. and Yoshida, K.},
      journal={IEEE Transactions on Robotics and Automation}, 
      title={Resolved motion rate control of space manipulators with generalized Jacobian matrix}, 
      year={1989},
      volume={5},
      number={3},
      pages={303--314},
      doi={10.1109/70.34766}
    }

## Basic Knowledge for SpaceDyn User

* We assume the system composed of _n+1_ bodies and connected by _n_ joints. Let the body 0 be a _reference body_. Multiple branches can attach on any single _body_, as far as the system keeps a topological tree configuration. There must be a single _joint_ between two bodies. We call a terminal point or the point of interest such as manipulator hand as endpoint. Each body, except body 0, can have one endpoint at maximum. In this document, the terms _body_ and _link_ are the same.

* The SpaceDyn allows force/torque input on (i) the centroid of the reference body, (ii) each endpoint, and (iii) each joint. The toolbox computes the position, velocity and acceleration of (1) the centroid of the reference body, (2) the centroid of each body, (3) each endopoint, and (4) each joint.

* Computation of input force/torque are open to user programming. You can arbitrary decide each joint as either active or passive one. If you give always zero torque, such as $`\tau_0 = 0`$, the corresponding joint behaves as a free joint. Or if you give such a torque as: 
```math
\tau_i = - K q_i - D \dot{q}_i
```
the joint behaves as a passive visco-elastic joint. You can treat even a flexible link, by modeling it as a discrete successive chain of rigid links connected by elastic joints. Of course, you can give any arbitrary control torque determined by your own control law, on all or arbitrary selected joints. 

* We know that the Denavit-Hartenberg (DH) notation is commonly used in the field of manipulator kinematics with the advantage of unique allocation of coordinate systems with minimum parameters, but we know that the DH sometimes locates the coordinate ofitin away from the location of an actual joint. From the dynamics point of view, the angular velocity and the inertia tensor should be defined around the corresponding joint axis or body centroid. We then do NOT use the DH notation but introduce a rule to define the coordinate system with more flexibility. 
!!! Note
    We do NOT use the DH notation in SpaceDyn.
  Our rule locates the origin of the coordinate systems with more flexibility. Our rule locates the origin of the frame on each joint and orients the primary axes so that the inertia tensor should be simpler, but admits three position and three orientation parameters among two successive coordinate systems.

* For the representation of attitude or orientation, we use 3 by 3 direction cosine matrices. The advantage of direction cosine is (1) singularity free, (2) we can easily defive Roll-Pitch-Yaw angles, Euler angles, or quartanions, and (3) it is easy to find the mathematical relationship with angular velocity. 

* On the other hand, we frequently need Roll-Pitch-Yaw (RPY) replresentation also. For example, in order to express the twisting angles between two coordinate systems, we consider $`\alpha`$ (roll) around $`x`$-axis, $`\beta`$ (pitch) around $`y`$-axis, then $`\gamma`$ (yaw) around $`z`$-axis.

* Weak points: The SpaceDyn is not good at dealing with kinematic constraints other than joint axes. It is also weak at dealing with the problems in which a contact point is dynamically changing. For those problems, a good user programming is required to model the constraint forces. 

