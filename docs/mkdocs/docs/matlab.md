# SpaceDyn as a MATLAB toolbox

<!-- 
SpaceDyn was originally created as MATLAB toolbox to handle the kinematic and dynamic analysis and simulation of articulated multi-body system with a free-floating base. Example of such systems include a satellite with mechanical appendages, a free-flying space robot, a wheeled mobile robot, and a walking robot, all of which makes motions in the environment with of without gravity. 

This toolbox can handle open chain systems with topological tree configuration. A parallel manipulator, for example, then cannot be supported directly. A walking robot contacting on the ground with more than two legs or limbs at a time seems to form a closed chain including the ground, however, we can handle such a system with a proper model of ground contact at each contact point. Parallel manipulators can be treated with virtual cut of a kinematic chain and corresponding virtual force model. 

Some academic papers regarding this toolbox is published by Kazuya Yoshida and his co-authors. For the technical points of this software, please consult with those publications as well as the following chapters of this document. 

We hope that you could find this toolbox useful for your application. 
-->

SpaceDyn was originally created as MATLAB toolbox [1]. We developed this toolbox motivated and inspired by [Robotics toolbox](https://petercorke.com/toolboxes/robotics-toolbox/) developed by Peter I. Corke. We took one m-file (```cross.m```) and use it as the original is, but SpaceDyn as a whole, does not have compatibility with the Peter's Robotics toolbox unfortunately.

!!! Note
    The MATLAB toolbox of Spacedyn requires MATLAB 5.0 or higher. (It uses functions that are not supported in 4.0 or lower.)

## Usage

Just clone the repository, and copy all the MATLAB functions in ```SpaceDyn/src/matlab``` to your workspace, and use them to model your robot and calculate its kinematics and dynamics. Each-by-each m-file usage is detailed in the following **quick reference** shown later.

## Technical Note

* For mathematical symbols, we give the name of variables with more than two letters. For example, vector _r_ and _c_ are coded by ```RR``` and ```cc```, respectively This is to avoid the confusion with control variables such as _i_, _j_, _k_, or _l_, _m_, _n_, etc, which are frequently used as iteration or array counters. However, there is an exeption: the symbol ```q``` is used for the joint variable vector _q_.

* Since MATLAB doesn't allow 0 for array index, we use ```R0``` and ```c0``` instead of ```RR[0]``` and ```cc[0]```, for example. 

* For the representation of attitude or orientation, we use 3 by 3 direction cosine matrices, coded with a symbol ```A```. For example, ```A0``` is the direction cosines to represent the attitude of the body 0. For the other bodies, a matrix ```AA``` is used. 

* For Roll-Pitch-Yaw (RPY) angles, we use the symbol ```Q```. For example, in order to express the twisting angles between two coordinate systems, we consider $`\alpha`$ (roll) around $`x`$-axis, $`\beta`$ (pitch) around $`y`$-axis, then $`\gumma`$ (yaw) around $`z`$-axis. The set of these angles are coded by ```Qi```.

* We use both the input variables and the global variables to pass the values to m-file functions. The input variables inside the braces are the variables changing time to time, such as joint angles, positions, orientations, and so on. The global variables are the ones holding constant once the model is given, such as topological description matrices, kinematics and dynamic parameters. 

<!-- * We assume the system composed of _n+1_ bodies and connected by _n_ joints. Let the body 0 be a _reference body_. Multiple branches can attach on any single _body_, as far as the system keeps a topological tree configuration. There must be a single _joint_ between two bodies. We call a terminal point or the point of interest such as manipulator hand as endpoint. Each body, except body 0, can have one endpoint at maximum. In this document, the terms _body_ and _link_ are the same.--> 

<!-- * The toolbox allows force/torque input on (i) the centroid of the reference body, (ii) each endpoint, and (iii) each joint. The toolbox computes the position, velocity and acceleration of (1) the centroid of the reference body, (2) the centroid of each body, (3) each endopoint, and (4) each joint. -->

## Quick Reference

The interactive quick reference guide written in html is already available in our [GitHub directory](https://github.com/Space-Robotics-Laboratory/SpaceDyn/tree/main/src/matlab/spacedyn_reference). 

Please download the repository ```spacedyn_reference``` and open ```index.html``` in web browser to get started. 

## Original User Manual

The original user manual is also avaliable.

#### [**The Spacedyn - a MATLAB Toolbox for Space and Mobile Robots**](http://www.astro.mech.tohoku.ac.jp/spacedyn/doc.pdf)

We heavily recommend the user to read this manual, which describe the development philosophy and programing rules, to understand the code, deeply.
