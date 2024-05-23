%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    version 1.0 // Oct.4, 1999, Last modification by K.Yoshida
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%CZ	A Principal Coordinate Transformation Matrix, Dirction Cosines; C3
%
%	CZ(theta) returns a 3x3 transformation representing a 
%	rotation of theta about the Z axis.
%
%	See also CX, CY.

function Cz = cz( theta )


Cz = [  cos(theta)  sin(theta)  0;
       -sin(theta)  cos(theta)  0;
                 0           0  1 ];

%%%EOF
