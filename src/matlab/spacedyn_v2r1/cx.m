%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    version 1.0 // Oct.4, 1999, Last modification by K.Yoshida
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%CX	A Principal Coordinate Transformation Matrix, Dirction Cosines; C1
%
%	CX(theta) returns a 3x3 transformation representing a 
%	rotation of theta about the X axis.
%
%	See also CY, CZ.
%

function Cx = cx(theta)


Cx =[ 1           0          0;
      0  cos(theta) sin(theta);
      0 -sin(theta) cos(theta) ];


%%%EOF
