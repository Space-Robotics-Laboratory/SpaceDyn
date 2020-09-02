%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    version 1.0 // Oct.4, 1999, Last modification by K.Yoshida
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%CY	A Principal Coordinate Transformation Matrix, Dirction Cosines; C2
%
%	CY(theta) returns a 3x3 transformation representing a 
%	rotation of theta about the Y axis.
%
%	See also CX, CZ.
%

function Cy = cy(theta)


Cy = [ cos(theta)  0  -sin(theta);
                0  1            0;
       sin(theta)  0   cos(theta) ];


%%%EOF
