%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	May 3, 1998, Last modification by K.Yoshida
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% TR2DIFF	Calculate difference between two transformation matrices
%		in the form of role/pitch/yaw angles.
%		This computation is true only for small angles,
%		and falls on singularity at the angle pi (180 deg).
%
% 1997.11.12 (C)Space Robotics Lab, by Koichi Fujishima
%

function diff = tr2diff( AA1, AA2 )

diff = [ cross(AA2(1:3,1), AA1(1:3,1)) + ...
         cross(AA2(1:3,2), AA1(1:3,2)) + ...
         cross(AA2(1:3,3), AA1(1:3,3)) ...
      ]/2;

%%%EOF
