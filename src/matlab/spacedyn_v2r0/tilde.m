%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	May 3, 1998, Last modification by K.Yoshida
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%	tilde.m : a skew-symmetric operator  Å@%%Å@scew?
%
%	B=A~
%
%	Copyright (C) Motiaki Shimizu 1996/10/9
%

function B = tilde(A)


B = [  0    -A(3)    A(2);
      A(3)    0     -A(1);
     -A(2)   A(1)     0   ];


%%%EOF
