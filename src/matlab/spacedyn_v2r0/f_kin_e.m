%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	February 4, 1999, Last modification by K.Yoshida
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% F_KIN_E	Foward Kinematics
%		Caluculate the Position/Orientation of the point
%		specified by connection vector 'joints.'
%		POS : 3x1 vector, ORI : 3x3 matrix.
%
%   (C)Space Robotics Lab
%   1998.1.12 A.Kurosu
%   2001.9.11 H.Hamano
%   2002.2.27 H.Nakanishi modified to version up
%

function [ POS_e , ORI_e ] = f_kin_e(LP, SV, joints)


% Check number of the corresponding joints
n = length(joints);
k = joints(n);

% Calcurate coordinate trasformation matrix of Effector
A_I_i = SV.AA(:,k*3-2:k*3);
A_i_EE = rpy2dc(LP.Qe(:,k))';
ORI_e = A_I_i*A_i_EE;

% Calculate position vector of Effector
POS_e = SV.RR(:,k) + A_I_i * LP.ce(:,k);


%%%EOF

