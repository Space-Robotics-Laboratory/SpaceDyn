%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	February 4, 1999, Last modification by K.Yoshida
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% F_KINE_J	Forward Kinematics:
%		Caluculate the Position/Orientation of the joints
%		corresponding the point specified by connection vector 'joints.'
%		POS : 3x1 vector, ORI : 3x3 matrix,
%		POS_j : 3xn, ORI : 3x3n,
%		where n : number of joints between the End-point to link 0.
%
%   (C)Space Robotics Lab
%   1998 1.9 A.Kurosu
%   2001.9.11 H.Hamano
%   2002.2.27 H.Nakanishi modified to version up
%

function [ POS_j , ORI_j ] = f_kin_j( LP, SV, joints )

global Ez


% Check the number of the corresponding joints
n = length( joints );


% Calculation of Orientation and Position of each joints
POS_j = [];
ORI_j = [];

for i = 1 : 1 : n
   
   PorR = ( LP.J_type(joints(i)) == 'P' );
   ORI_tmp = SV.AA(:,joints(i)*3-2:joints(i)*3);
   POS_tmp = SV.RR(:,joints(i)) + ORI_tmp * ( LP.cc(:,joints(i),joints(i)) - PorR*Ez * SV.q(joints(i)) );
   
   POS_j = [POS_j POS_tmp];
   ORI_j = [ORI_j ORI_tmp];
   
end


%%%EOF
