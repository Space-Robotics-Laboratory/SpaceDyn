%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	June 8, 1998, Last modification by K.Yoshida
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% CALC_GJ	Calculate the Generalized Jacobian.
%
%   CALC_GJ returns the generalized jacobian, GJ (6xn).
%
%   1998 (C)Space Robotics Lab, by Koichi Fujishima
%   2001.9.13 H.Hamano
%   2002.2.27 H.Hamano modified to version up
%

function [GJ,Addtherm] = calc_gj11( LP, SV, L, num_e, HH)

global Ez


% Calculate inertia matrices, HH
%HH = calc_hh( LP, SV );

% Find joint connection from the end-link to the 0-th link
joints = j_num( LP, num_e );

% calculate Jacobian and inertia matrices
Jm = calc_je( LP, SV, joints );

[ pe, tmp1 ] = f_kin_e( LP, SV, joints );
Js = [   eye(3,3) -tilde(pe-SV.R0);
       zeros(3,3)         eye(3,3) ];

Hs = HH(1:6,1:6);
Hm = HH(1:6,7:6+LP.num_q);

% Calculate the Generalized Jacobian
GJ = Jm - Js*inv(Hs)*Hm;
Addtherm = Js*inv(Hs)*L;


%%% EOF
