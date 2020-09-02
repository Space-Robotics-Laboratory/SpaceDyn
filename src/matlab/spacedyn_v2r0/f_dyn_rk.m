%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    version 1 // Oct 4, 1999, Last modification by K.Yoshida
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% F_DYN_RK	Compute the Forward Dynamics of the system
%		with Integration by Runge-Kutta method
%		(The fomula of delta_CO = tilde(w0')*C0 is used.
%		But its integration is subject to error accumulation,
%		yielding non-normal, non-orthogonal C0.
%		f_dyn_rk2 is recommended for large attitude motion.) 
%
%		A constant time step is used.
%
%		[R0,A0,v0,w0,q,qd]
%		   = f_dyn_rk(R0,A0,q,v0,w0,q,qd,F0,T0,Fe,Te,tau)
%
%		NOTE: v0,w0,vd0,wd0 are in the inertia frame.
%
%   2002.2.27 H.Nakanishi modified to version up
%

function SV = f_dyn_rk( LP, SV )

global Ez Gravity d_time

% Definition

C0 = SV.A0';

% 1st Step

SV = f_dyn(LP, SV);

k1_R0 = d_time * SV.v0;
k1_C0 = d_time * tilde(SV.w0)'*C0;
k1_A0 = k1_C0';
k1_q  = d_time * SV.qd;

k1_v0 = d_time * SV.vd0;
k1_w0 = d_time * SV.wd0;
k1_qd = d_time * SV.qdd;

SV.R0 = SV.R0 + k1_R0/2;
SV.A0 = SV.A0 + k1_A0/2;
SV.v0 = SV.v0 + k1_v0/2;
SV.w0 = SV.w0 + k1_w0/2;
SV.q = SV.q + k1_q/2;
SV.qd = SV.qd + k1_qd/2;

% 2nd Step

SV = f_dyn(LP, SV);

k2_R0 = d_time * ( SV.v0 + k1_v0/2 );
k2_C0 = d_time * ( tilde( SV.w0 + k1_w0/2 )' * C0 );
k2_A0 = k2_C0';
k2_q  = d_time * ( SV.qd + k1_qd/2 );

k2_v0 = d_time * SV.vd0;
k2_w0 = d_time * SV.wd0;
k2_qd = d_time * SV.qdd;

SV.R0 = SV.R0 + k2_R0/2;
SV.A0 = SV.A0 + k2_A0/2;
SV.v0 = SV.v0 + k2_v0/2;
SV.w0 = SV.w0 + k2_w0/2;
SV.q = SV.q + k2_q/2;
SV.qd = SV.qd + k2_qd/2;

% 3rd Step

SV = f_dyn(LP, SV);

k3_R0 = d_time * ( SV.v0 + k2_v0/2 );
k3_C0 = d_time * ( tilde( SV.w0 + k2_w0/2 )'* C0 );
k3_A0 = k3_C0';
k3_q  = d_time * ( SV.qd + k2_qd/2 );

k3_v0 = d_time * SV.vd0;
k3_w0 = d_time * SV.wd0;
k3_qd = d_time * SV.qdd;

SV.R0 = SV.R0 + k3_R0;
SV.A0 = SV.A0 + k3_A0;
SV.v0 = SV.v0 + k3_v0;
SV.w0 = SV.w0 + k3_w0;
SV.q = SV.q + k3_q;
SV.qd = SV.qd + k3_qd;


% 4th Step

SV = f_dyn(LP, SV);
k4_R0 = d_time * ( SV.v0 + k3_v0 );
k4_C0 = d_time * ( tilde( SV.w0 + k3_w0 )'* C0 );
k4_A0 = k4_C0';
k4_q  = d_time * ( SV.qd + k3_qd );

k4_v0 = d_time * SV.vd0;
k4_w0 = d_time * SV.wd0;
k4_qd = d_time * SV.qdd;

% Solution

SV.R0 = SV.R0 + ( k1_R0 + 2*k2_R0 + 2*k3_R0 + k4_R0 )/6;
C0_next = C0 + ( k1_C0 + 2*k2_C0 + 2*k3_C0 + k4_C0 )/6;
SV.A0 = C0_next';
SV.q  = SV.q + ( k1_q + 2*k2_q + 2*k3_q + k4_q )/6;
SV.v0 = SV.v0 + ( k1_v0 + 2*k2_v0 + 2*k3_v0 + k4_v0 )/6;
SV.w0 = SV.w0 + ( k1_w0 + 2*k2_w0 + 2*k3_w0 + k4_w0 )/6;
SV.qd = SV.qd + ( k1_qd + 2*k2_qd + 2*k3_qd + k4_qd )/6;

%%% EOF
