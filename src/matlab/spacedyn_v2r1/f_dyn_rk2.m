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
%
%		(The Rodorigues fomula for infinitesimal rotation
%		is used to update from A0 to A0_n.
%		This seems practically better than the fomula
%		delta_C0 = tilde(w0')*C0.)
%
%		A constant time step is used.
%
%		[R0,A0,v0,w0,q,qd]
%		   = f_dyn_rk2(R0,A0,q,v0,w0,q,qd,F0,T0,Fe,Te,tau)
%
%		NOTE: v0,w0,vd0,wd0 are in the inertia frame.
%
%       Edited by Haakon Fjelberg May 23, 1999
%
%   1998.9.2  written by K.FUJISHIMA
%   2002.2.27  H.Nakanishi modified to version up
%   2024.5.22  Updated by A.Uchida and M.Imai
%

function SV = f_dyn_rk2(LP, SV)

global Ez Gravity d_time


% 1st Step
SV.Qtn0 = dc2qtn(SV.A0');
SV.dQtn0 = w2dQtn(SV.w0, SV.Qtn0);
SVt = f_dyn(LP, SV);

k1_R0 = d_time * SV.v0;
k1_Qtn0 = d_time * SV.dQtn0;
k1_v0 = d_time * SVt.vd0;
k1_w0 = d_time * SVt.wd0;
k1_q  = d_time * SV.qd;
k1_qd = d_time * SVt.qdd;

SVt.R0 = SV.R0 + k1_R0/2;
SVt.Qtn0 = SV.Qtn0 + k1_Qtn0/2;
SVt.Qtn0 = SVt.Qtn0 / norm(SVt.Qtn0);
SVt.A0 = qtn2dc(SVt.Qtn0)';
SVt.v0 = SV.v0 + k1_v0/2;
SVt.w0 = SV.w0 + k1_w0/2;
SVt.q  = SV.q +  k1_q/2;
SVt.qd = SV.qd + k1_qd/2;


% 2nd Step
SVt.dQtn0 = w2dQtn(SVt.w0, SVt.Qtn0);
SVt = f_dyn(LP, SVt);

k2_R0 = d_time * SVt.v0;
k2_Qtn0 = d_time * SVt.dQtn0;
k2_v0 = d_time * SVt.vd0;
k2_w0 = d_time * SVt.wd0;
k2_q  = d_time * SVt.qd;
k2_qd = d_time * SVt.qdd;

SVt.R0 = SV.R0 + k2_R0/2;
SVt.Qtn0 = SV.Qtn0 + k2_Qtn0/2;
SVt.Qtn0 = SVt.Qtn0 / norm(SVt.Qtn0);
SVt.A0 = qtn2dc(SVt.Qtn0)';
SVt.v0 = SV.v0 + k2_v0/2;
SVt.w0 = SV.w0 + k2_w0/2;
SVt.q  = SV.q  + k2_q/2;
SVt.qd = SV.qd + k2_qd/2;

% 3rd Step
SVt.dQtn0 = w2dQtn(SVt.w0, SVt.Qtn0);
SVt = f_dyn(LP, SVt);

k3_R0 = d_time * SVt.v0;
k3_Qtn0 = d_time * SVt.dQtn0;
k3_v0 = d_time * SVt.vd0;
k3_w0 = d_time * SVt.wd0;
k3_q  = d_time * SVt.qd;
k3_qd = d_time * SVt.qdd;

SVt.R0 = SV.R0 + k3_R0;
SVt.Qtn0 = SV.Qtn0 + k3_Qtn0;
SVt.Qtn0 = SVt.Qtn0 / norm(SVt.Qtn0);
SVt.A0 = qtn2dc(SVt.Qtn0)';
SVt.v0 = SV.v0 + k3_v0;
SVt.w0 = SV.w0 + k3_w0;
SVt.q  = SV.q  + k3_q;
SVt.qd = SV.qd + k3_qd;

% 4th Step
SVt.dQtn0 = w2dQtn(SVt.w0, SVt.Qtn0);
SVt = f_dyn(LP, SVt);

k4_R0 = d_time * SVt.v0;
k4_Qtn0 = d_time * SVt.dQtn0;
k4_v0 = d_time * SVt.vd0;
k4_w0 = d_time * SVt.wd0;
k4_q  = d_time * SVt.qd;
k4_qd = d_time * SVt.qdd;

% Solution
SV.R0 = SV.R0 + ( k1_R0 + 2*k2_R0 + 2*k3_R0 + k4_R0 )/6;
SV.Qtn0 = SV.Qtn0 + ( k1_Qtn0 + 2*k2_Qtn0 + 2*k3_Qtn0 + k4_Qtn0 )/6;
SV.Qtn0 = SV.Qtn0 / norm(SV.Qtn0);
SV.A0 = qtn2dc(SV.Qtn0)';
SV.v0 = SV.v0 + ( k1_v0 + 2*k2_v0 + 2*k3_v0 + k4_v0 )/6;
SV.w0 = SV.w0 + ( k1_w0 + 2*k2_w0 + 2*k3_w0 + k4_w0 )/6;
SV.q  = SV.q  + ( k1_q  + 2*k2_q  + 2*k3_q  + k4_q  )/6;
SV.qd = SV.qd + ( k1_qd + 2*k2_qd + 2*k3_qd + k4_qd )/6;


%%% EOF
