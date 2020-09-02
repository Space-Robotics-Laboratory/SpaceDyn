%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INT_SR	Integration of the system motion by Simpson's rule
%
%		INT_SR returns velocity & position vectors 
%		of each link of the system in the inertia frame.
%		It uses the following simple Simpson's rule.
%
%       f(t+dt) = f(t) + ( dt/6 ) * [ g(t) + 4*g(t+dt/2) + g(t+dt)]
%       where df(t)/dt = g(t)
%
%       g(t+dt/2) and g(t+dt) are updated inside this function.
%
%		A constant time step is used.
%
%		[R0,A0,v0,w0,q,qd]=int_sr(R0,A0,v0,w0,vd0,wd0,q,qd,qdd)
%	
%		NOTE:v0,w0,vd0,wd0 are in the inertia frame.
%
%   2011.8. N.Uyama wrote 1st version
%


function  SV = int_sr( LP , SV )

global d_time

% set time step to half
d_time_original = d_time;
d_time = d_time / 2;

% 1st step

SVt = SV;

g1_v0  = SVt.v0;
g1_Cd0 = (tilde(SVt.w0))' * SVt.A0';
g1_vd0 = SVt.vd0;
g1_wd0 = SVt.wd0;
g1_qd  = SVt.qd;
g1_qdd = SVt.qdd;

SVt.R0 = SVt.R0 + g1_v0 * d_time;
SVt.A0 = SVt.A0 + g1_Cd0' * d_time;
SVt.v0 = SVt.v0 + g1_vd0 * d_time;
SVt.w0 = SVt.w0 + g1_wd0 * d_time;
SVt.q  = SVt.q  + g1_qd * d_time;
SVt.qd = SVt.qd + g1_qdd * d_time;

% 2nd step

SVt = f_dyn(LP,SVt);

g2_v0  = SVt.v0;
g2_Cd0 = (tilde(SVt.w0))' * SVt.A0';
g2_vd0 = SVt.vd0;
g2_wd0 = SVt.wd0;
g2_qd  = SVt.qd;
g2_qdd = SVt.qdd;

SVt.R0 = SVt.R0 + g2_v0 * d_time;
SVt.A0 = SVt.A0 + g2_Cd0' * d_time;
SVt.v0 = SVt.v0 + g2_vd0 * d_time;
SVt.w0 = SVt.w0 + g2_wd0 * d_time;
SVt.q  = SVt.q  + g2_qd * d_time;
SVt.qd = SVt.qd + g2_qdd * d_time;

% 3rd step

SVt = f_dyn(LP,SVt);

g3_v0  = SVt.v0;
g3_Cd0 = (tilde(SVt.w0))' * SVt.A0';
g3_vd0 = SVt.vd0;
g3_wd0 = SVt.wd0;
g3_qd  = SVt.qd;
g3_qdd = SVt.qdd;

% reset time step
d_time = d_time_original;

% Return integration results
SV.R0 = SV.R0 + d_time / 6 * (g1_v0  + 4*g2_v0  + g3_v0);
SV.A0 = SV.A0 + d_time / 6 * (g1_Cd0 + 4*g2_Cd0 + g3_Cd0)';
SV.v0 = SV.v0 + d_time / 6 * (g1_vd0 + 4*g2_vd0 + g3_vd0);
SV.w0 = SV.w0 + d_time / 6 * (g1_wd0 + 4*g2_wd0 + g3_wd0);
SV.q  = SV.q  + d_time / 6 * (g1_qd  + 4*g2_qd  + g3_qd);
SV.qd = SV.qd + d_time / 6 * (g1_qdd + 4*g2_qdd + g3_qdd);

%%%EOF
