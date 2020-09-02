%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    version 1 // Oct 4, 1999, Last modification by K.Yoshida
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% INT_EU2	Integration of the system motion by a simple Euler method
%
%		INT_EU returns velocity & position vectors 
%		of each link of the system in the inertia frame.
%		It uses a simple, most primitive  Euler integration.
%
%		(The Rodorigues fomula for infinitesimal rotation
%		is used to update from A0 to A0_n.
%		This seems practically better than the fomula
%		delta_C0 = tilde(w0')*C0.)
%
%		A constant time step is used.
%
%		[R0,A0,v0,w0,q,qd]=int_eu2(R0,A0,v0,w0,vd0,wd0,q,qd,qdd)
%	
%		NOTE:v0,w0,vd0,wd0 are in the inertia frame.
%
%   2002.2.27 A.Irobe modified to version up
%


function  SV = int_eu2( LP , SV )

global d_time


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integlation of q and qd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Single body
if LP.num_q == 0
   q_n = [];
   qd_n = [];
   
% Multi bodies
else
   qd_n = SV.qd + SV.qdd * d_time;
   q_n  = SV.q  + SV.qd  * d_time;
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integration of v0, R0, and A0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R0_n = SV.R0 + SV.v0  * d_time;
A0_n = aw( SV ) * SV.A0;
v0_n = SV.v0 + SV.vd0 * d_time;
w0_n = SV.w0 + SV.wd0 * d_time;

% outputs

SV.R0 = R0_n;
SV.A0 = A0_n;
SV.v0 = v0_n;
SV.w0 = w0_n;
SV.q  = q_n;
SV.qd = qd_n;

%%%EOF
