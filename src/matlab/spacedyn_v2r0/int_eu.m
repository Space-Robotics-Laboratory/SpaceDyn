%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    version 1 // Oct 4, 1999, Last modification by K.Yoshida
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% INT_EU	Integration of the system motion by a simple Euler method
%
%		INT_EU returns velocity & position vectors 
%		of each link of the system in the inertia frame.
%		It uses a simple, most primitive  Euler integration.
%
%		(The fomula of delta_CO = tilde(w0')*C0 is used.
%		But its integration is subject to error accumulation,
%		yielding non-normal, non-orthogonal C0.
%		int_eu2 is recommended for large attitude motion.) 
%
%		A constant time step is used.
%
%		[R0,A0,q,v0,w0,qd]=int_eu(R0,A0,q,v0,w0,qd,vd0,wd0,qdd)
%	
%		NOTE:v0,w0,vd0,wd0 are in the inertia frame.
%
%   2002.2.28 H.Hamano modified to version up
%


function SV = int_eu( LP, SV )

global d_time


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
v0_n = SV.v0 + SV.vd0 * d_time;
w0_n = SV.w0 + SV.wd0 * d_time;

% Note that C0 is the direction cosines of body 0
% C0 = ( A0 )^T

C0 = ( SV.A0 )';

% dC0 is the time derivative of C0
dC0 = (tilde(SV.w0))' * C0;
C0_n = C0 + dC0 * d_time;

% outputs
SV.R0 = R0_n;
SV.A0 = ( C0_n )';
SV.v0 = v0_n;
SV.w0 = w0_n;
SV.q  = q_n;
SV.qd = qd_n;


%%%EOF
