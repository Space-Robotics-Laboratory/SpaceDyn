% init_SV.m
%
% initialization of state valuables.
%
% 2002.3.1   Hiroshi Hamano
%

function SV = init_SV( LP )

SV.R0 = [ 0 0 0 ]';
SV.Q0 = [ 0 0 0 ]';
SV.A0 = rpy2dc(SV.Q0)';
SV.AA = zeros(3,3*LP.num_q);
SV.RR = zeros(3,LP.num_q);

SV.v0 = [ 0 0 0 ]';
SV.w0 = [ 0 0 0 ]';
SV.vd0 = [ 0 0 0 ]';
SV.wd0 = [ 0 0 0 ]';

SV.q = zeros(LP.num_q,1);
SV.qd = zeros(LP.num_q,1);
SV.qdd = zeros(LP.num_q,1);

SV.vv = zeros(3,LP.num_q);
SV.ww = zeros(3,LP.num_q);
SV.vd = zeros(3,LP.num_q);
SV.wd = zeros(3,LP.num_q);

SV.F0 = [ 0 0 0 ]';
SV.T0 = [ 0 0 0 ]';
SV.Fe = zeros(3,LP.num_q);
SV.Te = zeros(3,LP.num_q);
SV.tau = zeros(LP.num_q,1);


%%% EOF
