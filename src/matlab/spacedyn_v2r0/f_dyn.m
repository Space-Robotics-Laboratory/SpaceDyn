%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   version 1 // Oct 4, 1999, Last modification by K.Yoshida
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% F_DYN		Caluclation of Forward Dynamics
%
%		F_DYN returns accelation vectors, vd0,wd0 and qdd.
%
%   1997.11.25 T.hiraoka
%   2002.2.27  H.Nakanishi modified to version up
%

function SV = f_dyn( LP, SV )

global Ez Gravity


% Calculation of coordinate transfromation matrices
SV = calc_aa(LP,SV);

% Calculation of position vectors
SV = calc_pos(LP,SV);

% Calculation of inertia matrices, HH
HH = calc_hh(LP,SV);

% Calculation of velocty dependent term, Force0
% This is obtained by the RNE inverse dynamics computation with
% the accerelations and external forces zero.
vd0_tmp = SV.vd0;
wd0_tmp = SV.wd0;
qdd_tmp = SV.qdd;
Fe_tmp = SV.Fe;
Te_tmp = SV.Te;

SV.vd0 = zeros(3,1);
SV.wd0 = zeros(3,1);
SV.qdd = zeros(LP.num_q,1);
SV.Fe  = zeros(3,LP.num_q);
SV.Te  = zeros(3,LP.num_q);

Force0 = r_ne(LP,SV);

% Force = forces on the generalized coordinate.
% Force_ex = forces on the end points.
Force = zeros(6+ LP.num_q,1);
Force_ex = zeros(6+ LP.num_q,1);

% F0, T0 are forces on the centroid of the 0-th body.
Force(1:3) = SV.F0;
Force(4:6) = SV.T0;

% If Multi body system, tau is a joint torque.
if ( LP.num_q ~= 0 )
   Force(7:LP.num_q+6) = SV.tau;
end

% Calculate external forces %

% If single body system, no external forces.
if LP.num_q == 0
   % Note that the body 0 cannot have an endpoint.
   Fx   = zeros(3,1);
   Tx   = zeros(3,1);
   taux = [];
   
% Multi body system
else
   Fx    = zeros(3,1);
   Tx    = zeros(3,1);
   taux  = zeros(LP.num_q,1);
   F_ex  = zeros(6+LP.num_q,LP.num_q);%add  2011/3/7
   
   E_3 = eye(3,3);
   O_3 = zeros(3,3);
   num_e = 1;
   
   for i = 1 : LP.num_q
      
      if LP.SE(i)==1
         joints = j_num(LP,num_e);
         tmp = calc_je(LP,SV,joints);
         JJ_tx_i = tmp(1:3,:);
         JJ_rx_i = tmp(4:6,:);
         
         num_e = num_e + 1;
         
         A_I_i = SV.AA(:,i*3-2:i*3);
         Re0i = SV.RR(:,i) - SV.R0 + A_I_i*LP.ce(:,i);
         
         Me_i = [         E_3      O_3;
                  tilde(Re0i)      E_3;
                     JJ_tx_i'  JJ_rx_i' ];
         F_ex(:,i) = Me_i * [ Fe_tmp(:,i) ; Te_tmp(:,i) ];%
         
      end
      
   end
   
   for i = 1 : LP.num_q
      
      Fx   = Fx   + F_ex(1:3,i);
      Tx   = Tx   + F_ex(4:6,i);
      taux = taux + F_ex(7:6+LP.num_q,i);
      
   end
   
end

Force_ex(1:3) = Fx;
Force_ex(4:6) = Tx;
Force_ex(7:6+LP.num_q) = taux;

% Calculation of the acclelation
a_Force = Force - Force0 + Force_ex;

% Acc = HH\a_Force;
Acc = HH\a_Force;

vd0_tmp = Acc(1:3);
wd0_tmp = Acc(4:6);
qdd_tmp = Acc(7:6+LP.num_q);

if LP.num_q == 0
   SV.qdd=[];
end

SV.vd0 = vd0_tmp;
SV.wd0 = wd0_tmp;
SV.qdd = qdd_tmp;
SV.Fe = Fe_tmp;
SV.Te = Te_tmp;


%%%EOF
