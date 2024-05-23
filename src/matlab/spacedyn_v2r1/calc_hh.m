%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	February 4, 1999, Last modification by K.Yoshida
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% CALC_HH	Calculate the Inertia Matrices H.
%
%   CALC_HH returns the inertia matrices HH (6+n)x(6+n).
%
%
%   (C)Space Robotics Lab, by Koichi Fujishima
%   2001.9.12 modified by H.Hamano
%   2002.2.27 H.Hamano modified to version up
%

function HH = calc_hh( LP, SV )

global Ez


% Calculation of partial translational & rotational jacobian
JJ_t = calc_jt( LP, SV );
JJ_r = calc_jr( LP, SV );


% Calculation of HH matrix
wE = LP.mass * eye(3,3);

JJ_tg = zeros(3,LP.num_q);
HH_w = zeros(3,3);
HH_wq = zeros(3,LP.num_q);
HH_q  = zeros(LP.num_q,LP.num_q);


% If a Single body,
if LP.num_q == 0
   HH_w = SV.A0*LP.inertia0*SV.A0';
   HH = [         wE  zeros(3,3);
          zeros(3,3)        HH_w ];
   
   
% Multi body system
else
   
   % Calculation of the position of gravity centroid, Rg
   Rm = zeros(3,1);
   
   for i = 1 : LP.num_q
      
      Rm = Rm + LP.m(i) * SV.RR(:,i);
      
   end
   
   Rm = Rm + LP.m0 * SV.R0;
   Rg = Rm / LP.mass;
   
   wr0g = (Rg-SV.R0) * LP.mass;

   
   for i = 1 : LP.num_q
      
      r0i = SV.RR(:,i) - SV.R0;
      A_I_i = SV.AA(:,i*3-2:i*3);
      JJ_tg = JJ_tg + LP.m(i)*JJ_t(:,(i-1)*LP.num_q+1:i*LP.num_q);
      
      HH_w  = HH_w + A_I_i*LP.inertia(:,i*3-2:i*3)*A_I_i' ...
         + LP.m(i)*(tilde(r0i))'*(tilde(r0i));
      HH_wq = HH_wq ...
         + (A_I_i*LP.inertia(:,i*3-2:i*3)*A_I_i')*JJ_r(:,(i-1)*LP.num_q+1:i*LP.num_q) ...
         + LP.m(i) * (tilde(r0i)) * JJ_t(:,(i-1)*LP.num_q+1:i*LP.num_q);
      HH_q  = HH_q ...
         + JJ_r(:,(i-1)*LP.num_q+1:i*LP.num_q)'* ...
         (A_I_i*LP.inertia(:,i*3-2:i*3)*A_I_i')* ...
         JJ_r(:,(i-1)*LP.num_q+1:i*LP.num_q) ...
         + LP.m(i) * JJ_t(:,(i-1)*LP.num_q+1:i*LP.num_q)'*JJ_t(:,(i-1)*LP.num_q+1:i*LP.num_q);
      
   end
   
   HH_w = HH_w + SV.A0*LP.inertia0*SV.A0';
   HH = [          wE  tilde(wr0g)'  JJ_tg;
          tilde(wr0g)          HH_w  HH_wq;
               JJ_tg'        HH_wq'   HH_q ];
   
end


%%% EOF
