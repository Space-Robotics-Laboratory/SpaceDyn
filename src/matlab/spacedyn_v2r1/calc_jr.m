%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	February 4, 1999, Last modification by K.Yoshida
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CALC_JR	Rotational Jacobians w.r.t. link centroid
%
%   1997.6.16 T.hiraoka
%   2001.9.12 H.Hamano
%   2002.2.27 H.Hamano modified to version up
%

function JJ_r = calc_jr( LP, SV )

global Ez


% Calculation of translational jacobians
JJ_r = zeros(3,LP.num_q*LP.num_q);

% If a Single body,
if LP.num_q == 0
   
   JJ_r = [];

% If a Multi body system,
else
   
   for i = 1 : LP.num_q
      
      % Rotational joint
      if LP.J_type(i) == 'R'
         JJ_r(:,(i-1)*LP.num_q+i) = SV.AA(:,i*3-2:i*3)*Ez;
         
      % Prismatic joint
      else
         JJ_r(:,(i-1)*LP.num_q+i) = [ 0 0 0 ]' ;
         
      end
      
      j = LP.BB(i);
      
      if j ~= 0
         JJ_r(:,(i-1)*LP.num_q+1:(i-1)*LP.num_q+i-1) = JJ_r(:,(j-1)*LP.num_q+1:(j-1)*LP.num_q+i-1);
         
      end
      
   end
   
end


%%%EOF
