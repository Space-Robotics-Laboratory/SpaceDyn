%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	February 4, 1998, Last modification by K.Yoshida
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% CALC_JRE	Translational Jacobians w.r.t.
%           for the point specified by connection vector 'joints'
%
%   1998.1.12 A.Kurosu
%   2001.9.11 H.Hamano
%   2002.2.27 H.Hamano modified to version up
%

function JJ_re = calc_jre( LP, SV, joints )

global Ez


% Check number of joints
n = length(joints);


% If a Single body,
if LP.num_q == 0
   
   JJ_re = [];
   
% If a Multi body system,
else
   
   JJ_re = [];
   
   for i = 1 : 1 : n
      
%      A_I_i = SV.AA(:,joints(i)*3-2:joints(i)*3);
      
      % Rotational joint
      if LP.J_type(joints(i)) == 'R'
%         JJ_re = [ JJ_re A_I_i*Ez ];
         JJ_re = [ JJ_re SV.AA(:,joints(i)*3-2:joints(i)*3)*Ez ];
         
      % Prismatic joint
      else
         JJ_re = [ JJ_re [ 0 0 0 ]' ];
         
      end
      
   end
   
end


%%%EOF
