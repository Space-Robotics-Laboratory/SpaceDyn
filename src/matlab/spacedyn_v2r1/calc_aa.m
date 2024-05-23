%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   version 1 // Oct 4, 1999, Last modification by K.Yoshida
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% CALC_AA 	Calculate the coordinate transfrom matrices in the Robotics convention.
%
%   CALC_TRN( A0 , q ) returns the coordinate tranformation matrices AA.
%   AA is a collection of A_I_1, A_I_2, ... A_I_n.
%
%   2002.2.27 H.Hamano modified to version up
%

function SV = calc_aa( LP, SV )


% If a Single body,
if LP.num_q == 0
   
   SV.AA = [];
   
% If a Multi body system,
else
   
   % Calculation of coordinate transformation matrices
   A_I_0 = SV.A0;
   
   for i = 1 : LP.num_q
      
      % Check the link connection: Is the lower one of this link, 0 ?
      if LP.BB(i) == 0
         
         % Current (i-th) link connects to the 0-th link
         
         % Rotational joint
         if LP.J_type(i) == 'R'
            A_0_i = (rpy2dc( LP.Qi(1,i) , LP.Qi(2,i) , LP.Qi(3,i)+SV.q(i) )');
            
         % Prismatic joint
         else
            A_0_i = (rpy2dc( LP.Qi(:,i) )');
            
         end
         
         SV.AA(:,i*3-2:i*3) = A_I_0*A_0_i;
         
      else
         
         % Current (i-th) link doesn't connect to the 0-th link
         
         % Rotational joint
         if LP.J_type(i) == 'R'
            A_BB_i = (rpy2dc( LP.Qi(1,i), LP.Qi(2,i), LP.Qi(3,i)+SV.q(i) )');
            
         % Prismatic joint
         else
            A_BB_i = (rpy2dc( LP.Qi(:,i) )');
            
         end
         
         SV.AA(:,i*3-2:i*3) = SV.AA(:,LP.BB(i)*3-2:LP.BB(i)*3)*A_BB_i;
         
      end
      
   end
   
end


%%%EOF
