%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Sept. 16, 1998, Last modification by K.Yoshida
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% I_KINE	Inverse Kinematics
%		Calculate joint angle solution (q_sol) for
%		given end point position/orientation.
%		It takes an iterative approach from a specified
%		initial posture (q_init). 
%               
%   (C)Space Robotics Lab
%   1998.1.13 Akihide Kurosu
%   2001.9.17 H.Hamano
%   2002.2.27 H.Nakanishi modified to version up
%

function q_sol = i_kine( LP, SV, POS_e, ORI_e, q_init, num_e )

global Ez


% Check matrix size of POS_e and ORI_e
% if (size(POS_e)~=[3,1]) || (size(ORI_e)~=[3,3])
%    error('Dimensiones of Position and Orientation does not meet the requirement');
% end

% Set some values
loop_limit = 1000;
norm_limit = 1e-8;
nm    = 1;
count = 0;
gain  = 0.1*ones(6,1);
SV.q = q_init;

% Start convergent calculation
while ( nm > norm_limit )
   
   % Calculate the orientation/position of all link
   SV = calc_aa(LP, SV);
   SV = calc_pos(LP, SV);
   
   % Calculate the joint connection from base to endpoint
   joints = j_num( LP, num_e );
   
   % Present endpoint position/orientation (forward kinematics)
   [now_p, now_o] = f_kin_e( LP, SV, joints );
   
   % Calculate the error between present and goal position/orientation
   err = [];
   err_p = POS_e - now_p;
   err_o = tr2diff(ORI_e, now_o);
   err = [ err_p; err_o ];
   err = err.*gain;
   
   % Caluculate the Jacobian matrix
   Jacob = calc_je( LP, SV, joints );
   
   % Calculate joint angular velocity using Jacobian
   SV.qd = pinv( Jacob ) * err;
   
   % Next Joint angle
   SV.q = SV.q + SV.qd;
   
   % Calculate the norm of joint velocity
   nm = norm(SV.qd);
   
   % Check the loop condition
   count = count + 1;
   if ( count > loop_limit )
      error('Solution does not converge');
      break;
   end
   
end

q_sol = SV.q;

%%%EOF
