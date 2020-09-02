%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   version 1 // Oct 4, 1999, Last modification by K.Yoshida
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% R_NE		Inverse Dynamic computation by the Recursive Newton-Euler method
%
%		R_NE returns a generalized force, which consists of
%		the reaction forces FF0, TT0 on the link 0, and torque tau
%		of each joint.
%
%
%   1997.11.25 T.Hiraoka
%   1998.2.13  K.Fujishima
%   2000.9.6   H.Hamano
%   2002.2.28  A.Irobe modified to version up
%

function Force = r_ne( LP,SV )

global Ez Gravity


% Calculation of coordinate transfromation matrices
A_I_0 = SV.A0;

% Calculation of velocity vectors vv,ww
% NOTE:	vv,ww are in the Inertial frame

SV = calc_vel( LP,SV );

% Calculation of acceletion vectors vd,wd
% NOTE:	vd,wd are in the Inertial frame

SV = calc_acc( LP,SV );


% Calculation of inertia force & torque of link 0
% NOTE:	FF,TT(FF0,TT0) are in the Inertial frame.

FF0 = LP.m0 * (SV.vd0-Gravity);
TT0 = (A_I_0*LP.inertia0*A_I_0')*SV.wd0 ...
   + cross( SV.w0 , ((A_I_0*LP.inertia0*A_I_0')*SV.w0) );

% Calculation of inertia force & torque of link i
% from link 1 to n
% Single or multi body ?

if LP.num_q == 0
   % If a Single body,
   FF = [];
   TT = [];
   
   % If a Multi body system,
else
   FF = zeros(3,LP.num_q);
   TT = zeros(3,LP.num_q);
   for i = 1 : LP.num_q
      A_I_i  = SV.AA(:,i*3-2:i*3);
      in_i = LP.inertia(:,i*3-2:i*3);
      FF(:,i) = LP.m(i) * (SV.vd(:,i)-Gravity) ;
      TT(:,i) = (A_I_i*in_i*A_I_i')*SV.wd(:,i) ...
         + cross( SV.ww(:,i) , ((A_I_i*in_i*A_I_i')*SV.ww(:,i)) );
   end
   
end


% Equilibrium of forces & torques on each link
% On the i-th link

Fj = zeros(3,LP.num_q);
Tj = zeros(3,LP.num_q);

if (LP.num_q ~= 0)
   
   % Multi body system
   % from link n to 1
   
   for i = LP.num_q : -1 : 1
      
      F = zeros(3,1);
      T = zeros(3,1);
      
      for j=i+1:LP.num_q
         
         F =  F + LP.SS(i,j)*Fj(:,j);
         
      end
      
      Fj(:,i) = FF(:,i) + F + LP.SE(i)*SV.Fe(:,i);
      
      for j=i+1:LP.num_q
         
         A_I_i = SV.AA(:,i*3-2:i*3);
         T =  T ...
            + LP.SS(i,j)*( cross(A_I_i*(LP.cc(:,i,j)-LP.cc(:,i,i)+(LP.J_type(i)=='P')*Ez*SV.q(i)),Fj(:,j) ) + Tj(:,j) );
         
      end
      
      if LP.J_type(i) == 'R'
         % Rotational joint
         Tj(:,i) = TT(:,i) + T ...
            - cross( A_I_i*LP.cc(:,i,i),FF(:,i) ) ;
         
      else
         % Prismatic joint
         Tj(:,i) = TT(:,i) + T ...
            + cross( A_I_i*(Ez*SV.q(i))-A_I_i*LP.cc(:,i,i) , FF(:,i) );
         
      end
      
      Tj(:,i) = Tj(:,i) ...
         + LP.SE(i)*( cross(A_I_i*(LP.ce(:,i)-LP.cc(:,i,i)+(LP.J_type(i)=='P')*Ez*SV.q(i)),SV.Fe(:,i) ) + SV.Te(:,i) );
      
   end
   
   
   % Equilibrium on the link 0
   
   F = zeros(3,1);
   T = zeros(3,1);
   
   for i=1:LP.num_q
      
      if (LP.S0(i) ~= 0)
         F =  F + LP.S0(i)*Fj(:,i);
         
      end
      
   end
   
   FF0 = FF0 + F;
   
   for i=1:LP.num_q
      
      if (LP.S0(i) ~= 0)
         T = T + LP.S0(i)*( cross( (A_I_0*LP.c0(:,i)) , Fj(:,i) ) + Tj(:,i) );
         
      end
      
   end
   
   TT0 = TT0 + T;
   
end


% Calculation of torques of each joint

% Single body

if LP.num_q == 0
   SV.tau = zeros(0);
   
else
   % Multi body system
   
   for i = 1 : LP.num_q
      
      A_I_i = SV.AA(:,i*3-2:i*3);
      
      if LP.J_type(i) == 'R'
         % Rotational joint
         SV.tau(i,1) = Tj(:,i)'*(A_I_i*Ez);
         
      else
         % Prismatic joint
         SV.tau(i,1) = Fj(:,i)'*(A_I_i*Ez);
         
      end
      
   end
   
end


% Compose a generalized force

% Single or multi body ?

if LP.num_q == 0
   % Single body,
   Force = [ FF0' TT0' ]';
   
else
   % Multi body system,
   Force = [ FF0' TT0' SV.tau' ]';
   
end


%%%EOF
