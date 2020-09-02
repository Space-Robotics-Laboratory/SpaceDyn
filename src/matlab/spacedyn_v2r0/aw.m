%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    version 1.0 // Oct.4, 1999, Last modification by K.Yoshida
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	Aw( w0 ) returns a 3x3 transformation representing a 
%	rotation about the vector w0.
%
%   2002.2.27 H.Hamano modified to version up
%
%	See also CX, CY, CZ.

function E0 = aw( SV )

global d_time


if ( norm(SV.w0)==0 )
   E0 = eye(3);
   
else
   th = norm(SV.w0) * d_time;
   w = SV.w0 ./ norm(SV.w0);
   
   E0 =[ cos(th)+w(1)^2*(1-cos(th)) ...
         w(1)*w(2)*(1-cos(th))-w(3)*sin(th) ...
         w(3)*w(1)*(1-cos(th))+w(2)*sin(th);
         
         w(1)*w(2)*(1-cos(th))+w(3)*sin(th) ...
         cos(th)+w(2)^2*(1-cos(th)) ...
         w(3)*w(2)*(1-cos(th))-w(1)*sin(th);
         
         w(3)*w(1)*(1-cos(th))-w(2)*sin(th) ...
         w(3)*w(2)*(1-cos(th))+w(1)*sin(th) ...
         cos(th)+w(3)^2*(1-cos(th)) ];
   
end

%%%EO
