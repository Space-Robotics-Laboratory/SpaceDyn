%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    version 1 // Oct 4, 1999, Last modification by K.Yoshida
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% 	RPY2DC	Roll/pitch/yaw angles to direction cosine matrix
%
%		RPY2DC(R,P,Y) returns a 3x3 direction cosine matrix
%		for the specified roll/pitch/yaw angles.  
%		These correspond to rotations about the X, Y, Z axes respectively.
%
%		See also DC2RPY, EUL2DC.
%

function C = rpy2dc( roll , pitch , yaw )


if length(roll) == 3
   C = cz( roll(3) ) * cy( roll(2) ) * cx( roll(1) );

else
    C = cz( yaw ) * cy( pitch ) * cx( roll );

   
end


%%%EOF
