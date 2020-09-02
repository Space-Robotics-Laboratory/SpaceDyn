%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   version 1 // Oct 4, 1999, Last modification by K.Yoshida
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% 	EUL2DC	3-1-3 Euler angles to transformation matrix
%
%		EUL2DC(Z1,X2,Z3) returns a Direction Cosine matrix 
%		for the specified Euler angles.  
%		These correspond to rotations about the
%		Z, X, Z axes respectively.
%
%		See also DC2EUL, RPY2DC.
%

function C = eul2dc( phi , theta , psi )


if length( phi ) == 3
   C = cz( phi(3) ) * cx( phi(2) ) * cz( phi(1) );
   
else
   C = cz( psi ) * cx( theta ) * cz( phi );
   
end

%%% EOF