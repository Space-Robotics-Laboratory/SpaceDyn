%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    version 1 // Oct 4, 1999, Last modification by K.Yoshida
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	DC2EUL	Direction Cosine matrix to 3-1-3 EULER angles
%
%	[A B C] = DC2EUL(C) returns a vector of 3-1-3 Euler angles corresponding
%	from the direction cosine matrix C.
%
%   1997 11/19 T.Hiraoka
%
%	See also  EUL2DC, DC2RPY
%

function euler = dc2eul(C)

	euler = zeros(3,1);

	if ( abs( C(1,3) ) < eps & abs( C(2,3) ) < eps )

		euler(3) = 0;
		euler(2) = atan2( C(2,3) , C(3,3) );
		euler(1) = atan2( C(1,2) , C(1,1) );

	else

		euler(3) = atan2( C(1,3) , C(2,3) );
		s3 = sin(euler(3));
		c3 = cos(euler(3));
		euler(2) = atan2( s3*C(1,3)+c3*C(2,3) , C(3,3) );
		euler(1) = atan2( C(3,1) , -C(3,2) );

	end


%%%EOF
