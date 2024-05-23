%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    version 1 // Oct 4, 1999, Last modification by K.Yoshida
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% 	QTN2DC	Convert quaternion to direction cosine matrix,
%
%		C0 = QTN2DC( QTN )
%
% Note: Here the "Quaternion" specifically means Euler Symmetric Parameters.
%       See p.414-415 in "Spacecraft Attitude Determination and Control,"
%       Ed. by J. R. Werts, Kluwer Academic Publishers.
%
%   2020.9.29 W. Ribeiro
%   2024.5.22 A.Uchida and M.Imai
%

function C0 = qtn2dc( qtn )

% input check-out
if ~all(size( qtn ) ~= [4, 1])
    error('Cannot compute the rotation matrix\n');
end

qx = qtn(1, 1);
qy = qtn(2, 1);
qz = qtn(3, 1);
qw = qtn(4, 1);

C0 = [ qx^2 - qy^2 - qz^2 + qw^2       2 * (qx * qy + qz * qw)       2 * (qx * qz - qy * qw);
         2 * (qx * qy - qz * qw)    -qx^2 + qy^2 - qz^2 + qw^2       2 * (qy * qz + qx * qw);
         2 * (qx * qz + qy * qw)       2 * (qy * qz - qx * qw)    -qx^2 - qy^2 + qz^2 + qw^2];

end
%%% EOF