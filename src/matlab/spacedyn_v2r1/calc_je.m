%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	February 4, 1998, Last modification by K.Yoshida
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% CALC_JE	Calculation of the Jacobian Matrix
%           for the endpoint given by connection vector 'joints.'
%
%   Space Robotics Lab (C)
%   1997.1.14 A.Kurosu
%   2002.9.11 H.Hamano
%   2002.2.27 H.Hamano modified to version up
%

function Jacobian = calc_je( LP, SV, joints )

global Ez


% number of links from base to endpoint.
n = length(joints);

% Calculation of Jacobian
JJ_te = zeros(3,LP.num_q);
JJ_re = zeros(3,LP.num_q);

JJ_te = calc_jte( LP, SV, joints );
JJ_re = calc_jre( LP, SV, joints );
JJ = [ JJ_te; JJ_re ];

% Compose the Jacobian using the corresponding joints.
Jacobian = zeros(6,LP.num_q);

for i = 1 : 1 : n
   Jacobian(1:6 , joints(i)) = JJ(1:6 , i);
end


%%%EOF
