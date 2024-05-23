%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    version 1 // Oct 4, 1999, Last modification by K.Yoshida
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% 	W2DQTN	Calclulate devision of quaternion,
%
%		dQTN = W2DQTN( W, QTN )
%
% Note: Here the "Quaternion" specifically means Euler Symmetric Parameters.
%       See p.414-415 in "Spacecraft Attitude Determination and Control,"
%       Ed. by J. R. Werts, Kluwer Academic Publishers.
%
%   2024.5.22 A.Uchida and M.Imai
%

function dQtn = w2dQtn(w, Qtn)

if ~all( size(w) == [3, 1])
  error("Cannot calculate dQtn: wrong w size")
end
if ~all( size(Qtn) == [4, 1])
  error("Cannot calculate dQtn: wrong Qtn size")
end
if (vecnorm(Qtn) - 1) >= .1
  error("Cannot calculate dQtn: Qtn is not normalized")
end

qx = Qtn(1, 1);
qy = Qtn(2, 1);
qz = Qtn(3, 1);
qw = Qtn(4, 1);

wx = w(1, 1);
wy = w(2, 1);
wz = w(3, 1);

dQtn = .5 * [ qw*wx + qz*wy - qy*wz;
             -qz*wx + qw*wy + qx*wz;
              qy*wx - qx*wy + qw*wz;
             -qx*wx - qy*wy - qz*wz];

end
%%% EOF