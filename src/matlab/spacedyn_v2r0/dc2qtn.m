%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    version 1 // Oct 4, 1999, Last modification by K.Yoshida
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% 	DC2QTN	Convert direction cosine matrix to quaternion,
%
%		QTN = DC2QTN( C0 )
%
%		ex. Q0 = dc2qtn( C0 );
%		ex. Q0 = dc2qtn( A0' );
%
% Note: Here the "Quaternion" specifically means Euler Symmetric Parameters.
%       See p.414-415 in "Spacecraft Attitude Determination and Control,"
%       Ed. by J. R. Werts, Kluwer Academic Publishers.
%
%   1998.6.11 K.FUJISHIMA
%   1998.9.16 Rewritten by K.Yoshida
%

function qtn = dc2qtn( C0 )

% input check-out

	if [3,3]~=size( C0 )
		error('Cannot compute the quaternion\n');
	end
	check = det(C0);
	if abs(check-1)>0.1,
		fprintf('Warning: mat2qtn: input is not a normal matrix.\n');
	end

% start
	
	qtn    = zeros(4,1);

	qtn(4) = sqrt( 1 + C0(1,1) + C0(2,2) + C0(3,3) ) / 2;

	if qtn(4) > 0.001

% In case qtn(4) is not small,

		qtn(1:3) = [ C0(2,3)-C0(3,2);
		             C0(3,1)-C0(1,3);
	        	     C0(1,2)-C0(2,1) ] / (4*qtn(4));
	else

% In case qtn(4) is small,

		qtn(1) = sqrt( 1 + C0(1,1) - C0(2,2) - C0(3,3) ) / 2;
		if qtn(1) > 0.001,
			qtn(2:4) = [ C0(1,2)+C0(2,1);
		            	     C0(3,1)+C0(1,3);
	        	    	     C0(3,2)+C0(2,3) ] / (4*qtn(1));
		else
			det(C0)
			error('Cannot compute the quaternion. Check if the direction cosine is a normal orthogonal matrix.\n');
		end

	end

% precision check-out
	
	check = qtn'*qtn;
	if abs(check-1)>0.1
		fprintf('Warning: mat2qtn: precision in attitude computation is not guaranteed.\n')
	end

%%% EOF





