%CROSS	Vector cross product
%
%	CROSS(V1, V2) returns the vector cross product V1 x V2
%
%	Copyright (C) Peter Corke 1990
%

function n = cross(u, v)
	n = zeros(3,1);
	n(1) = u(2)*v(3) - u(3)*v(2);
	n(2) = u(3)*v(1) - u(1)*v(3);
	n(3) = u(1)*v(2) - u(2)*v(1);
