% File: euler_matrix.m
% Author: Jeremy Engels
% Date: 15 November 2019
% Description: returns the 3x3 Euler rotation matrix for transformation of
%   angle phi about axis n (n = 1,2,3)

function R = euler_matrix(n, phi)
    if n == 1
        R = [1 0 0; 0 cos(phi) sin(phi); 0 -sin(phi) cos(phi)];
    elseif n == 2
        R = [cos(phi) 0 -sin(phi); 0 1 0; sin(phi) 0 cos(phi)];
    elseif n == 3
        R = [cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1];
    else
        error('n must be between 1 and 3');
    end
end