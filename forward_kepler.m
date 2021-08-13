% File: forward_kepler.m
% Author: Jeremy Engels
% Date: 13 November 2019
% Description: calculates true anomaly from time since perigee via kepler's
%   Method, with Newton's Method

function theta = forward_kepler(t,T,e,tol)

% part 1
Me = 2*pi*t/T;

% part 2
N = length(Me);
actualE = zeros(N,1);

for p = 1:N
    E = zeros(1000,1);
    E(1) = 0;
    E(2) = Me(p);
    if Me(p) < tol
        E(2) = 1;
    end
    i = 2;

    while abs(E(i) - E(i-1)) > tol
        i = i + 1;
        E(i) = Me(p) + e*sin(E(i-1));
        if i == 10000
            error("Newton's Method does not converge within the given tolerance after 10000 iterations");
        end
    end

    actualE(p) = E(i);
end

% part 3
theta = 2*atan(tan(actualE/2)*((1-e)/(1+e))^(-1/2));

end

