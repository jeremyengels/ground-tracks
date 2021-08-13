% File: backward_kepler.m
% Author: Jeremy Engels
% Date: 13 November 2019
% Description: calculates time-since-perigee from true anomaly via backward
%   kepler's Method

function t = backward_kepler(theta,T,e)

E = 2*atan(((1-e)/(1+e))^(1/2)*tan(theta/2));
Me = E - e*sin(E);
t = Me*T/(2*pi);

end