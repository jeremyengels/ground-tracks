% File: ground_track.m
% Author: Jeremy Engels
% Date: 10 November 2019
% Description: returns a ground track plot for the input orbital elements
%   for P periods. a should have units km and all the angles should have
%   units of degrees.

function ground_track(a,e,i,omega,Omega,theta0,P)

D2R = pi/180;
R2D = 180/pi;
H2S = 3600;

i = i*D2R;
omega = omega*D2R;
Omega = Omega*D2R;
theta0 = theta0*D2R;

% calculatio n setup
dt = 1;
tol = 0.0001;

% calculate various parameters
RE = 6378;
mu = 3.98e5;                    % km^3 / s^2 
h = sqrt(a*mu*(1-e^2));         % km^2 / s
T = 2*pi*sqrt(a^3/mu);          % sec
t = 0:dt:P*T;
N = length(t);

% incorporate J2 effects
J2 = 0.00108263;
Omega_dot = -1.5*sqrt(mu)*J2*RE^2*cos(i)/((1-e^2)^2*a^(7/2));
Omega = Omega + t*Omega_dot;
omega_dot = -1.5*sqrt(mu)*J2*RE^2*(5/2*sin2(i) - 2)/((1-e^2)^2*a^(7/2));
omega = omega + t*omega_dot;

% Kepler's Method to find theta

t0 = backward_kepler(theta0,T,e);

theta = forward_kepler(t+t0,T,e,tol);

% solve orbital equation
r = h^2/mu * 1./(1+e*cos(theta));
r_PF = zeros(3,N);

for j = 1:N
    r_PF(:,j) = r(j) * [cos(theta(j)); sin(theta(j)); 0];
end

% transform PF to GEC
r_GEC = zeros(3,N);
for j = 1:N
    r_GEC(:,j) = (euler_matrix(3,Omega(j)))'*(euler_matrix(1,i))'*...
        (euler_matrix(3,omega(j)))'*r_PF(:,j);
end

% calculate right ascension and declination angles
delta = (pi/2 - acos(r_GEC(3,:)'./r(:)))';

r_proj = zeros(3,N);
alpha = zeros(N,1);

for j = 1:N
    r_proj(:,j) = r_GEC(:,j) - [0; 0; r_GEC(3,j)];
    r_proj_x = r_proj(1,j)'/norm(r_proj(:,j));
    r_proj_y = r_proj(2,j)'/norm(r_proj(:,j));
    alpha(j) = acos(r_proj_x);
    if r_proj_y >= 0
        alpha(j) = acos(r_proj_x);
    else
        alpha(j) = pi - acos(r_proj_x);
    end
end

alpha = 0.5*unwrap(2*alpha);

% calculate ground track
omegaE = 15.04*D2R/H2S;
thetaE = omegaE*t';

alpha_GT = alpha - thetaE;
delta_GT = delta;

alpha_GT = alpha_GT - 2*pi*ceil(0.5*floor(alpha_GT/pi));

alpha_GT = alpha_GT*R2D;
delta_GT = delta_GT*R2D;

% plot without animations
figure(1)
geoplot(delta_GT,alpha_GT,'r.')
grid off
geolimits([0 30],[-180 180])
geobasemap satellite


end