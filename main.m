% File: main.m
% Author: Jeremy Engels
% Date: 12 August 2021
% Description: Calls ground_track.m for a Molniya orbit

clc; clear; close all;

a = 26550;
e = 0.74;
i = 63.4349;
omega = 270;
Omega = 0;
theta0 = 175;

ground_track(a,e,i,omega,Omega,theta0,5);