%***** RUN 2D MODEL FROM IMAGE ***********************************

% clear workspace
clear all; close all; %clc;

% load model setup from image, interpolate to target grid size
W       = 16e3;     % domain width (must correspond to width of image) [m]
Nx      = 200;      % target no. of columns in x-direction
h       = W/Nx;     % grid spacing based on image width and target grid size
n_units = 9;        % number of rock units contained in image

% units = value of each pixel (colour)
% D = original depth
% Nz = target no. of rows in z-direction
[units,D,Nz] = ModelFromImage('section.tiff',n_units,W,Nx);

% material properties for each rock unit (update based on your calibration)

matprop = [
% unit  conductivity  density  heat capacity  heat production
   1	    3.6788	    2697.6	    1000	    4.172
   2	    1	        2000	    1000	    1
   3	    3.2197	    2703.5	    1000	    5.575
   4	    1	        2000	    1000	    1
   5	    1	        2000	    1000	    1
   6	    0.924	    2083.1	    1000	    1
   7	    1	        2000	    1000	    1
   8	    0.919	    1905.9	    1000	    1
   9	    1e-6        1000	    1000	    0];  % air/water

% get coefficient fields based on spatial distribution of rock units from image

rho    = reshape(matprop(units,3),Nz,Nx); % density
Cp     = reshape(matprop(units,4),Nz,Nx); % specific heat capacity
kT     = reshape(matprop(units,2),Nz,Nx); % conductivity
Hr     = reshape(matprop(units,5),Nz,Nx); % heat rate

% calculate heat diffusivity [m2/s]
k0 = kT ./ (rho .* Cp);

% set model parameters


Ttop  = 0;            % surface temperature
Tbot  = 175;           % top/base T-gradient
nop   = 100;          % output figure produced every 'nop' steps
yr    = 3600*24*365;  % seconds per year [s]
tend  = 1e4*yr;       % stopping time [s]
%*****  RUN MODEL
run('./transect_2D.m');