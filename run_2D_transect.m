%***** RUN 2D MODEL FROM IMAGE ***********************************

% clear workspace
clear all; close all; %clc;

% load model setup from image, interpolate to target grid size

W       = 16e3;      % domain width (must correspond to width of image) [m]
Nx      = 200;       % target no. of columns in x-direction
h       = W/Nx;      % grid spacing based on image width and target grid size
n_units = 9;         % number of rock units contained in image
test = 'no';         % test simulation or not
verification = true; % close boundaries, no heat source, test energy conservation

% units = value of each pixel (colour)
% D = original depth
% Nz = target no. of rows in z-direction
[units,D,Nz] = ModelFromImage('section.tiff',n_units,W,Nx);

% material properties for each rock unit taken from [1] provided excel file,
% [2] Rybach and Cermak [1982], [3] Waples and Waples [2004], [4] British
% Geological Survey, [5] Encyclopedia Britannica, [6] Sammalj¨arvi et al.
% [2017], [7] Yu et al. [1993] and [8] Gao et al. [2021] 

matprop = [
% unit  conductivity  density  heat capacity  heat production  porosity

  1	    3.68	    2676	    1172 	      4.172e-6         1.03   % Granite phase 1 [1][3][7]
  2	    2.465	    2700	    979	          2.8e-6           1.2    % Basement gneiss, ChatGPT used to convert conductivity units [1][2][5][6]
  3	    3.2197	    2703.5	    1172	      5.575e-6         1.03   % Granite phase 2 [1][3][7]
  4	    0.77	    1942.3	    740	          0.75e-6          25.3   % Sand [1][2][3][4]
  5	    0.77	    2648	    740	          0.95e-6          32.5   % Gravel [1][2][3][4]
  6	    0.924	    2081.7	    860	          1.43e-6          0.55   % Clay, mudstone [1][2][3][4]
  7	    1.67	    1916	    910	          0.91e-6          17     % Silt [1][2][3][8]
  8	    0.919	    1909.78	    740	          0.75e-6          21.2   % Mud, silt, sand [1][2][3]
  9	    1e-6        1000	    1000	      0                100];  % air/water

% get coefficient fields based on spatial distribution of rock units from image

switch test % test case with constant coefficients

    case 'no'

        rho    = reshape(matprop(units,3),Nz,Nx);     % density
        Cp     = reshape(matprop(units,4),Nz,Nx);     % specific heat capacity
        kT     = reshape(matprop(units,2),Nz,Nx);     % conductivity
        Hr     = reshape(matprop(units,5),Nz,Nx);     % heat rate
        phi    = reshape(matprop(units,6),Nz,Nx)/100; % porosity percentage

        % weight variables for porosity, assuming pores are filled with
        % air. ChatGPT aided in this conversion.

        rho = (1 - phi) .* rho + phi .* 1000;
        Cp  = (1 - phi) .* Cp + phi .* 1000;   
        kT  = (1 - phi) .* kT + phi .* 1e-6;  
        Hr  = (1 - phi) .* Hr;          

    case 'yes'

        rho    = 2400*ones(Nz,Nx); % density
        Cp     = 1000*ones(Nz,Nx); % specific heat capacity
        kT     = ones(Nz,Nx);      % conductivity
        Hr     = ones(Nz,Nx);      % heat rate
       
end


% set model parameters
dTdz = [0, 35/1000];  % set boundary condition
T0  = 10;             % surface temperature [degree C]
Tair = 10;            % air temperature [degree C]
nop   = 10;          % output figure produced every 'nop' steps
yr    = 3600*24*365;  % seconds per year [s]
tend  = 1e6*yr;       % stopping time [s]
CFL   = 1/5;          % Time step limiter


%*****  RUN MODEL
run('./transect_2D.m');

