%***** 2D diffusion model of heat transport *******************************

%*****  Initialise Model Setup

% create x-coordinate vectors
xc = h/2:h:W-h/2;      % x-coordinate vector for cell centre positions [m]
zc = h/2:h:D-h/2;      % z-coordinate vector for cell centre positions [m]
xf = 0:h:W;            % x-coordinate vectore for cell face positions [m]
zf = 0:h:D;            % z-coordinate vectore for cell face positions [m]
[Xc,Zc] = meshgrid(xc,zc);  % create 2D coordinate arrays

% set up index array for boundary conditions
ix = [ 1,1:Nx,Nx ];  % closed/insulating sides
iz = [ 1,1:Nz,Nz+1 ];  % closed/insulating top, flux grad at bottom

% set initial condition for temperature at cell centres
T   = T0 + dTdz(2).*Zc;  % initialise T array on linear gradient

%*****  Solve Model Equations

t = 0;  % initial time [s]
tau = 0;  % initial time step count
dt = 1;

while t <= tend

    % increment time and step count
    t = t+dt;
    tau = tau+1;

    % 4th-order Runge-Kutta time integration scheme
            
    dTdt1 = diffusion(T,k0,h,ix,iz);
    dTdt2 = diffusion(T+dTdt1/2*dt,k0,h,ix,iz);
    dTdt3 = diffusion(T+dTdt2/2*dt,k0,h,ix,iz);
    dTdt4 = diffusion(T+dTdt3  *dt,k0,h,ix,iz);

    T = T + (dTdt1 + 2*dTdt2 + 2*dTdt3 + dTdt4)/6 * dt + (Hr ./ rho ./ Cp);

 
    % plot model progress every 'nop' time steps
    if ~mod(tau,nop)
        makefig(xc,zc,T);
    end

   
end

%*****  Utility Functions  ************************************************

% Function to calculate diffusion rate
function [dTdt] = diffusion(f, dTdzt, dTdzb, k0, h, ix, iz)

% calculate heat flux by diffusion
qx = - diff(k0(:, ix)) .* diff(f(:, ix))/h;
qz = - diff(k0(iz, :) .* diff(f(iz, :))/h;

% calculate heat flux by diffusion
qz(1, :) = k0(1, :) .* dTdzt
qz(end, :) = k0(end, :) .* dTdzb

% calculate flux balance for rate of change
dTdt = - (diff(qx)/h + diff(qz)/h);

end

% Function to make output figure
function makefig(x,z,T)

clf; 

% plot temperature in subplot 1
imagesc(x,z,T); axis equal tight; colorbar; hold on
contour(x,z,T,[100,150,200],'k');
ylabel('z [m]','FontSize',15)
title('Temperature [C]','FontSize',17)

drawnow;

end