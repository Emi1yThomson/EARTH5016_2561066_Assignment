%***** 2D diffusion model of heat transport *******************************

%*****  Initialise Model Setup

% create x-coordinate vectors
xc = h/2:h:W-h/2;           % x-coordinate vector for cell centre positions [m]
zc = h/2:h:D-h/2;           % z-coordinate vector for cell centre positions [m]
xf = 0:h:W;                 % x-coordinate vectore for cell face positions [m]
zf = 0:h:D;                 % z-coordinate vectore for cell face positions [m]
[Xc,Zc] = meshgrid(xc,zc);  % create 2D coordinate arrays

% set up index array for boundary conditions
ix = [ 1,1:Nx,Nx ];  % closed/insulating sides
iz = [ 1,1:Nz,Nz ];  % closed/insulating top, flux grad at bottom

% set initial condition for temperature at cell centres
T   = T0 + dTdz(2).*Zc;  % initialise T array on linear gradient

% set up condition for air
air = units == 9;

%*****  Solve Model Equations

% calculate heat diffusivity [m2/s]
k0 = kT ./ (rho .* Cp);

t = 0;  % initial time [s]
tau = 0;  % initial time step count
dt = CFL * (h/2)^2/max(k0, [], 'all');


if verification 
    % calculate thermal capacitance of domain
    Cthsum = sum(rho(:) .* Cp(:) *h*h);

    Esum = []; % empty list for energy of domain
    t_vals = []; % empty list for time
end

while t <= tend

    % increment time and step count
    t = t+dt;
    tau = tau+1;

    
    if ~verification
        % reset air section
        T(air) = Tair;
    end

    % 4th-order Runge-Kutta time integration scheme
            
    dTdt1 = diffusion(T,dTdz,k0,h,ix,iz,verification);
    dTdt2 = diffusion(T+dTdt1/2*dt, dTdz, k0,h,ix,iz,verification);
    dTdt3 = diffusion(T+dTdt2/2*dt, dTdz,k0,h,ix,iz,verification);
    dTdt4 = diffusion(T+dTdt3  *dt, dTdz,k0,h,ix,iz,verification);

    if verification
        t_vals(end+1) = t;      
        T = T + (dTdt1 + 2*dTdt2 + 2*dTdt3 + dTdt4)/6 * dt;     
        Esum(end+1) = sum(T(:)) * Cthsum; % calculate total energy of domain
    end

    if ~verification
        T = T + (dTdt1 + 2*dTdt2 + 2*dTdt3 + dTdt4)/6 * dt + (Hr ./ rho ./ Cp);

        % plot model progress every 'nop' time steps
        if ~mod(tau,nop)
            makefig(xc,zc,T,t,yr);
        end
    end

end

if verification

    % plot energy over time
    clf; 
    plot(t_vals/yr, Esum, 'b-');
    xlim([min(t_vals/yr) max(t_vals/yr)]);
    xlabel('Time [yr]', 'FontSize', 15, 'FontName','Times New Roman');
    ylabel('Energy [J]', 'FontSize', 15, 'FontName','Times New Roman');
    title('Energy of Domain with Time', 'FontSize', 15,'FontName','Times New Roman');
end

%*****  Utility Functions  ************************************************

% Function to calculate diffusion rate
function [dTdt] = diffusion(f, dTdz, k0, h, ix, iz, verification)

% average k0 values to get cell face values
kx = k0(:, ix(1:end-1)) + k0(:, ix(2:end))/2;
kz = k0(iz(1:end-1), :) + k0(iz(2:end), :)/2;

% calculate heat flux by diffusion
qx = - kx .* diff(f(:, ix), 1, 2)/h;
qz = - kz .* diff(f(iz, :), 1, 1)/h;

if ~verification
    % set boundary conditions
    qz(end, :) = -kz(end, :) .* dTdz(2);
end

% calculate flux balance for rate of change
dTdt = -(diff(qx, 1, 2)/h + diff(qz, 1, 1)/h);

end


% Function to make output figure
function makefig(x,z,T,t,yr)

clf; 

% plot temperature
imagesc(x,z,T); axis equal tight; colorbar; hold on
ylabel('Depth [m]','FontSize',15, 'FontName','Times New Roman');
xlabel('Width [m]','FontSize',15, 'FontName','Times New Roman');
ylabel(colorbar, 'Temperature [\circC]','FontSize',15, 'FontName','Times New Roman')
title(['Temperature; Time = ',num2str(round(t/yr)),' yr'],'FontSize',17, 'FontName','Times New Roman');
[C,h] = contour(x,z,T, [45,90,135], 'k');
clabel(C,h,'Fontsize', 12, 'Color', ['r'], 'FontName','Times New Roman');

drawnow;

end