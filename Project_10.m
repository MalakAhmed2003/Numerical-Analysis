clc;
clear all; 

% Problem:
          
   _______B______
% |              |
% |              |
% |              |
% | A            |C
% |              |
% |              |_____D_______
% |                            |
% |                            | E
% |___________F________________|

% You have an L-shaped thin plate with the following material properties and dimensions:
% - Thermal conductivity (k) = 4000 W/m·K
% - Specific heat capacity (c) = 700 J/kg·K
% - Density (ρ) = 4717 kg/m^3
% The dimensions of the plate are as follows:
% - A = 1m, F = 1m, B = 0.2m, E = 0.6m, C = 0.4m, D = 0.8m

% Boundary conditions:
% - The top edge (B) is maintained at a constant temperature of 150°C.
% - The right edge (E) is maintained at a constant temperature of 50°C.
% - The other edges (A, F, C, D) are insulated, meaning no heat can flow across these edges.
% - The initial temperature of all elements in the plate is 0°C.
% - The upper and lower surfaces are also insulated, so the solution is 2-dimensional only.

% Task:
% (a) Solve for the steady-state temperature distribution of the plate:
% - Solve for the steady-state temperature map, where the temperature does not change with time.
% - Plot the steady-state temperature distribution using a color plot or image.

% (b) Solve for the time-varying temperature solution for the plate:
% - Start at t = 0, with all nodes (except the fixed-temperature nodes) initialized at 0°C.
% - Run the simulation until the temperature converges to the steady state.
% - Plot the temperature distribution at four time instances: t = 25s, t = 150s, t = 800s, and at the final time when the solution has converged.


%% Part a
% We can firstly solve the problem as if it is a square and then when
% updating the values of the temperature, we skip he empty region
% The plate width 
Lx= 1; 
% The plate length
Ly= 1; 
% The number of nodes in x direction (x_resolution)
Nx=40; 
% The number of nodes in y direction (y_resolution)
Ny=40; 
d_x = Lx/Nx;
d_y = Ly/Nx;
% Initial temperature in all nodes
T_initial= 0; 
% The constant temperature on the top side (at y=Ny+1)
T_B = 150 ; 
% The constant temperature on the right side (at x=Nx+1)
T_E = 50 ;
% Tolerence for steady state section (correct up to 5 significant figures)
tolerence_ss=0.5*(10^(2-5)); 
% The length of the B boundary
B = .2;
% The length of the B boundary
E = .6;
start_point = 2;
% Calculating the C_boundary position
extra_bound_x_pos = ceil(start_point + B/(d_x));
% Calculating the D_boundary position
extra_bound_y_pos = ceil(start_point + E/(d_y));

% The iteration counter
k=1;
% To make sure that the updating loop will not run infinitely
max_iter = 1000;
% Initializing the error for the first iteration to run
err_SS_max(k)=1; 
% Max T to set axes limits in plotting
T_max=max([T_B T_E T_initial]); 
% Min T to set axes limits in plotting
T_min=min([T_B T_E T_initial]); 

% Initializing the temperature matrix based on the boundary conditions

CurrTemp_ss=zeros(Nx+start_point,Ny+start_point);  NextTemp=zeros(Nx+start_point,Ny+start_point);
CurrTemp_ss(:,1)=T_initial;                        NextTemp(:,1)=T_initial;  
% B_side initilization
CurrTemp_ss(:,Ny+start_point-1)=T_B;               NextTemp(:,Ny+start_point-1)=T_B;
CurrTemp_ss(:,Ny+start_point)=T_B;                 NextTemp(:,Ny+start_point)=T_B; 
% E_side initilization
CurrTemp_ss(Nx+start_point-1,:)=T_E;               NextTemp(Nx+start_point-1,:)=T_E;
CurrTemp_ss(Nx+start_point,:)=T_E;                 NextTemp(Nx+start_point,:)=T_E; 


while err_SS_max(k)>=tolerence_ss 
% Looping over the points in the grid(Plate)
for i=start_point:Nx 
for j=start_point:Ny
     if i > extra_bound_x_pos && j> extra_bound_y_pos
     else
        % Bottom edge is insulated (F)
        if j==start_point 
            NextTemp(i,j)= (CurrTemp_ss(i+1, j)+CurrTemp_ss(i-1, j)+ 2*CurrTemp_ss(i, j+1))/4;
        % Left edge is insulated (A)
        elseif i==start_point 
            NextTemp(i,j)=(2*CurrTemp_ss(i+1, j)+CurrTemp_ss(i, j+1)+CurrTemp_ss(i, j-1))/4;
        % C edge is insulated 
        elseif i == extra_bound_x_pos && j >= extra_bound_y_pos 
            NextTemp(i,j)=(2*CurrTemp_ss(i-1, j)+CurrTemp_ss(i, j+1)+CurrTemp_ss(i, j-1))/4;
        % D edge is insulated
        elseif j == extra_bound_y_pos && i >= extra_bound_x_pos 
            NextTemp(i,j)=(CurrTemp_ss(i+1, j)+CurrTemp_ss(i-1, j)+2*CurrTemp_ss(i, j-1))/4;
        else % Inner points normal update
            NextTemp(i, j) = (CurrTemp_ss(i+1, j)+CurrTemp_ss(i-1, j)+CurrTemp_ss(i, j-1)+CurrTemp_ss(i, j+1))/4;
        end
     end
end
end
% Updating the iteration k
k=k+1; 
% Calculating error
err_SS_max(k)=abs(max(max((NextTemp-CurrTemp_ss)./ NextTemp))); 
% Updating the current temperature matrix
CurrTemp_ss=NextTemp; 
end

% Plotting 
% Generate the x and y coordinates values (From points numbers --> coordinates)
x = (0:Nx+1) * d_x;  
y = (0:Ny+1) * d_y;  

% Creating a mesh grid for plotting
[X, Y] = meshgrid(x, y);  % X and Y are now 2D grids

% Set the points in the rectangle defined by extra_bound_x_pos and
% extra_bound_y_pos to NaN (the part of the plate that is empty, no metal)
% since we do not need to plot it
CurrTemp_plot = CurrTemp_ss(2:Nx+2, 2:Ny+2);
CurrTemp_plot(extra_bound_x_pos+1:Nx+2, extra_bound_y_pos+1:Ny+2) = NaN;
% Plotting the temperature distribution
figure(1);
hold on
title(sprintf('Temperature at steady state'));
surf(X, Y, CurrTemp_plot)  
cb = colorbar;
caxis([T_min T_max]); 
% View from the top
view(90, -90);  
xlim([0 Lx]); xlabel('Length');
ylim([0 Ly]); ylabel('Width');
zlim([T_min T_max]); zlabel('Temperature for the steady state solution');

%% Part b:
% Material name
name=('Graphene'); 
% Thermal conductivity value (W/m.K)
conductivity = 4000; 
% Specific heat value (J/kg K)
specific_heat = 700;
% Density value (kg/m^3)
density = 4717;
% The plate width 
Lx= 1; 
% The plate length
Ly= 1; 
% The number of nodes in x direction (x_resolution)
Nx=40; 
% The number of nodes in y direction (y_resolution)
Ny=40; 
d_x = Lx/Nx;
dx = d_x;
d_y = Ly/Nx;
dy = d_y;
% Initial temperature in all nodes
T_initial= 0; 
% The constant temperature on the top side (at y=Ny+1)
T_B = 150 ; 
% The constant temperature on the right side (at x=Nx+1)
T_E = 50 ;
% The length of the B boundary
B = .2;
% The length of the B boundary
E = .6;
start_point = 2;
% Calculating the C_boundary position
extra_bound_x_pos = ceil(start_point + B/(d_x));
% Calculating the D_boundary position
extra_bound_y_pos = ceil(start_point + E/(d_y));

% tolerence for numerical simulation (0.1 deg Celesius difference between
% the time domain final solution and the steady-state solution)
tolerence = .1; 

% k, thermal diffusivity, is the thermal conductivity divided by density and specific heat capacity
k = conductivity/(specific_heat*density); 
% Max T to set axes limits in plotting
T_max=max([T_B T_E T_initial]); 
% Min T to set axes limits in plotting
T_min=min([T_B T_E T_initial]); 
% Setting dt via the stability criterion
dt=0.125*(((Lx/Nx)^2+(Ly/Ny)^2)/k)







% Initializing temperature matrices
% set max iterations 75,000 due to memory limitations, max time steps
T=zeros(Nx+2,Ny+2,75000); 
T(:,1,:)=T_initial;
T(:,Ny+1,:)=T_B;
T(:,Ny+2,:)=T_B; 
T(Nx+1,:,:)=T_E;
T(Nx+2,:,:)=T_E; 
T(:,:,1)=T_initial;

iter=1;
% Initialising error for the first iteration
err_E_max(iter)=100; 


while err_E_max(iter)>=tolerence && iter <11
for i=2:Nx
for j=2:Ny

if i > extra_bound_x_pos && j> extra_bound_y_pos 
     else
        if j==start_point % Bottom edge is insulated (F)
             T(i,j,iter+1)=T(i,j,iter)+(dt*k*(T(i-1,j,iter)-2*T(i,j,iter)+T(i+1,j,iter))/dx^2)+(dt*k*(-2*T(i,j,iter)+2*T(i,j+1,iter))/dy^2);
        elseif i==start_point % Left edge is insulated (A)
            T(i,j,iter+1)=T(i,j,iter)+dt*k*((-2*T(i,j,iter)+2*T(i+1,j,iter))/dx^2)+(dt*k*(T(i,j-1,iter)-2*T(i,j,iter)+T(i,j+1,iter))/dy^2);
        elseif i == extra_bound_x_pos && j >= extra_bound_y_pos % C edge is insulated 
            T(i,j,iter+1)=T(i,j,iter)+dt*k*((-2*T(i,j,iter)+2*T(i-1,j,iter))/dx^2)+(dt*k*(T(i,j-1,iter)-2*T(i,j,iter)+T(i,j+1,iter))/dy^2);
        elseif j == extra_bound_y_pos && i >= extra_bound_x_pos % D edge is insulated
             T(i,j,iter+1)=T(i,j,iter)+(dt*k*(T(i-1,j,iter)-2*T(i,j,iter)+T(i+1,j,iter))/dx^2)+(dt*k*((-2*T(i,j,iter)+2*T(i,j-1,iter))/dy^2));
        else % Inner points normal update
            T(i,j,iter+1)=T(i,j,iter)+(dt*k*(T(i-1,j,iter)-2*T(i,j,iter)+T(i+1,j,iter))/dx^2)+(dt*k*(T(i,j-1,iter)-2*T(i,j,iter)+T(i,j+1,iter))/dy^2);
        end

     end
end
end
iter=iter+1;
% Calculating the error
err_E_max(iter) = max(max(abs((T(:,:,iter) - CurrTemp_ss))));



% Testing solution convergence
if round(err_E_max(iter),5)==round(err_E_max(iter-1),5) && err_E_max(iter)~= 0
errordlg('The solution is not converging, Please choose a larger tolerence','Tolerence Error');
close(message)
return
end

end

% Plotting
x=zeros(1,Nx+2);
y=zeros(1,Ny+2); 
for i = 1:Nx+2
x(i) =(i-1)*dx;
end
for i = 1:Ny+2
y(i) =(i-1)*dy;
end

T=T(:,:,1:iter); % delete the unused preallocated zero layers
SStime=iter*dt; % steady state time
SSsteps=iter; %total number of steady state steps

% Animated plot
f2=figure(2);
ax2=axes(f2);
CurrTemp_plot = T(2:Nx+2, 2:Ny+2, :);
CurrTemp_plot(extra_bound_x_pos+1:Nx+2, extra_bound_y_pos+1:Ny+2, :) = NaN;
for j=floor(linspace(1,SSsteps,250))
surf(ax2,x,y,CurrTemp_plot(:,:,j))
title(ax2,sprintf('Temperature at time : %i seconds ',round(j*dt)))
cb=colorbar(ax2);
caxis(ax2,[T_min T_max]);
view(ax2,90,-90);
xlim(ax2,[0 Lx+dx]); xlabel('Length');
ylim(ax2,[0 Ly+dy]); ylabel('Width');
zlim(ax2,[T_min T_max]); zlabel('Temperature for the time varying solution');
drawnow
end

%%
err_E_max(iter)
