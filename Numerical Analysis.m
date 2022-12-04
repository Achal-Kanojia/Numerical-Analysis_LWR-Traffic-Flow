syms t x rho
clear all
clc

%% building a roadmap with length, demand, supply, Vmax parameters%%
roadmap.length(1) = 10;        % in km
roadmap.Demand(1) = @(rho) (80.*rho).*(rho<=40) + (3200).*(rho>40);
roadmap.Supply(1) = @(rho) (3200).*(rho<=40) + (20.*(100-rho)+2000).*(rho>40);
roadmap.Vmax(1) = 90;          % in km/hr

%% Plotting demand-supply curves %%
figure (1)
hold on
ezplot(roadmap.Demand,[0 200])
ezplot(roadmap.Supply,[0 200])
hold off
xlabel('Density (in veh/km)')
ylabel('Flow (in veh/hr)')
title('Function for Demand and Supply')
legend('Demand','Supply',"Location","best")

%% Initial and Boundary conditions %%
rho_0 = @(x) 20.*(x<=0.5) + 100.*(x>0.5) ;   % Initial density in veh/km
Upstream_Demand = @(t) 1800;                 % Upstream demand in veh/hr
Downstream_Supply = @(t) 3000;               % Downstream supply in veh/hr
Delta_x = 0.1;                               % Stepsize in space in km
T = 0.5;                                     % Total time in hr

%% Calling on Godunov function with given conditions %%
[rho,Delta_t]=Godunov(roadmap,rho_0,Upstream_Demand,Downstream_Supply,Delta_x,T) ;


%--------------------------------------------------------------------------------------%


%%%%%%%% Godunov Function %%%%%%%%%
function [rho,Delta_t]=Godunov(roadmap,rho_0,Upstream_Demand,Downstream_Supply,Delta_x,T)
% Inputs = length of the section, demand, supply function, max velocity,
% boundary condition, spatial stepsize and Time period
% Output = Density at each step, time stepsize satisfying CFL.

number_section = length(roadmap);       % total number of sections
V_max = 0;                              % Maximum Velocity
X = [];                                 % Matrix for defining avg location of each section
L_Section = zeros(1,number_section);    % Length of each section
L_Variable = 0;

for i = 1:number_section
    V_max = max( V_max, roadmap(i).Vmax) ;     % calculating max velocity
    % Computing total length of network by adding length in every section
    L_Section(i) = roadmap(i).length;
    Previous_Section = floor(L_Variable/Delta_x)+1 ;
    Current_Section = floor((L_Variable+L_Section(i))/Delta_x) ;
    X(2,Previous_Section:Current_Section) = i.*ones(1,Current_Section-Previous_Section+1) ;
    L_Variable = L_Variable + L_Section(i);
end

L = sum(L_Section);                         % Network's total length
X(1,:) = Delta_x/2:Delta_x:L-Delta_x/2;     % Network cells

%% Courant–Friedrichs–Lewy(CFL) Condition %%
security_factor = 1.5;         % SF>=1
Delta_t = Delta_x / (security_factor*V_max);   % Timestep satisfying CFL
disp(Delta_t)

T_size = length(0:Delta_t:T);       % Size of T
X_size = length(X(1,:));            % Size of X
rho = zeros(T_size,X_size);         % zero matrix for rho
rho(1,:) = rho_0(X(1,:));           % Assigning boundary values
temp = Delta_t / Delta_x ;

for i=1:T_size-1                    % Time loop
    
    for j = 1                       % Space loop
        index_j = X(2,j);           % index of cell
        index_j_next = X(2,j+1);    % index of next cell
        Supply = roadmap(index_j).Supply ;
        Demand = roadmap(index_j).Demand ;
        Supply_next = roadmap(index_j_next).Supply;
        inwards_flow = min( Upstream_Demand(i*Delta_t), Supply(rho(i,j)) );
        outwards_flow = min( Demand(rho(i,j)), Supply_next(rho(i,j+1)) );
        rho(i+1,j) = rho(i,j) + temp * (inwards_flow-outwards_flow);
        inwards_flow = outwards_flow;
    end
    
    for j= 2:X_size-1
        index_j = X(2,j);
        index_j_next = X(2,j+1);
        Demand = roadmap(index_j).Demand ;
        Supply_next = roadmap(index_j_next).Supply ;
        outwards_flow = min( Demand(rho(i,j)), Supply_next(rho(i,j+1)) );
        rho(i+1,j) = rho(i,j) + temp * (inwards_flow-outwards_flow);
        inwards_flow = outwards_flow;
    end
    
    for j = X_size
        index_j = X(2,j);
        Demand = roadmap(index_j).Demand ;
        outwards_flow = min( Demand(rho(i,j)), Downstream_Supply(i*Delta_t) );
        rho(i+1,j) = rho(i,j) + temp * (inwards_flow-outwards_flow);
    end
end

figure
[X,Y]=meshgrid(0:Delta_t:T,Delta_x/2:Delta_x:L-Delta_x/2);
surf(X,Y,rho','EdgeColor','none');
view(2)
axis tight
colormap(jet)
fig = colorbar ;
fig.Label.String = 'Density (veh/km)';
xlabel('Time (hr)','Fontsize',14)
ylabel('Space (km)','Fontsize',14)
end