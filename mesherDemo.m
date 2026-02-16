
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mesherDemo.m
% Demo on constructing O-grid for NACA 0012 Airfoil inside structural
% domain. Based on Karman-Trefftz transformation of airfoil. 
% Author: Yiren Shen 
% Date: July-12-2023
% Version: 1.0
% Date: July-19-2023
% Version: 1.1: produce tet mesh
% Version: 1.5: refactorized utility functions, export laminated 2d strain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% import utillity function
util = surf3dutil();


% cal NACA 0012 Nodes 
num_nodes_circular = 51;
[x_up, y_up, x_lo, y_lo] = util.drawNACA([0 0 1 2], num_nodes_circular);
% define z 

x_all = [flip(x_lo); x_up(2:end)];
y_all = [flip(y_lo); y_up(2:end)];

z_all = util.constructZ(x_all, y_all);



% plot figures 
PLOT = 0;
if PLOT 
    figure 
    plot(x_up,y_up, '.--', x_lo, y_lo, '.--')
    xlim([0,1])
    ylim([-0.2, 0.2])
    xlabel('x')
    ylabel('y')
end

% convert nodes into conformal coordinate

% THESE constants are for NACA 0012 Airfoil only.
% TODO: function to calculate csts from airfoil profile.



% calculate rho and tau based on airfoil
[rho_LE, tau_TE] = util.airfoilPropCal(z_all);
tau = tau_TE;
z1 = 1.0089304115; % put z1 at TE of airfoil
z2 = 0.5*rho_LE; % put z2 at 1/2 leading edge radii 
P = pi/(2*pi - tau);
zeta1 = 0.77043505; % defined based on TE location 
zeta2 = 1-real((1-zeta1)*((1-z1)/(1-z2))^(-P));



z_ellipse = util.drawEllipse(z_all, 0.25, 0 )';
% fit parameters of e_ellipse into karman trefftz shape 
[rho_LE_ellipse, tau_TE_ellipse] = util.airfoilPropCal(z_ellipse);
P_ellipse = pi/(2*pi - tau_TE_ellipse);

z1_ellipse = max(real(z_ellipse));
z2_ellipse = min(real(z_ellipse)) + rho_LE_ellipse/2;
eps = 0.01;
zeta1_ellipse = 0.77043505;
zeta2_ellipse = 1-real((1-zeta1_ellipse)*((z1_ellipse+eps-z1_ellipse)/(z1_ellipse+eps-z2_ellipse))^(-P_ellipse));

%%% TODO: The ellipse conformal mapping is buggy. Need to look into zeta1
%%% and zeta 2 calculation. 

% 
% tau = 0.2818725;        %TE angle, for NACA0012 ONLY! 
% P = pi/(2*pi - tau);
% z1 = 1.0089304115;
% z2 = 0.0079337;
% zeta1 = 0.77043505;
% zeta2 = 0.0079337; 
% zetac = 0.485916;




%% Put airfoil into conformal coordinate 
zeta_all = zeros(length(z_all),1);

for ii = 1:length(zeta_all)
    cur_z = z_all(ii);
    cur_zeta = util.forwardKarmanTrefftz(cur_z, z1, z2, zeta1, zeta2, P);
    zeta_all(ii) = cur_zeta;
end

% plot in the conformal plane 
PLOT = 0;


if PLOT 
    figure 
    plot(real(zeta_all),imag(zeta_all), '.--')
    axis equal
    xlabel('x')
    ylabel('y')
    title('NACA0012 in \zeta plane')
end



%% put ellipse into conformal coordinate: 

zeta_ellipse_all = zeros(length(z_ellipse),1);


for ii = 1:length(zeta_ellipse_all)
    cur_z = z_ellipse(ii);
    cur_zeta = util.forwardKarmanTrefftz(cur_z, z1_ellipse, z2_ellipse, zeta1_ellipse, zeta2_ellipse, P_ellipse);
    zeta_ellipse_all(ii) = cur_zeta;
end

% plot in the conformal plane 
PLOT = 0;


if PLOT 
    figure 
    plot(real(zeta_ellipse_all),imag(zeta_ellipse_all), '.--')
    axis equal
    xlabel('x')
    ylabel('y')
    title('central ellipse in \zeta plane')
end





%% put back into physical coordinate: 

%% NACA
z_maped = zeros(length(zeta_all),1);
for ii = 1:length(zeta_all)
    cur_zeta = zeta_all(ii);
    cur_z = util.backwardKarmanTrefftxz(cur_zeta, z1, z2, zeta1, zeta2, P);
    z_maped(ii) = cur_z;
end


% plot in the conformal plane 
PLOT = 0;

if PLOT 
    figure 
    plot(real(z_maped),imag(z_maped), '.--')
    axis equal
    xlabel('x')
    ylabel('y')
    title('NACA0012 in z plane')
end


%% Ellipse 

z_ellipse_mapped = zeros(length(zeta_ellipse_all),1);
for ii = 1:length(z_ellipse_mapped)
    cur_zeta = zeta_ellipse_all(ii);
    cur_z = util.backwardKarmanTrefftxz(cur_z, z1_ellipse, z2_ellipse, zeta1_ellipse, zeta2_ellipse, P_ellipse);
    z_ellipse_mapped(ii) = cur_z;
end


% plot in the conformal plane 
PLOT = 0;

if PLOT 
    figure 
    plot(real(z_ellipse_mapped),imag(z_ellipse_mapped), '.--')
    axis equal
    xlabel('x')
    ylabel('y')
    title('ellipse in z plane')
end



%% create O-shaped series 
scaleFactor = util.chebyshev(51, [0,0.9]);
scaleFactor = scaleFactor(1:ceil(length(scaleFactor)/2));
laminarThickness = 0.005;
% Append fixed width extrusion to chebyshev extrusion
scaleFactor = [linspace(1,0.9,(0.1/laminarThickness))';scaleFactor(4:end)];

num_layer_r = length(scaleFactor);


PLOT = 0; % if plot 2D airfoil contour for debugging.

z_store = [];

if PLOT
    figure 
    hold on 
end
for ii = 1:length(scaleFactor)
    curScale = scaleFactor(ii);
    zeta_scaled = zeta_all;
    z_maped = zeros(length(zeta_all),1);
    for ii = 1:length(zeta_scaled)
        cur_zeta = zeta_scaled(ii);
        % scaling z1, z2 and P using current scale spacing.
        cur_z = util.backwardKarmanTrefftxz(cur_zeta, z1*curScale, z2*curScale, zeta1, zeta2, P*(1/curScale)^(1/8))+0.5*(1-curScale);
        z_maped(ii) = (1/(1/curScale)^(1/2))*cur_z;
    end
    if PLOT
        plot(real(z_maped),imag(z_maped), '.--');
    end
    z_store = [z_store; z_maped];
end

num_nodes_planar = size(z_store,1); % calculate number of nodes in a plane 

%% construct OH Grid 

l_max = 5;
l_loc = [0:1/25:l_max];
num_layer_z = length(l_loc);

all_z_loc = [];
% assemble node locations

for ii = 1:length(l_loc)
    cur_l = l_loc(ii);
    all_z_loc = [all_z_loc; [z_store, ones(length(z_store),1)*cur_l] ];
end

% extract xyz 
x = real(all_z_loc(:,1));
y = imag(all_z_loc(:,1));
z = real(all_z_loc(:,2));

PLOT = 0;

if PLOT
    figure 
    scatter3(x,y,z,'k.')
    axis equal
end



%% Calculate Mesh Connectivity 




element_connectivity_mat = util.buildElement(num_nodes_circular, num_nodes_planar, num_layer_r, num_layer_z);





%%% alternatively, build tetahedron element: 
tet_element_connectivity_mat = util.buildTetElement(element_connectivity_mat);










%% Unit Test for Karman-Treffz Shape


function nvsshape
% plot Karman-Treffz shape with different power
phi = linspace(0, 2*pi, 360);
zeta_real = -1 + 2* cos(phi);
zeta_img = 2*sin(phi);

zeta = zeta_real + 1i * zeta_img;


figure 
hold on 
subplot(1,2,1)
plot(real(zeta), imag(zeta),'.-')
xlabel('real \zeta')
ylabel('imag \zeta')
hold off

% The exponential P in Jameson's Paper
n = [-1, 0.1, 0.5, 1, 1.5, 2, 5];

subplot(1,2,2)
hold on
for ii = 1:length(n)
    curN = n(ii);
    z = curN*((zeta+1).^curN + (zeta-1).^curN)./ ...
          ( (zeta+1).^curN - (zeta-1).^curN );
    plot(real(z), imag(z),'.')
end
xlabel('real z')
ylabel('imag z')
legend('n=-1','n=0.1','n=0.5','n=1','n=1.5','n=2', 'n=5')

hold off

end


