%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mesherDemo_oml.m
% Demo on constructing O-grid for NACA 0012 Airfoil inside structural
% domain. 
% Creare only the OML surface grid at mid-plane for multiple layers
% Author: Yiren Shen 
% Date: August-28-2023
% Version: 0.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% import utillity function
util = surf3dutil();




% cal NACA 0012 Nodes 
num_nodes_circular = 33;
num_layer_r = 1;
num_layer_z = 21; 




[x_up, y_up, x_lo, y_lo] = util.drawNACA([0 0 1 2], ceil(num_nodes_circular/2));

% modify TE of generated data 
x_up(1) = 1.01;
x_lo(1) = 1.01;

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

% concatnate nodes, 0th node start from TE

flipped_xup = flip(x_up);
flipped_yup = flip(y_up);

x_circ = [x_lo; flipped_xup(2:end)];
y_circ = [y_lo; flipped_yup(2:end)];

% remove the repeated tail 
x_circ = x_circ(1:end-1);
y_circ = y_circ(1:end-1);

num_nodes_circular = length(x_circ);
num_nodes_planar = num_nodes_circular;




% construct OH grid
node_xyz = [];
l_max = 4;
l_loc = linspace(0, l_max, num_layer_z);

node_z = [];
node_x = [];
node_y = [];

for ii = 1:numel(l_loc)
    cur_z = repmat(l_loc(ii), num_nodes_planar,1);
    node_z = [node_z; cur_z];
    node_x = [node_x; x_circ];
    node_y = [node_y; y_circ];
end


node_xyz = [node_x, node_y, node_z];

if PLOT
    figure 
    hold on
    scatter3(node_x,node_y,node_z,'k.')
    for ii = 1:length(node_x)
        text(node_x(ii), node_y(ii), node_z(ii), num2str(ii));
    end
    
    axis equal
end

% build element connectivity 

element_connectivity_mat = [];
for ii = 1:num_layer_z -1
    for jj = 1:num_nodes_planar
        if jj <= num_nodes_planar -1 
            cur_node_idx = num_nodes_planar*(ii-1) + (jj-1) +1;
            node_right = cur_node_idx +1 ;
            node_north = cur_node_idx + num_nodes_planar;
            node_north_right = node_north +1;
            element_connectivity_mat = [element_connectivity_mat;
                                [cur_node_idx, node_right, node_north_right,node_north]];
        else
            cur_node_idx = num_nodes_planar*(ii-1) + (jj-1) + 1;
            node_right = num_nodes_planar*(ii-1) +1;
            node_north = cur_node_idx + num_nodes_planar;
            node_north_right = node_right + num_nodes_planar;
            element_connectivity_mat = [element_connectivity_mat;
                                [cur_node_idx, node_right, node_north_right,node_north]];
        end

    end
end










