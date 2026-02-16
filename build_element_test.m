
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: build_element_test.m
% Author: Yiren Shen
% Description: This MATLAB script generates a hexahedral mesh and
%              an option to convert it into a tetrahedral mesh for a 3d 
%               cantilever beam. The script constructs nodes and elements 
%               based on the specified parameters. The code is used to a)
%               provide a unit test problem, and b) examine the element
%               connectivity solver.
% Date: July-19, 2023
% Version: 1.1 July-24-23:: refactorized code, use util class.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% import util functions. 
util = surf3dutil();


%%% Parameters 
% outer most layer 1
num_edge_nodes  = 11;
% number of layers in centric and z direction 
num_layer_r = 6;
num_layer_z = 12;
height_z = 6; % total length in z direction
%%%

% 5*5 nodes
x0 = [-2;-1;0;1;2;2;2;2;2;1;0;-1;-2;-2;-2;-2];
y0 = [2;2;2;2;2;1;0;-1;-2;-2;-2;-2;-2;-1;0;1];

% define n*n nodes with x,y between [-2,2] 
x0 = [linspace(-2,2,num_edge_nodes)';ones(num_edge_nodes-2,1)*2; ...
        linspace(2,-2,num_edge_nodes)';ones(num_edge_nodes-2,1)*-2];

y0 = [ones(num_edge_nodes-2,1)*2; linspace(2,-2,num_edge_nodes)'; ...
        ones(num_edge_nodes-2,1)*-2; linspace(-2,2,num_edge_nodes)'];
y0 = [y0(end); y0];
y0 = y0(1:end-1);

z0 = zeros(length(x0),1);


x = x0;
y = y0;
z = z0;
% inner 1 layer 1
x = [x; x0*0.9];
y = [y; y0*0.9];
z = [z; zeros(length(x0),1)];

% inner 2 layer 1 
x = [x; x0*0.8];
y = [y; y0*0.8];
z = [z; zeros(length(x0),1)];

% inner 3 layer 1
x = [x; x0*0.8^2];
y = [y; y0*0.8^2];
z = [z; zeros(length(x0),1)];

% inner 4 layer 1 
x = [x; x0*0.8^3];
y = [y; y0*0.8^3];
z = [z; zeros(length(x0),1)];

% inner 5 layer 1
x = [x; x0*0.8^4];
y = [y; y0*0.8^4];
z = [z; zeros(length(x0),1)];

% compute total number of nodes for a circle and one plane 
num_circular_nodes = 4*(num_edge_nodes-1);
num_planar_nodes = num_circular_nodes*num_layer_r;


% build index for laminated layers O rings
perimental_index_group_2d = reshape((1:num_planar_nodes),[num_circular_nodes,6])';
    %we know that there are 6 layers of 40 nodes each perimental rings
    % each row corresponds to a ring in O grid topology 


PLOT = 0;

if PLOT
    figure
    plot(x,y,'*-')
end

% build layer 2, 3, ... on z direction

x = repmat(x,num_layer_z,1);
y = repmat(y,num_layer_z,1);

z = zeros(length(x),1);

for ii = 1:num_layer_z
    z((ii-1)*num_planar_nodes+1:ii*num_planar_nodes) = (height_z/num_layer_z)*(ii-1);
end



PLOT = 0;
if PLOT 
    figure 
    scatter3(x,y,z,'filled')
    for ii = 1:length(x)
       % fprintf(num2str(ii))
        text(x,y,z,int2str(ii),'FontSize',14)
    end
end


% extract nodal layers: 
    %z is all nodal coordinate in format of z=[x+iy,z]
lam_layers_index = util.findLamLayersIndex(num_layer_z, num_circular_nodes, 6, [x+1i*y,z], perimental_index_group_2d);
nonZeroColumns = any(lam_layers_index, 1);
lam_layers_index = lam_layers_index(:,nonZeroColumns);


%% build elements


%%% build hexahedron element
element_connectivity_matrix = util.buildElement(num_circular_nodes, num_planar_nodes, num_layer_r, num_layer_z);

% plot the built element
PLOT = 0;
if PLOT
    figure 
    hold on
    scatter3(x,y,z,'filled')
    
    for ii = 1:size(element_connectivity_matrix,1)
        plot3(x(element_connectivity_matrix(ii,:)), y(element_connectivity_matrix(ii,:)), z(element_connectivity_matrix(ii,:)))
    end
    
    hold off 
end



%%% alternatively, build tetahedron element: 
tet_element_connectivity_mat = util.buildTetElement(element_connectivity_matrix);


% check the volume of built elements
volume_vec = util.checkTetVolume([x+1i*y,z], tet_element_connectivity_mat);

% for negative volume elements, re-order nodes index
tet_element_negative_vol_index = volume_vec < 0;
reordered_nodes = util.reorderNodes4PositiveVolume(...
    tet_element_connectivity_mat(tet_element_negative_vol_index,:), [x+1i*y,z]);
tet_element_connectivity_mat(tet_element_negative_vol_index,:) = reordered_nodes;

% compute tet elment skewness: 
% recommend skewness < 0.85, closer to 0 indicates better cell quality.
tet_skewness = util.computeTetSkewness(tet_element_connectivity_mat, [x+1i*y,z]);

% plot the built element
PLOT = 1;
if PLOT
    figure 
    hold on
    scatter3(x,y,z,'filled')
    
    for ii = 1:size(tet_element_connectivity_mat,1)
        line = plot3(x(tet_element_connectivity_mat(ii,:)), y(tet_element_connectivity_mat(ii,:)), z(tet_element_connectivity_mat(ii,:)));
        alpha(line, 0.5);
    end
    
    hold off 
end










    