%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                 ROTATION COORDINATE DESCENT
%
%
%  This package contains the source code which implements the
%  Rotation Coordinate Descent (RCD and RCDL) in
%
%                 Rotation Coordinate Descent for 
%             Fast Globally Optimal Rotation Averaging
%            
%
%  The source code and demo are suplied for academic use only.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
close all;

format long;

% path
addpath 'scripts/';

N_a = [ 1000, 2000, 3000];
r = 50; Sigma = 0.1; Density = 0.4;
fig_p = [ 0, 500, 1000];
tic;
for i = 1 : numel(N_a)
    N = N_a(i); 
    
    fprintf('Status: Loading viewgraph : Number of cameras %d, Graph Density %d\n', N, Density);

    load (strcat('./data/N_', num2str(N), '_data.mat'));
    input_file =  strcat('./data/N_', num2str(N), '.txt');
    output_file = 'tmp_out.txt';

    fprintf('Status: Calling RCD\n');   
    
    % Run RCD    
    Rest_RCD = rotavg_rcd(N, input_file, output_file);
    
    fprintf('\n');
    fprintf('Status: Plotting camera poses computed using estimated rotations\n');
    figure('Renderer', 'painters', 'Position', [fig_p(i) 0 700 700]);        
      
    subplot(2,2,1); 
    visualise_camera_pose(Rinit, Cgt, 'Init');
    subplot(2,2,2); 
    visualise_camera_pose(Rest_RCD, Cgt, 'RCD');  
    plot_angular_error(Rinit, Rest_RCD, Rgt);
    sgtitle( {['Number of cameras: ', num2str(N)], 'Note : The cameras have been subsampled to 50 for a better visualisation.'}); 
    drawnow
    
    fprintf('---------------------------------------------------------------------------------------------------------------------------------------------\n');
    fprintf('\n'); 
end

fprintf('The RCD demo ends in %fs\n', toc);



