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
scenarios = {'torus3D' , 'grid3D'};

fig_pos = [ 0, 1000];
tic;
for i = 1 : numel(scenarios)
   
    scenario = scenarios{i};
    fprintf('Status: Running RCDL for SLAM viewgraph :%s\n', scenario);
    input_file = strcat('./data/', scenario, '.txt');
    output_file = 'tmp_out.txt';
    load (strcat('./data/', scenario, '_poses_demo.mat'));  
    
    N = size(poses,2);

    % Run RCDL
    Rest_RCDL = rotavg_rcdl( N, input_file, output_file);
    
    fprintf('\n');
    fprintf('Status: Plotting camera poses\n');
    figure('Renderer', 'painters', 'Position', [fig_pos(i) 500 600 600]);   
    subplot(1,2,1); visualise_camera_pose_g2o(Rinit, 'Init', poses, scenario);
    subplot(1,2,2); visualise_camera_pose_g2o(Rest_RCDL, 'RCDL', poses, scenario);  
    sgtitle({['RCDL on ', scenario], 'Note: We only solve for rotations.'}, 'FontSize', 16);
    drawnow
   
end
fprintf('Demo ends in %f\n', toc);
