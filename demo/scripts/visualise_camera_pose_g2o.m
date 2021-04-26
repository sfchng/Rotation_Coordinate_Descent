function visualise_camera_pose_g2o(Rest_RCDL, type, poses, scenario)

    N = size(Rest_RCDL,3);
    
    for i = 1 : N
        
        Tgt(:,i) = -poses(i).R' * poses(i).t;
        camera_location(:,i) = - (Rest_RCDL(:,:,i) * Rest_RCDL(:,:,1)')' * Tgt(:,i); 
        
    end
    
    if strcmp(type, 'Init')
        plot3(camera_location(1,:), camera_location(2,:), camera_location(3,:), 'r.');
        title('Initial 6DoF Poses', 'FontSize', 13);
    else
        plot3(camera_location(1,:), camera_location(2,:), camera_location(3,:), 'b.');
        title({'6DoF Poses with','RCDL Estimated Rotations'}, 'FontSize', 13);
    end
    
    if strcmp(scenario, 'grid3D')
        axis equal; view([55 30]);
        margin = 10;
        axis([-margin 20+margin -margin 20+margin 1 20+margin])  
        axis off
    end
    
    if strcmp(scenario, 'torus3D')
        axis equal; view([55 30]);
        axis equal;
        axis off;
    end    
  
end