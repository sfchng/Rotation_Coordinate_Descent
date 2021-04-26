function visualise_camera_pose( Rest, Cgt, type)

    N = size(Rest, 3);
    hold on;
    
    sample_size = N/50;
     
    if strcmp(type, 'Init')
        k = 1;
        for i = 1 : sample_size : N
            camera_pose = rigid3d(Rest(:,:,i)*Rest(:,:,1)', Cgt(:,i)');
            
            if k == 1 || mod(k, 10) == 0
                plotCamera('AbsolutePose', camera_pose, 'Size', 1.5, 'color', 'r', 'Label', num2str(i));
            else
                plotCamera('AbsolutePose', camera_pose, 'Size', 1.5, 'color', 'r'); 
            end
            k = k + 1;
        end
        title('Initial Rotations', 'FontSize', 13);
    else
        k = 1;
        for i = 1 : sample_size : N
            camera_pose = rigid3d(Rest(:,:,i)*Rest(:,:,1)', Cgt(:,i)');
            
            if k == 1 || mod(k, 10) == 0
                plotCamera('AbsolutePose', camera_pose, 'Size', 1.5, 'color', 'b', 'Label', num2str(i));
            else
                plotCamera('AbsolutePose', camera_pose, 'Size', 1.5, 'color', 'b'); 
            end
            k = k + 1; 
        end
        title([type, ' Estimated Rotations'], 'FontSize', 13); 
    end
     
end