function plot_angular_error(Rinit, Rest, Rgt)

    for i=1:size(Rgt,3)    
        err_est(i) = rot2angle( (Rest(:,:,i)*Rest(:,:,1)') * (Rgt(:,:,i)*Rgt(:,:,1)')');    
        err_init(i) = rot2angle( (Rinit(:,:,i)*Rinit(:,:,1)') * (Rgt(:,:,i)*Rgt(:,:,1)')');         
    end

    subplot(2,2,3);  plot(err_init, 'r');
    ylabel('Error (deg)'); xlabel('Camera ID');
    title([{'Angular Error of Initial '}, {'Rotations w.r.t. the Groundtruth'}], 'FontSize', 11);    
     
    subplot(2,2,4);  plot(err_est, 'b');
    yline(mean(err_est), 'k--', 'LineWidth', 1);
    text(10, mean(err_est+0.1), {'Mean', 'error:'}, 'FontSize', 10);
    ylabel('Error (deg)'); xlabel('Camera ID');
    title([{'Angular Error of RCD Estimated '}, {'Rotations w.r.t. the Groundtruth'}], 'FontSize', 11);    
    
end