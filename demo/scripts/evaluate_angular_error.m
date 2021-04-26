function evaluate_angular_error(Rest, Rgt, method)

    for i=1:size(Rgt,3)
        
        ang_err(i) = rot2angle( (Rest(:,:,i)*Rest(:,:,1)') * (Rgt(:,:,i)*Rgt(:,:,1)')');
        
    end
    
    %figure; plot(ang_err); title(strcat('Method: ', method, ' | Angular error w.r.t groundtruth'));
    
    fprintf('Method:%s , Mean angular error(deg) %f\n', method, mean(ang_err));
    fprintf('Method:%s , Max angular error(deg) %f\n', method, max(ang_err));    
    fprintf('Method:%s , Median angular error(deg) %f\n', method, median(ang_err));
end