function Rest = rotavg_rcd(N, input_file, output_file)
    
    % Invoke implementation in C
    cmd = sprintf('src/bin/rcd %s %s %d %d', input_file, output_file, 10000, 0);   
 
    unix(cmd);

    % Retrieve output
    fid = fopen(output_file, 'r');
    Rest = textscan(fid, '%f %f %f %f\n', N);
    Rest = quat2rotm(cell2mat(Rest));
    fclose(fid);
    
end