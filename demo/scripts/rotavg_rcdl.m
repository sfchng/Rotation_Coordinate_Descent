function Rest = rotavg_rcdl(N, input_file, output_file)

    % Invoke implementation in C
    fprintf("Status: Running RCDL in C++\n");
    cmd = sprintf('src/bin/rcdl %s %s L2 5 %d %d %d', input_file, output_file, 10000,  100000);        
    unix(cmd);

    % Retrieve output
    fid = fopen(output_file, 'r');
    Rest = textscan(fid, '%f %f %f %f\n', N);
    Rest = quat2rotm(cell2mat(Rest));
    fclose(fid);
    
end