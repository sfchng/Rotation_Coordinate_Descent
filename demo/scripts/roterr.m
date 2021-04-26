function [R] = roterr(err)
% rotation matrix of err ang error about a random axis of rotation

R = vrrotvec2mat([2*rand(1,3)-1, err*randn()] );
end