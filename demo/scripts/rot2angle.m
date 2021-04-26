function [theta] = rot2angle(R)
%finds the rotation angle of a rotation matrix
% theta=(180/pi)*acos(min(1, max(-1, 1-0.5*abs(trace(eye(3)-R)))));

theta=(180/pi)*acos(min(1, max(-1, 1-0.25*trace((eye(3)-R)'*(eye(3)-R))))); %%% it turns out that this formula is numerically more stable than the above formula

end