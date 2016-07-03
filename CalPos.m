function [value] = CalPos(forward,a)
%% 根据前向运动解计算工具球坐标系在基坐标系中的位置和姿态RPY�?
% 输出为６*１的相对于基坐标系坐标矩阵x,y,z,alpha,beta,gamma
if a==0
    syms site_x site_y site_z pose_1 pose_2 pose_3;
    site_x = forward(1,4);
    site_y = forward(2,4);
    site_z = forward(3,4);
    pose_1 = atan(forward(2,1)/forward(1,1));
    pose_2 = -atan(forward(3,1)/sqrt(forward(1,1)^2+forward(2,1)^2));
    pose_3 = atan(forward(3,2)/forward(3,3));
else
    site_x = forward(1,4);
    site_y = forward(2,4);
    site_z = forward(3,4);

% pose_1 = atan(-forward(1,3)/forward(2,3));
% pose_2 = atan((forward(1,3)*sin(pose_1)-forward(2,3)...
%     *cos(pose_1))/forward(3,3));
% pose_3 = atan(-(forward(1,2)*cos(pose_1)+forward(2,2)...
%     *sin(pose_1))/(forward(1,1)*cos(pose_1)+forward(2,1)...
%     *sin(pose_1)));
    pose_2 = -atan(forward(3,1)/sqrt(forward(1,1)^2+forward(2,1)^2));
    if abs(pose_2+pi/2) <= 0.000001
        pose_2 = -pi/2;
        pose_1 = 0;
        pose_3 = -atan(forward(1,2)/forward(2,2));
    elseif abs(pose_2-pi/2) <= 0.000001
        pose_2 = pi/2;
        pose_1 = 0;
        pose_3 = atan(forward(1,2)/forward(2,2));
    else
        pose_1 = atan(forward(2,1)/forward(1,1));
        pose_3 = atan(forward(3,2)/forward(3,3));
    end
end
value = [site_x site_y site_z pose_1 pose_2 pose_3];