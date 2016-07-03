%% 用于试验加入基坐标系到CMM坐标系转换矩阵校准的模型
%% 利用三坐标测量仪的机器人连杆参数标定算法仿真
%% 计算顺序为p,a,t,d
% By zhangyanbo@rokae.com
% 2016.04.12
clear; 
clc;

%% 测量系统专用
% 计算工具球坐标系相对基坐标系转换矩阵，此时工具球坐标为CMM坐标系中坐标
% v1 = sphere(2,:)-sphere(1,:); v1 = v1/norm(v1);
% v3 = cross(v1,sphere(3,:)-sphere(1,:)); v3 = v3/norm(v3);
% v2 = cross(v3,v1);
% I = (v1;v2;v3);
% I_0 = (v1_0;v2_0;v3_0);   % 基坐标系原点在CMM坐标系中的坐标向�?% transl = (x-x_0,y-y_0,z-z_0);
% M = I_0^-1*(I_0-transl);

%% 
% 实际参数矩阵
%  alpha    a    theta     d
%  0.005  0.033  0.003   0.385
%  0.005  0.343  0.003  (0.005) 
%  0.005 -0.032  0.003   0.005  
%  0.005  0.003  0.003   0.350  
%  0.005  0.003  0.003   0.005  
%  0.005  0.003 (0.003) (0.005)

% 采样点数
Sample = 10;

% 基坐标系到CMM坐标系的转换参数RPY名义�?
cmm_n = zeros(1,6);
cmm_n(1) = 0.5; %x_n
cmm_n(2) = 0.5; %y_n
cmm_n(3) = 0.5; %z_n
cmm_n(4) = 0.01;    %alpha_n
cmm_n(5) = 0.01;    %beta_n
cmm_n(6) = 0.01;    %gamma_n
% 实际�?偏差值均�?.001
cmm_t = zeros(1,6);
cmm_t(1) = 0.501;
cmm_t(2) = 0.501;
cmm_t(3) = 0.501;
cmm_t(4) = 0.011;
cmm_t(5) = 0.011;
cmm_t(6) = 0.011;

SphereTrans = [1 0 0 0;0 1 0 0;0 0 1 0.2;0 0 0 1];  %工具球装置相对工具转换矩�?

% 连杆长度 
L12z = 380; %�?��轴z方向上的长度
L12x = 30;  %�?��轴x方向上的长度
L23z = 340; %二三轴z方向上的长度
L34z = 35;  %三四轴z方向上的长度
L34x = 120; %三四轴x方向上的长度
L45x = 225; %四五轴x方向上的长度
L56x = 90;  %五六轴x方向上的长度

% 转为DH并转换单位，�?
d_n(1) = L12z/1000; %link offset                         0.38
a_n(1) = L12x/1000; %link length                         0.03
a_n(2) = L23z/1000; %                                    0.34
a_n(3) = -L34z/1000; % 沿X_3, 从Z_3移动到Z_4              -0.035
d_n(4) = (L34x+L45x)/1000; % 沿轴4，从X_3移动到X_4的距�?    0.345
d_n(6) = L56x/1000; % 沿轴6                               0.09

p_n(1) = -pi/2;
d_n(2) = 0;
p_n(2) = 0;
d_n(3) = 0;
p_n(3) = pi/2;
a_n(4) = 0;
p_n(4) = -pi/2;
d_n(5) = 0;
a_n(5) = 0;
p_n(5) = pi/2;
a_n(6) = 0;
p_n(6) = 0;

t_n = [-1 1 10 1 3 5];
% 角度误差
theta_e = [0.003,0.003,0.003,0.003,0.003,0.003];
result = zeros(1,27);

% 实际的Link对象，不会随迭代变化
L_t(1) = Link('d', d_n(1)+0.005, 'a', a_n(1)+0.003, 'alpha', p_n(1)+0.005);
L_t(2) = Link('d', d_n(2)+0.005, 'a', a_n(2)+0.003, 'alpha', p_n(2)+0.005);
L_t(3) = Link('d', d_n(3)+0.005, 'a', a_n(3)+0.003, 'alpha', p_n(3)+0.005);
L_t(4) = Link('d', d_n(4)+0.005, 'a', a_n(4)+0.003, 'alpha', p_n(4)+0.005);
L_t(5) = Link('d', d_n(5)+0.005, 'a', a_n(5)+0.003, 'alpha', p_n(5)+0.005);
L_t(6) = Link('d', d_n(6)+0.005, 'a', a_n(6)+0.003, 'alpha', p_n(6)+0.005);

% 改变每次的位�?%% pos(k,j)
pos = 2*rand(Sample,6)-1;

% 计算各个参数的偏�?
[dC] = Partial2(SphereTrans);

error = zeros(1,6);

% 对冗余参数进行估计，提高辨识精度
% d_n(2) = d_n(2)+0.006;
% d_n(6) = d_n(6)+0.006;
% theta_e(6) = theta_e(6)-0.004;

for k = 1:4
    coord0 = sqrt(error(1)^2+error(2)^2+error(3)^2);
    p_n(1) = p_n(1)+result(1);
    a_n(1) = a_n(1)+result(2);
    theta_e(1) = theta_e(1)-result(3);
    d_n(1) = d_n(1)+result(4);
    p_n(2) = p_n(2)+result(5);
    a_n(2) = a_n(2)+result(6);
    theta_e(2) = theta_e(2)-result(7);
    p_n(3) = p_n(3)+result(8);
    a_n(3) = a_n(3)+result(9);
    theta_e(3) = theta_e(3)-result(10);
    d_n(3) = d_n(3)+result(11);
    p_n(4) = p_n(4)+result(12);
    a_n(4) = a_n(4)+result(13);
    theta_e(4) = theta_e(4)-result(14);
    d_n(4) = d_n(4)+result(15);
    p_n(5) = p_n(5)+result(16);
    a_n(5) = a_n(5)+result(17);
    theta_e(5) = theta_e(5)-result(18);
    d_n(5) = d_n(5)+result(19);
    p_n(6) = p_n(6)+result(20);
    a_n(6) = a_n(6)+result(21);
    cmm_n(1) = cmm_n(1)+result(22);
    cmm_n(2) = cmm_n(2)+result(23);
    cmm_n(3) = cmm_n(3)+result(24);
    cmm_n(4) = cmm_n(4)+result(25);
    cmm_n(5) = cmm_n(5)+result(26);
    cmm_n(6) = cmm_n(6)+result(27);
    [result,error] = one2(dC,L_t,p_n,a_n,d_n,theta_e,Sample,pos,cmm_n,cmm_t,SphereTrans);
    coord1 = sqrt(error(1)^2+error(2)^2+error(3)^2);
%    result'
%    error'
    disp('位置误差精度变化');
    Res_Improve = coord0-coord1
end