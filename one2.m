%% 用于试验加入基坐标系到CMM坐标系转换矩阵校准的模型精度问题
%% 实现�?��迭代计算，共测量10个点，每个点测量3*4*3个数据，�?��测量360个数�?% 每次迭代共用�?��数据，不�?��反复测量
%% 计算顺序为p,a,t,d
% By zhangyanbo@rokae.com
% 2016.04.12
%clearvars -except result result1;
% clc;
function [result,error] = one2(dC,L_t,p_n,a_n,d_n,theta_e,Sample,pos,cmm_n,cmm_t,SphereTrans)

% 根据DH参数创建Link对象
%       theta    d        a             alpha
L_n(1) = Link('d', d_n(1), 'a', a_n(1), 'alpha', p_n(1));
L_n(2) = Link('d', d_n(2), 'a', a_n(2), 'alpha', p_n(2));
L_n(3) = Link('d', d_n(3), 'a', a_n(3), 'alpha', p_n(3));
L_n(4) = Link('d', d_n(4), 'a', a_n(4), 'alpha', p_n(4));
L_n(5) = Link('d', d_n(5), 'a', a_n(5), 'alpha', p_n(5));
L_n(6) = Link('d', d_n(6), 'a', a_n(6), 'alpha', p_n(6));

syms a1 p1 d1 t1 a2 p2 d2 t2 a3 p3 d3 t3 a4 p4 d4 t4 ...
    a5 p5 d5 t5 a6 p6 d6 t6 x y z a b g;
% 参数误差矩阵
%  alpha    a      d    theta
%  0.005  0.003  0.005  0.003
%  0.005  0.003 (0.005) 0.003
%  0.005  0.003  0.005  0.003
%  0.005  0.003  0.005  0.003
%  0.005  0.003  0.005  0.003
%  0.005  0.003 (0.005)(0.003)

% 参�?位姿
% qz = [0 0 0 0 0 0]; % 零角�?% q0 = [0 -pi/2 pi 0 0 0]; % 初始位型
% qstretch = q0 +[0 pi/2 -pi/2 0 0 0]; %伸直
% qt = q0 + [0 0 -pi/2 0 0 0]; % z方法伸直

%�?��接连杆，机器人本�?
xb06_n = SerialLink(L_n,'name', 'XB06_n', 'manufacturer', 'Rokae');    
xb06_t = SerialLink(L_t,'name', 'XB06_t', 'manufacturer', 'Rokae');    

% 工具球相对工具球坐标系的位置，假设３个球心组成一个直角，且球1在轴6�?
% sphere = zeros(3,3);
% sphere(1,:) = [0,0,0];  %球１坐标，坐标系原点
% sphere(2,:) = [0.05,0,0];
% sphere(3,:) = [0,0.05,0];

%xb06_t.display();    %显示D-H参数�?
Jc = zeros(6*Sample,27);
Ec = zeros(6*Sample,1);
    
for k = 1:Sample
    % 加上工具球装置，计算前向运动�?    % 改变每次的位�?%% pos(k,j)
    t_n(1) = -1+pi/9*pos(k,1);
    t_n(2) = 1+pi/5*pos(k,2);
    t_n(3) = 10-pi/4*pos(k,3);
    t_n(4) = 1+pi/6*pos(k,4);
    t_n(5) = 3-pi/9*pos(k,5);
    t_n(6) = 5-pi/8*pos(k,6);
    % CMM坐标系到基坐标系名义转换
    Tc_n = Rotz(cmm_n(4))*Roty(cmm_n(5))*Rotx(cmm_n(6))+...
        transl(cmm_n(1),cmm_n(2),cmm_n(3))-...
        [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1]; 
    % 名义前向运动�?    
    forward_n = Tc_n*xb06_n.fkine(t_n)*SphereTrans; %前向运动�?    % CMM坐标系到基坐标系实际转换
    Tc_t = Rotz(cmm_t(4))*Roty(cmm_t(5))*Rotx(cmm_t(6))+...
        transl(cmm_t(1),cmm_t(2),cmm_t(3))-...
        [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1];
    % 实际角度
    theta_t = t_n+theta_e;
    %�?��际前向运动解
    forward_t = Tc_t*xb06_t.fkine(theta_t)*SphereTrans; 
    % xb06_n.plot(q0,'jaxes');

    % 根据前向运动解计算工具球坐标系在基坐标系中的位置和姿�?    % 名义位姿
    [value_n] = CalPos(forward_n,1);
    % 实际位姿
    [value_t] = CalPos(forward_t,1);

    % 位姿误差
    error = value_t'-value_n';
    
    % 对偏导矩阵赋�?    
    J = double(vpa(subs(dC,{a1,p1,d1,t1,a2,p2,d2,t2,a3,p3,d3,t3,a4,p4,d4,t4,...
        a5,p5,d5,t5,a6,p6,d6,t6,x,y,z,a,b,g},{a_n(1),p_n(1),d_n(1),t_n(1),...
        a_n(2),p_n(2),d_n(2),t_n(2),a_n(3),p_n(3),d_n(3),t_n(3),...
        a_n(4),p_n(4),d_n(4),t_n(4),a_n(5),p_n(5),d_n(5),t_n(5),...
        a_n(6),p_n(6),d_n(6),t_n(6),cmm_n(1),cmm_n(2),cmm_n(3),cmm_n(4),cmm_n(5),cmm_n(6)})));
    
    Jc(6*k-5:6*k,:) = J;
    Ec(6*k-5:6*k,1) = error;
end

% result = (Jc'*Jc)^-1*Jc'*Ec
% 奇异值分解，解决矩阵近似奇异问题
% [U,S,V] = svd(Jc);
% result = V/S*U'*Ec;
[Q,R] = qr(Jc);
result = pinv(R)*Q'*Ec
% 求解参数误差，顺序为alpha,a,theta,d