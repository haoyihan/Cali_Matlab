%% 用于计算加入基坐标系到CMM坐标系转换矩阵校准的模型偏微�?
function [dC] = Partial2(SphereTrans)
%% 计算转换矩阵的偏微分矩阵
% 输出为转换矩阵关于各个关节的四个参数的偏微分
%% 偏微分计算顺序为p,a,t,d,冗余参数为d2，d6和t6
syms a1 p1 d1 t1;
syms a2 p2 d2 t2;
syms a3 p3 d3 t3;
syms a4 p4 d4 t4;
syms a5 p5 d5 t5;
syms a6 p6 d6 t6;
syms x y z a b g;
% 分别为a，alpha，d，theta

%�?��建从CMM坐标系到基坐标系的转换矩�?
T0 = Rotz(a)*Roty(b)*Rotx(g)+transl(x,y,z)-...
    [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1]; 
T1 = transl(0,0,d1)*Rotz(t1)*transl(a1,0,0)*Rotx(p1);
T2 = transl(0,0,d2)*Rotz(t2)*transl(a2,0,0)*Rotx(p2);
T3 = transl(0,0,d3)*Rotz(t3)*transl(a3,0,0)*Rotx(p3);
T4 = transl(0,0,d4)*Rotz(t4)*transl(a4,0,0)*Rotx(p4);
T5 = transl(0,0,d5)*Rotz(t5)*transl(a5,0,0)*Rotx(p5);
T6 = transl(0,0,d6)*Rotz(t6)*transl(a6,0,0)*Rotx(p6);
T7 = SphereTrans;

% 计算CMM坐标系到工具球坐标系的转换矩�?
T = T0*T1*T2*T3*T4*T5*T6*T7;

% 转化为（x,y,z,psi,theta,phi）后分别求偏�?
[C] = CalPos(T,0);

% 计算对各个参数的偏导
% m×n,m是坐标参数，n是每个连杆的DH参数
dC = jacobian(C,[p1;a1;t1;d1;p2;a2;t2;p3;a3;t3;d3;p4;a4;t4;d4;p5;a5;t5;d5;p6;a6;x;y;z;a;b;g]);
end