function [Angle, Transfer] = ReadData(AngleFile,PoseFile)
%读取角度值和对应的姿态值并简单处理
ex1 = importdata(AngleFile);
Angle = ex1.data;	//角度数据,(10个姿态*6个关节角度值)
ex2 = importdata(PoseFile);
DataCoord = ex2.data;	//姿态原始数据,(3个球*4个球面点)*(10个姿态*3个坐标xyz)
Coord = zeros(10,3,4,3);
BallCenter = zeros(10,3,3);
A = zeros(3,3);
b = zeros(3,1);
Transfer = zeros(10,4,4);	//转换矩阵,(10个姿态*(4*4矩阵))
for i = 1:12
	for j = 1:30
	//姿态，球，点，坐标
		Coord(ceil(j/3),ceil(i/4),mod((i+3),4)+1,mod((j+2),3)+1) = DataCoord(i,j);
	end
end
//计算工具球坐标系与CMM坐标系转换矩阵
for i = 1:10
	for j = 1:3
		b(1) = (Coord(i,j,1,1)^2 + Coord(i,j,1,2)^2 + Coord(i,j,1,3)^2 - Coord(i,j,2,1)^2 - Coord(i,j,2,2)^2 - Coord(i,j,2,3)^2)^2/2.0;
		b(2) = (Coord(i,j,1,1)^2 + Coord(i,j,1,2)^2 + Coord(i,j,1,3)^2 - Coord(i,j,3,1)^2 - Coord(i,j,3,2)^2 - Coord(i,j,3,3)^2)^2/2.0;
		b(3) = (Coord(i,j,1,1)^2 + Coord(i,j,1,2)^2 + Coord(i,j,1,3)^2 - Coord(i,j,4,1)^2 - Coord(i,j,4,2)^2 - Coord(i,j,4,3)^2)^2/2.0;
		A = [Coord(i,j,1,1) - Coord(i,j,2,1) Coord(i,j,1,2) - Coord(i,j,2,2) Coord(i,j,1,3) - Coord(i,j,2,3);...
			 Coord(i,j,1,1) - Coord(i,j,3,1) Coord(i,j,1,2) - Coord(i,j,3,2) Coord(i,j,1,3) - Coord(i,j,3,3);...
			 Coord(i,j,1,1) - Coord(i,j,4,1) Coord(i,j,1,2) - Coord(i,j,4,2) Coord(i,j,1,3) - Coord(i,j,4,3)];
		BallCenter(i,j,:) = A\b;
	end
	v1 = BallCenter(i,2,:) - BallCenter(i,1,:); 
	v1 = v1/norm(v1);
	v3 = cross(v1,sphere(3,:)-sphere(1,:)); 
	v3 = v3/norm(v3);
	v2 = cross(v3,v1);
	Transfer(i,:,:) = [v1(1) v1(2) v1(3) BallCenter(i,1,1);...
					   v2(1) v2(2) v2(3) BallCenter(i,1,2);...
					   v3(1) v3(2) v3(3) BallCenter(i,1,3);...
					   0 0 0 1];
end