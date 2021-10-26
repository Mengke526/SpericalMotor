%% 球形电机驱动电流计算
%% 1.数据清空
clc,clear
%% 2.程序计时开始
tic
%% 3.球形电机当前欧拉角设置---------------------
alpha = 0*(pi/180);
beta = 0*(pi/180);
gamma = 0*(pi/180);
%% 4.球形电机期望力矩大小设置-------------------
Td = [0;0;0];
%% 5.旋转矩阵R及其偏导a、b、c的计算
Sa = sin(alpha); Ca = cos(alpha);
Sb = sin(beta); Cb = cos(beta);
Sg = sin(gamma); Cg = cos(gamma);
R = [Cg*Cb,-Sg*Cb,Sb;
    Sg*Ca+Cg*Sb*Sa,Cg*Ca-Sg*Sb*Sa,-Cb*Sa;
    Sg*Sa-Cg*Sb*Ca,Cg*Sa+Sg*Sb*Ca,Cb*Ca];
a = [0,0,0;
    -Sg*Sa+Cg*Sb*Ca,-Cg*Sa-Sg*Sb*Ca,-Cb*Ca;
    Sg*Ca+Cg*Sb*Sa,Cg*Ca-Sg*Sb*Sa,-Cb*Sa];
b = [-Cg*Sb,Sg*Sb,Cb;
    Cg*Cb*Sa,-Sg*Cb*Sa,Sb*Sa;
    -Cg*Cb*Ca,Sg*Cb*Ca,-Sb*Ca;];
c = [-Sg*Cb,-Cg*Cb,0;
    Cg*Ca-Sg*Sb*Sa,-Sg*Ca-Cg*Sb*Sa,0;
    Cg*Sa+Sg*Sb*Ca,-Sg*Sa+Cg*Sb*Ca,0];
%% 6.PM、EM数量及位置参数的设置计算-----------
numPM = 24; % PM的数量
numEM = 12;% EM的数量
ThetaPM = (pi/180).*[105.*ones(1,numPM/2),75.*ones(1,numPM/2)]; % PM在球坐标系中的仰角Theta
PhiPM = 15*pi/180 + [2*pi/(numPM/2).*linspace(0,numPM/2-1,numPM/2),2*pi/(numPM/2).*linspace(0,numPM/2-1,numPM/2)]; %PM在球坐标系中的方位角Phi
ThetaEM = (pi/180).*[115.*ones(1,2*numEM/3),90.*ones(1,numEM/3)]; % EM在球坐标系中的仰角Theta
PhiEM = [2*pi/(2*numEM/3).*linspace(0,2*numEM/3-1,2*numEM/3),2*pi/(2*numEM/3).*linspace(0,numEM/3-1,numEM/3)+pi/(2*numEM/3)]; % EM在球坐标系中的方位角Phi
epi = zeros(3,numPM); % PM初始（电机欧拉角为0时）单位方向向量
for i = 1:numPM
    epi(:,i) = [sin(ThetaPM(i))*cos(PhiPM(i));sin(ThetaPM(i))*sin(PhiPM(i));cos(ThetaPM(i))];
end
eri = R*epi; % PM单位方向向量
esj = zeros(3,numEM); % EM单位方向向量
for j = 1:numEM
    esj(:,j) = [sin(ThetaEM(j))*cos(PhiEM(j));sin(ThetaEM(j))*sin(PhiEM(j));cos(ThetaEM(j))];
end
% 第i个PM与第j个EM的sigma角计算------------
sigma = zeros(numPM,numEM);
for i = 1:numPM
    for j = 1:numEM
        sigma(i,j) = (180/pi)*acos(dot(eri(:,i),esj(:,j)));
    end
end
%% 7.PM磁化方向及大小的设置--------------------------
uM0 = 1.465; % 剩余磁化强度 T
lambda = uM0.*[(-1).^(linspace(1,numPM/2,numPM/2)),(-1).^(linspace(numPM/2,numPM-1,numPM/2))]; % PM的磁化方向
%% 8.K的计算--------------------------------------
Thetaij = zeros(numPM,numEM);
K = zeros(3,numEM);
for j = 1:numEM
    for i = 1:numPM
        Thetaij(i,j) = lambda(i)*dfLambdaFit(sigma(i,j))/sqrt(1-dot(eri(:,i),esj(:,j))^2);
        K(1,j) = K(1,j)+Thetaij(i,j)*(a*epi(:,i))'*esj(:,j);
        K(2,j) = K(2,j)+Thetaij(i,j)*(b*epi(:,i))'*esj(:,j);
        K(3,j) = K(3,j)+Thetaij(i,j)*(c*epi(:,i))'*esj(:,j);
    end
end
%% 9.对电机功率要求最小的电流最优解------------
u = K'/(K*K')*Td;
u = u/2;
u0 = [u(1:8);-u(5:8);-u(1:4);u(9:12);-u(9:12)];
%% 10.程序计时结束----------------------------------
toc
