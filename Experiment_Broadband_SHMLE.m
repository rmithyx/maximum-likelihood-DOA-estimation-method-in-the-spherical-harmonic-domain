%%%%%%%%%%%%%%%%%%球阵列最大似然(SHMLE)声源定位算法%%%%%%%%%%%%%%%%%%%
% 传声器摆放位置为：           5----------8
%                          6/---------7/ |
%                          |  |       |  |
%                          |  1-------|--4
%                          2/---------3 /
%
%          坐标系：         ^
%                       z  |
%                          |
%                          | ------>    y
%                         / o
%                        /
%                   x   ~
%%
clear
MicPos = [ -0.08   -0.08   -0.08 %传声器单元直角坐标系位置
    0.08   -0.08   -0.08
    0.08    0.08   -0.08
    -0.08    0.08   -0.08
    -0.08   -0.08    0.08
    0.08   -0.08    0.08
    0.08    0.08    0.08
    -0.08    0.08    0.08
    ];
a = 0.08*sqrt(3);%球阵列半径a
[Mic_Phi,Mic_Theta,Mic_R] = cart2sph(MicPos(:,1),MicPos(:,2),MicPos(:,3));%直角坐标系调整到球坐标系
Mic_Theta = pi/2-Mic_Theta;%仰角调整
M = 8;%传声器数目
c = 340;%声速
K = 2048;%数据帧长度
N = 1;%球阵列展开阶数
%%
load('w_cal.mat');%传声器校准文件
% [data,fs] = audioread('ListeningRoom_Whitenoise.wav');%视听室实验数据
[data,fs] = audioread('AnechoicChamber_Whitenoise.wav');%消声室实验数据
for num = 1:M
    data(:,num) = filter(w_cal(:,num),1,data(:,num));%分别对每个传声器的数据进行校准
end
FrameNumber = floor(length(data(:,1))/K);%将接收到的数据分帧
FrameFlag = zeros(FrameNumber,1);%数据有效性标记
cnt = 0;
for num = 1:FrameNumber
    if(sum(data((num-1)*K+1:num*K,1).^2)>10*1e-5)%判定数据是否有效
        FrameFlag(num) = 1;
        cnt = cnt+1;
        flags(cnt) = (num-1)*K+1;%有效帧数据起始位置
    else
        FrameFlag(num) = 0;
    end
end
%%
Fu = round(K*c*N/(fs*2*pi*a));%计算频率上下限 ka ~ [0.5 1]
Fl = round(Fu/2);
%%
X = zeros(K,M);%将信号变换到频域
for m = 1:M
    x_p(:,m) = data([1:K]+flags(4),m);
    X(:,m) = fft(x_p(:,m));
end
%%
Bn = zeros((N+1)^2,Fu);%计算bn(ka)
for k = 1:Fu
    ka = 2*pi*k/K*fs/c*a;
    for n = 0:N
        Bn(n^2+1:(n+1)^2,k) = 4*pi*(1j)^n*SphBesselj(n,ka);
    end
end
%%
Y_nm = zeros((N+1)^2,M);%计算球谐函数
for n = 0:N
    for m = 1:M
        Y_nm(n^2+1:(n+1)^2,m) = SphHarmonic(n,Mic_Theta(m),Mic_Phi(m));
    end
end
X_nm = zeros(Fu,(N+1)^2);%不同频率球谐变换结果
for k = 1:Fu
    X_nm(k,:) = 4*pi/M*X(k+1,:)*Y_nm';
end
X_nm = X_nm.';

%%
theta = (0:3:180)/180*pi;%划分空间网格
phi = (0:3:360)/180*pi;
for num1 = 1:length(theta)
    num1
    for num2 = 1:length(phi)
        temp = 0;
        P_nm = zeros((N+1)^2,1);
        for n = 0:N 
            P_nm(n^2+1:(n+1)^2,1) = SphHarmonic(n,theta(num1),phi(num2));
        end
        
        for k = Fl:Fu
            d_nm = diag(Bn(:,k))*conj(P_nm);%计算dnm
            temp = temp+norm(X_nm(:,k)-d_nm*pinv(d_nm)*X_nm(:,k))^2;%论文公式(19)
        end
        Out(num1,num2) = temp;
    end
end

Out1 = -10*log10(Out);
imagesc(phi/pi*180,theta/pi*180,Out1-max(max(Out1)))
colorbar;set(gca,'ydir','normal')
% caxis([-17 0])
set(gca,'Fontsize',18)
xlabel('{\it\phi} [deg]');ylabel('{\it\theta} [deg]');zlabel('[dB]')






