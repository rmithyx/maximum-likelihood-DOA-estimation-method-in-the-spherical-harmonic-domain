%%%%%%%单频球谐域最大似然(SHMLE)算法%%%%%%%%%%%%%%%%%%%%
clear
a = 0.139;%球阵列半径
c = 340;%声速
ka = 3.142;%球Bessel函数零点对应频率
Theta_l = 90/180*pi;%声源位置
Phi_l = 180/180*pi;
load('position32.mat')%传声器位置，半径为1
M = 32;%传声器数目
N = 4;%展开阶数
MicPos = position32*a;
[Mic_Phi,Mic_Theta,Mic_R] = cart2sph(MicPos(:,1),MicPos(:,2),MicPos(:,3));%直角坐标系调整到球坐标系
Mic_Theta = pi/2-Mic_Theta;%仰角调整
%%
bn = zeros(51,1);%计算bn(ka)
for n = 0:50
    bn(n+1,1) = 4*pi*(1j)^n*SphBesselj(n,ka);
end
%%
p = zeros(M,1);%单频时，M个传声器的声压
for m = 1:M
    for n = 0:50
        Bn = bn(n+1)*eye(2*n+1);
        p(m,1) = p(m,1)+(SphHarmonic(n,Theta_l,Phi_l))'*...
            Bn*SphHarmonic(n,Mic_Theta(m),Mic_Phi(m));
    end
end
p = awgn(p,20,'measured');%添加噪声
%%
Y_nm = zeros((N+1)^2,M);%计算球谐函数
for num = 1:M
    for n = 0:N
        Y_nm((n^2+1:(n+1)^2),num) = SphHarmonic(n,Mic_Theta(num),Mic_Phi(num));
    end
end
p_nm = 4*pi/M*conj(Y_nm)*p;%信号从频域变换到球谐域
%%
Bn = zeros((N+1)^2,1);
for n = 0:N
    Bn(n^2+1:(n+1)^2) = bn(n+1);
end
L = 1;%声源数目 
a_nm = p_nm./Bn;
S = a_nm*a_nm';
[V,D] = eig(S);
[Y,I] = sort(diag(D));
E = V(:,I(1:end-L));%计算噪声子空间
%%
theta = (0:3:180)/180*pi;
phi = (0:3:360)/180*pi;
P_music = zeros(length(theta),length(phi));%MUSIC谱
for num1 = 1:length(theta)
    num1
    for num2 = 1:length(phi)
        y_nm = zeros(1,(N+1)^2);
        for n = 0:N
            y_nm(1,n^2+1:(n+1)^2) = SphHarmonic(n,theta(num1),phi(num2)).';
        end
        P_music(num1,num2) = 1/(y_nm*E*E'*y_nm');
    end
end
out = 10*log10(abs(P_music));
out = out-max(max(out));
[x,y] = find(out==max(max(out)));

figure
imagesc(phi/pi*180,theta/pi*180,out);
h = colorbar('fontsize',16);
set(get(h,'Ylabel'),'String','[dB]','Fontsize',16,'Fontname','arial')
set(gca,'yDir','normal')
set(gca,'Fontsize',16)
xlabel('{\it\phi} [deg]');ylabel('{\it\theta} [deg]');
