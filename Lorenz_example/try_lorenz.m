%绘制 Lorenz 吸引子在 x-z 平面上的投影图
clc;clear;close all;
global c
c=2;
[T,Y]=ode45(@chao_SimpleLorenz,0:0.01:500,[1 2 3 ]);
figure
plot(Y(40000:end,1),Y(40000:end,3),'b')
xlim([-15 15])
xlabel('\itx')
ylabel('\itz')
grid on;
set(gca,'linewidth',0.5,'fontsize',12,'fontname','Times');             %在时间范围内使用 “ode45” 函数求解 “chao_SimpleLorenz” 的微分方程，并绘制 Lorenz 吸引子的投影图。
%--------------------------------------------------------------------------
clc;clear
global c
c=2;                                                                   %使用 “lyapunov” 函数计算 Lorenz 系统的 Lyapunov 指数，并将结果存储在 “Res” 变量中。最后打印输出了最终的 Lyapunov 指数。

[T,Res]=lyapunov(3,@SimLorenz_ly,@ode45,0,0.5,500,[ 0.1 0 0.1],5);          %最后打印输出了最终的 Lyapunov 指数。
Ly=Res(end,:)
%----------------------Lyapunov 指数用于描述混沌系统中的指数敏感性-------------------------------
clc
clear
global c
% C=-2:0.05:8;L=length(C);

L=250;
C=linspace(-2,8,L);
ly=zeros(3,L);
Com=zeros(1,L);
figure
for i=1:L
    c=C(i);
    
    [T,Y]=ode45(@chao_SimpleLorenz,0:0.01:100,[1 2 3]);
    data=Y(5000:end,3);
    Com(i)=COFuZadu(data,15);
    for j=2:(length(data)-1)
        if data(j)>data(j-1)&&data(j)>data(j+1)
            plot(c,data(j),'.r','markersize',1);
            hold on;
            if j==20
                break;
            end
        end
    end
    
   [T,Res]=lyapunov(3,@SimLorenz_ly,@ode45,0,0.5,100,[ 0.1 0 0.1],5);
    ly(:,i)=Res(end,:);
    disp(i)				
end
xlabel('\itc')
ylabel('{\itz}_{max}')
grid on;
set(gca,'linewidth',0.5,'fontsize',12,'fontname','Times');
%------------------绘制参数 c 对系统行为的影响---------------------

figure
plot(C,ly(1,:),'r','linewidth',1)
hold on
plot(C,ly(2,:),'k','linewidth',1)
plot(C,ly(3,:),'b','linewidth',1)
ylim([-15,2])
xlabel('\itc')
ylabel('Lyapunov')
grid on;
set(gca,'linewidth',0.5,'fontsize',12,'fontname','Times');


figure
plot(C,Com,'r','linewidth',1)
xlabel('\itc')
ylabel('C_0')
grid on;
set(gca,'linewidth',0.5,'fontsize',12,'fontname','Times');


