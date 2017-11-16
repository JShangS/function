%--------MTI归一化幅度响应曲线-----------
clc;
clear;
Tr=0.001;                 %三角波调制信号周期
f=0:1:2000;
y1=abs(sin(pi*f*Tr));              %k1为一次对消
plot(f,y1,'r');
title('MTI归一化幅度响应曲线')
xlabel('频率/Hz')
ylabel('归一化后的幅度')
text(730,0.8,'k=1,一次对消(红色)')
hold on;
y2=abs(sin(pi*f*Tr).*sin(pi*f*Tr));     %k=2为二次对消
plot(f,y2,'b');
text(230,0.4,'k=2,二次对消(蓝色)')
y3=abs(sin(pi*f*Tr).*sin(pi*f*Tr).*sin(pi*f*Tr));
plot(f,y3,'g');               %k=3为三次对消
text(500,0.2,'k=3,三次对消(绿色)')


