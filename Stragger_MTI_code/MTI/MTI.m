%--------MTI��һ��������Ӧ����-----------
clc;
clear;
Tr=0.001;                 %���ǲ������ź�����
f=0:1:2000;
y1=abs(sin(pi*f*Tr));              %k1Ϊһ�ζ���
plot(f,y1,'r');
title('MTI��һ��������Ӧ����')
xlabel('Ƶ��/Hz')
ylabel('��һ����ķ���')
text(730,0.8,'k=1,һ�ζ���(��ɫ)')
hold on;
y2=abs(sin(pi*f*Tr).*sin(pi*f*Tr));     %k=2Ϊ���ζ���
plot(f,y2,'b');
text(230,0.4,'k=2,���ζ���(��ɫ)')
y3=abs(sin(pi*f*Tr).*sin(pi*f*Tr).*sin(pi*f*Tr));
plot(f,y3,'g');               %k=3Ϊ���ζ���
text(500,0.2,'k=3,���ζ���(��ɫ)')


