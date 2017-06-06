function [ CRB_result ] = JS_CRB_PolyPhase( I, N, Ts, SNR )
%CRB_ �˴���ʾ�йش˺�����ժҪ
% % ����ʽ��λ�źŵ�CRB
% % I���״�,0:���ࡢ1:����Ƶ�ʡ�2:��Ƶ�ʡ�....
% % N������ʱ���źŵ���
% % Ts�� �������
% % SNR�������
% % CRB_result�� ÿһ���ǲ�ͬ������µĵ���������������
D = diag (Ts.^(0:I));
n0 = -N/2 + 1/2;
n =n0 : N/2 - 1/2;
H = zeros(I,I);
for k = 0:I
    for l = 0:I
        H(k+1,l+1) = sum(n.^(k+l));
    end
end
SNR = 10.^(-SNR./10);
inv_D = inv(D);
inv_H = inv(H);
CRB_result = zeros(I+1, length(SNR));
for i = 1:length(SNR)
    CRB_result(:,i) = diag(0.5 * SNR(i) * inv_D * inv_H * inv_D);
end
end

