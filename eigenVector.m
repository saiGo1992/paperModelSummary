%����ʸ��ͼ
% ����Ƶ�ʽ�Ͳ�ʸ����������ʸ�������Ͳ�ʸ�Ĺ�ϵ
function eigenVector(omega,waveVector,M,m,K1,K2,w_0_loss,w_0_gain)
for i=1:length(omega)
    if imag(omega(i))~=0
       continue;
    end
    %��Ч������ϵ������
    M_eff_loss2 = M+m*w_0_loss/(w_0_loss^2-omega(i)^2);
    M_eff_gain2 = M+m*w_0_gain/(w_0_gain^2-omega(i)^2);
    matrix = [K1+K2-omega(i)^2*M_eff_loss2,      -(K2*exp(-1i*waveVector(i))+K1);
               K1+K2*exp(1i*waveVector(i)),      -(K1+K2-omega(i)^2*M_eff_gain2)];
    [V,D]=(eig(matrix));
    double(V);
    plot(waveVector(i),real(V(1,:)),'g.');
    hold on;
end














%�½���ͼ����
% figure;
% for c=1:s(2)
%     %���ݹ�һ������
% %     if real(kw(1,c))~=0 && imag(kw(2,c))~=0
% %        fan_real=sqrt(real(kw(1,c))^2+real(kw(2,c))^2);
% %        fan_imag=sqrt(imag(kw(1,c))^2+imag(kw(2,c))^2);
% %        real(kw(1,c))/fan_real
%        plot(fr(1,c)/pi,real(kw(1,c)),'g.');
%        hold on;
% %        plot(fr(1,c)/pi,imag(kw(1,c)),'r.');
%        hold on;
% %     else
% %        plot(fr(1,c)/pi,real(kw(1,c)),'b.');
% %        hold on;
% %        plot(fr(1,c)/pi,imag(kw(1,c)),'r.');
% %        hold on;
% end
% %����ͼ��ӱ���
% title(['˥��С���λ��']);
% %�½���ͼ����
% figure;
% for c=1:s(2)-100
%     %���ݹ�һ������
% %     if real(kw(1,c))~=0 && imag(kw(2,c))~=0
%        fan_real=sqrt(real(kw(1,c))^2+real(kw(2,c))^2);
%        fan_imag=sqrt(imag(kw(1,c))^2+imag(kw(2,c))^2);
%        plot(fr(2,c)/pi,real(kw(2,c))/fan_real,'b.');
%        hold on;
%        plot(fr(2,c)/pi,imag(kw(2,c))/fan_imag,'r.');
%        hold on;
% %     else
% %        plot(fr(2,c)/pi,real(kw(2,c)),'b.');
% %        hold on;
% %        plot(fr(2,c)/pi,imag(kw(2,c)),'r.');
% %        hold on;
% end
% %����ͼ��ӱ���
% title(['����С���λ��']);
