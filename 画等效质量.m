clear all;
clc;
syms x;
% parametric input
K1=1.25;
K2=1.25;
M=2;
m=0.5;
alpha=0.3;%ɨ�費ͬ��alpha����gain/loss����
G=0.3; %ɨ�費ͬ��G�����ڵ��ɾ���ϵ��
mag=[];
magf=[];
kw=[];
fr=[];

w_0_loss = sqrt(2*G*(1.0+alpha*1i)/m);
w_0_gain = sqrt(2*G*(1.0-alpha*1i)/m);
M_eff_loss = M+m*w_0_loss/(w_0_loss^2-x^2);
M_eff_gain = M+m*w_0_gain/(w_0_gain^2-x^2);

x=0:0.001:5;
figure;
plot(x,real(M+m*w_0_loss./(w_0_loss^2-x.^2)),'b*');
hold on;
plot(x,imag(M+m*w_0_loss./(w_0_loss^2-x.^2)),'c*');
%plot(x,real(M+m*(1.0-alpha*1i)*w_0_gain./(w_0_gain^2-x.^2)),'r*');
% figure;
% plot(x,real(M+m*w_0_gain./(w_0_gain^2-x.^2)),'b*');
% hold on;
% plot(x,imag(M+m*w_0_gain./(w_0_gain^2-x.^2)),'r*');
%ɨ�費ͬ��alpha����gain/loss����
% for a = 1:length(alpha)
%     %ɨ�費ͬ��G�����ڵ��ɾ���ϵ��
%     for g=1:length(G)
%         %resonant frequency
%         w_0_loss = sqrt(2*G(g)*(1.0+alpha(a)*i)/m);
%         w_0_gain = sqrt(2*G(g)*(1.0-alpha(a)*i)/m);
%         %��Ч����
%         M_eff_loss = M+m*w_0_loss/(w_0_loss^2-x^2);
%         M_eff_gain = M+m*w_0_gain/(w_0_gain^2-x^2);
%         %ɨ�貨ʸ
%         w=0:0.005*pi:pi;
%         %�鲿��Ϊ0��Ƶ�ʼ���
%         jg=0;
%         %�鲿Ϊ0��Ƶ�ʼ���
%         jf=0;
%         %�½���ͼ����
% %         figure;
%         for h=1:length(w)
%             %dispersion relation
%             omega=solve((K1+K2-x^2*M_eff_loss)*(K1+K2-x^2*M_eff_gain)-2*K1*K2*cos(w(h))-K1^2-K2^2,x);%����һ����λ�������Ӧ�Ľ�Ƶ��
%             for d=1:length(omega)
%                 %��¼����0�Ľ�(ʵ�����鲿�ֿ���)�Ͷ�Ӧ��Ƶ��
% %                 if double(real(omega(d)))>=0 
% %                    kw = [kw;real(omega(d))];
% %                    magf = [magf;imag(omega(d))];
% %                    plot(w(h)/(pi),real(omega(d)),'b*-');
% %                    %if imag(omega(d))>0
% %                        plot(w(h)/(pi),imag(omega(d)),'r*-');
% %                    %end                    
% %                end
%                 
%                 hold on;
%             end
%         end
%         %����ͼ��ӱ���
%         title(['alpha=',num2str(alpha(a)),'��G=',num2str(G(g))]);
%     end
% end

%��Ч��������(������)

% plot(omega(d),imag(M+m*w_0_loss/(w_0_loss^2-omega(d)^2)),'r*');



