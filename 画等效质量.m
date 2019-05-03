clear all;
clc;
syms x;
% parametric input
K1=1.25;
K2=1.25;
M=2;
m=0.5;
alpha=0.3;%扫描不同的alpha，即gain/loss因子
G=0.3; %扫描不同的G，即内弹簧劲度系数
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
%扫描不同的alpha，即gain/loss因子
% for a = 1:length(alpha)
%     %扫描不同的G，即内弹簧劲度系数
%     for g=1:length(G)
%         %resonant frequency
%         w_0_loss = sqrt(2*G(g)*(1.0+alpha(a)*i)/m);
%         w_0_gain = sqrt(2*G(g)*(1.0-alpha(a)*i)/m);
%         %等效质量
%         M_eff_loss = M+m*w_0_loss/(w_0_loss^2-x^2);
%         M_eff_gain = M+m*w_0_gain/(w_0_gain^2-x^2);
%         %扫描波矢
%         w=0:0.005*pi:pi;
%         %虚部不为0的频率集合
%         jg=0;
%         %虚部为0的频率集合
%         jf=0;
%         %新建画图窗口
% %         figure;
%         for h=1:length(w)
%             %dispersion relation
%             omega=solve((K1+K2-x^2*M_eff_loss)*(K1+K2-x^2*M_eff_gain)-2*K1*K2*cos(w(h))-K1^2-K2^2,x);%代入一个相位，求出对应的角频率
%             for d=1:length(omega)
%                 %记录大于0的解(实部和虚部分开画)和对应的频率
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
%         %给绘图添加标题
%         title(['alpha=',num2str(alpha(a)),'　G=',num2str(G(g))]);
%     end
% end

%等效质量曲线(有两个)

% plot(omega(d),imag(M+m*w_0_loss/(w_0_loss^2-omega(d)^2)),'r*');



