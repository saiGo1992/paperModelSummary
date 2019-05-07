%K增大带隙增大
%G增大带隙减小?
%M增大带隙增大
%目标:mass_factor增大，带隙先减小后增大
clear all;
clc;
syms x;
% parametric input

%质量比例
mass_factor=0.19;
alpha=0.3;%扫描不同的alpha，即gain/loss因子
%外弹簧K比例
K_factor=1;
K1=1.25;
K2=1.25;
M=4;
m=2;
%0.227~0.228之间局域共振带隙简并
G=0.3; %扫描不同的G，即内弹簧劲度系数
mag=[];%保存所有频率解的实部
magf=[];%保存所有频率解的虚部

%扫描不同的K
for e=1:length(K_factor)
    K2=K1*K_factor(e);
    %扫描不同的质量比例
    for c=1:length(mass_factor)
        m=M*mass_factor(c);
        %扫描不同的alpha，即gain/loss因子
        for a = 1:length(alpha)
            %扫描不同的G，即内弹簧劲度系数
            for g=1:length(G)
                %resonant frequency
                w_0_loss = sqrt(2*G(g)*(1.0+alpha(a)*1i)/m);
                w_0_gain = sqrt(2*G(g)*(1.0-alpha(a)*1i)/m);
                %等效质量
                M_eff_loss = M+m*w_0_loss/(w_0_loss^2-x^2);
                M_eff_gain = M+m*w_0_gain/(w_0_gain^2-x^2);
                %扫描波矢
                w=linspace(-pi,pi);
                %新建画图窗口
%                 figure;
                for h=1:length(w)
                    %dispersion relation
                    omega=vpasolve((K1+K2-x^2*M_eff_loss)*(K1+K2-x^2*M_eff_gain)-2*K1*K2*cos(w(h))-K1^2-K2^2,x);%代入一个相位，求出对应的角频率
                    index=find(real(double(omega))>=0);%获取实部大于0的解的下标
                    omega=omega(index); %根据下标得到对应的元素
                    temp = ones(size(omega))*w(h);
    %                 mag=[mag,temp];
    %                 magf=[magf,omega];
                    for d=1:length(omega)
                        %画实部
                        plot(w(h)/pi,real(omega(d)),'b*','LineWidth',2,'MarkerSize',2);
                        hold on;
                        %画虚部
%                         plot(w(h)/pi,imag(omega(d)),'r*','LineWidth',2,'MarkerSize',2);
%                         hold on;

                        %画等效质量
%                         plot(real(M+m*w_0_loss/(w_0_loss^2-omega(d)^2)),real(omega(d)),'c*','LineWidth',2,'MarkerSize',2);
%                         hold on;
%                         plot(imag(M+m*w_0_loss/(w_0_loss^2-omega(d)^2)),real(omega(d)),'y*','LineWidth',2,'MarkerSize',2);
%                         hold on;
                          %解出对应的位移

                    end
                end
                %给绘图添加标题
                title(['alpha=',num2str(alpha(a)),'  G=',num2str(G(g)),'  KFactor=',num2str(K1),'  massFactor=',num2str(mass_factor(c)),' M=',num2str(M)]);
    %             title(['alpha=',num2str(alpha(a)),'  massFactor=',num2str(mass_factor(c))])
    %             title(['alpha=',num2str(alpha(a)),'  K=',num2str(K1)]);
            end
        end
    end
end
% for i=1:length(magf(:,1))
%     scatter(mag(i,:),real(magf(i,:)),'bo');
%     hold on;
%     scatter(mag(i,:),imag(magf(i,:)),'ro');
%     hold on;
% end





%    q = acos((2*K-w(h)^2*(M+m*w_0_loss^2/(w_0_loss^2-w(h)^2)))*(2*K-w(h)^2*(M+m*w_0_gain^2/(w_0_gain^2-w(h)^2)))/(2*K^2)-1)/(2*d);
%     if imag(q) ~= 0 
%         jg=jg+1;
%         mag(jg)=imag(q);
%         magf(jg)=w(h);
%     else
%         jf=jf+1;
%         kw(jf)=q;
%         fr(jf)=w(h);
%     end

%全部波矢的实部和虚部分别取出并画图
% jg=jg+1;
% mag(jg)=imag(q);
% magf(jg)=w(h);
% jf=jf+1;
% kw(jf)=real(q);
% fr(jf)=w(h);

%画图，色散曲线
% scatter(kw,fr,'*');
% hold on;
% scatter(mag,magf,'*','r');