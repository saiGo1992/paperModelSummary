%K增大带隙增大
%G增大带隙减小?
%M增大带隙增大
%目标:mass_factor增大，带隙先减小后增大
clear all;
clc;
syms x;
% parametric input

%质量比例
mass_factor=0.15;
alpha=0.31;%扫描不同的alpha，即gain/loss因子
%外弹簧K比例
K_factor=1;
K1=2;
K2=2;
M=7;
m=0.6;
%0.227~0.228之间局域共振带隙简并
G=0.5; %扫描不同的G，即内弹簧劲度系数
mag=[];%按列保存所有频率
%magf=[];%保存所有频率解的虚部

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
                w=0;
                %新建画图窗口
%                 figure;
                for h=1:length(w)
                    %dispersion relation
                    omega=vpasolve((K1+K2-x^2*M_eff_loss)*(K1+K2-x^2*M_eff_gain)-2*K1*K2*cos(w(h))-K1^2-K2^2,x);%代入一个相位，求出对应的角频率
                    index=find(real(double(omega))>=0);%获取实部大于0的解的下标
                    omega=omega(index); %根据下标得到对应的元素
                    omega=unique(omega);%除掉相同元素
%                    temp = ones(size(omega))*w(h);
    %                 mag=[mag,temp];
    %                 magf=[magf,omega];
                    omega=sort(omega);%排序，默认是升序
                    for d=1:length(omega)
                        %画实部
%                         plot(w(h)/pi,real(omega(d)),'m*','LineWidth',2,'MarkerSize',2);
%                         hold on;
                        %画虚部
%                         plot(w(h)/pi,imag(omega(d)),'r*','LineWidth',2,'MarkerSize',2);
%                         hold on;
                        %画等效质量
%                         plot(real(M+m*w_0_loss/(w_0_loss^2-omega(d)^2)),real(omega(d)),'c*','LineWidth',2,'MarkerSize',2);
%                         hold on;
%                         plot(imag(M+m*w_0_loss/(w_0_loss^2-omega(d)^2)),real(omega(d)),'y*','LineWidth',2,'MarkerSize',2);
%                         hold on;
                        %创建矩阵                        
                        matr=[2*K1-omega(d)^2*(M+m*w_0_loss/(w_0_loss^2-omega(d)^2)), -(exp(-1i*w(h))+1)*K1;
                              (exp(1i*w(h))+1)*K1, -(2*K1-omega(d)^2*(M+m*w_0_gain/(w_0_gain^2-omega(d)^2)))];
                        [V,D]=eig(matr);
                        
                        %画绝对值的比值
                        %V=abs(V);
                        %画虚部的比值
%                         V
%                         V=imag(V)
%                         Ratio=V(1,:)./V(2,:);
                        %画实部的比值
%                         V=real(V);
%                         Ratio=V(1,:)./V(2,:);
%                         %开始画
%                         for i=length(Ratio)
%                             if abs(Ratio(i))>200
%                                 continue;
%                             end
%                             if d == 2
%                                 plot(Ratio(i),real(omega(d)),'b*','LineWidth',2,'MarkerSize',2);
%                             end
%                             hold on;
%                         end
                        
                        %直接画本征矢量
                        V
                        V=real(V);
                        for i=length(V(1,:))
%                             if d == 2
                                plot(V(1,i),real(omega(d)),'b*','LineWidth',2,'MarkerSize',2);
                                hold on;
                                plot(V(2,i),real(omega(d)),'r*','LineWidth',2,'MarkerSize',2);
                                hold on;
%                                 plot(V(1,i)/V(2,i),real(omega(d)),'b*','LineWidth',2,'MarkerSize',2);
%                                 hold on;
%                             end
                        end
                        
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



