%K�����϶����
%G�����϶��С?
%M�����϶����
%Ŀ��:mass_factor���󣬴�϶�ȼ�С������
clear all;
clc;
syms x;
% parametric input

%��������
mass_factor=0.15;
alpha=0.31;%ɨ�費ͬ��alpha����gain/loss����
%�ⵯ��K����
K_factor=1;
K1=2;
K2=2;
M=7;
m=0.6;
%0.227~0.228֮��������϶��
G=0.5; %ɨ�費ͬ��G�����ڵ��ɾ���ϵ��
mag=[];%���б�������Ƶ��
%magf=[];%��������Ƶ�ʽ���鲿

%ɨ�費ͬ��K
for e=1:length(K_factor)
    K2=K1*K_factor(e);
    %ɨ�費ͬ����������
    for c=1:length(mass_factor)
        m=M*mass_factor(c);
        %ɨ�費ͬ��alpha����gain/loss����
        for a = 1:length(alpha)
            %ɨ�費ͬ��G�����ڵ��ɾ���ϵ��
            for g=1:length(G)
                %resonant frequency
                w_0_loss = sqrt(2*G(g)*(1.0+alpha(a)*1i)/m);
                w_0_gain = sqrt(2*G(g)*(1.0-alpha(a)*1i)/m);
                %��Ч����
                M_eff_loss = M+m*w_0_loss/(w_0_loss^2-x^2);
                M_eff_gain = M+m*w_0_gain/(w_0_gain^2-x^2);
                %ɨ�貨ʸ
                w=0;
                %�½���ͼ����
%                 figure;
                for h=1:length(w)
                    %dispersion relation
                    omega=vpasolve((K1+K2-x^2*M_eff_loss)*(K1+K2-x^2*M_eff_gain)-2*K1*K2*cos(w(h))-K1^2-K2^2,x);%����һ����λ�������Ӧ�Ľ�Ƶ��
                    index=find(real(double(omega))>=0);%��ȡʵ������0�Ľ���±�
                    omega=omega(index); %�����±�õ���Ӧ��Ԫ��
                    omega=unique(omega);%������ͬԪ��
%                    temp = ones(size(omega))*w(h);
    %                 mag=[mag,temp];
    %                 magf=[magf,omega];
                    omega=sort(omega);%����Ĭ��������
                    for d=1:length(omega)
                        %��ʵ��
%                         plot(w(h)/pi,real(omega(d)),'m*','LineWidth',2,'MarkerSize',2);
%                         hold on;
                        %���鲿
%                         plot(w(h)/pi,imag(omega(d)),'r*','LineWidth',2,'MarkerSize',2);
%                         hold on;
                        %����Ч����
%                         plot(real(M+m*w_0_loss/(w_0_loss^2-omega(d)^2)),real(omega(d)),'c*','LineWidth',2,'MarkerSize',2);
%                         hold on;
%                         plot(imag(M+m*w_0_loss/(w_0_loss^2-omega(d)^2)),real(omega(d)),'y*','LineWidth',2,'MarkerSize',2);
%                         hold on;
                        %��������                        
                        matr=[2*K1-omega(d)^2*(M+m*w_0_loss/(w_0_loss^2-omega(d)^2)), -(exp(-1i*w(h))+1)*K1;
                              (exp(1i*w(h))+1)*K1, -(2*K1-omega(d)^2*(M+m*w_0_gain/(w_0_gain^2-omega(d)^2)))];
                        [V,D]=eig(matr);
                        
                        %������ֵ�ı�ֵ
                        %V=abs(V);
                        %���鲿�ı�ֵ
%                         V
%                         V=imag(V)
%                         Ratio=V(1,:)./V(2,:);
                        %��ʵ���ı�ֵ
%                         V=real(V);
%                         Ratio=V(1,:)./V(2,:);
%                         %��ʼ��
%                         for i=length(Ratio)
%                             if abs(Ratio(i))>200
%                                 continue;
%                             end
%                             if d == 2
%                                 plot(Ratio(i),real(omega(d)),'b*','LineWidth',2,'MarkerSize',2);
%                             end
%                             hold on;
%                         end
                        
                        %ֱ�ӻ�����ʸ��
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
                %����ͼ��ӱ���
                title(['alpha=',num2str(alpha(a)),'  G=',num2str(G(g)),'  KFactor=',num2str(K1),'  massFactor=',num2str(mass_factor(c)),' M=',num2str(M)]);
    %             title(['alpha=',num2str(alpha(a)),'  massFactor=',num2str(mass_factor(c))])
    %             title(['alpha=',num2str(alpha(a)),'  K=',num2str(K1)]);
            end
        end
    end
end



