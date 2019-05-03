clear all;
clc;
syms x;
% parametric input
K=1.25;
M=2.0;
m=0.5;
alpha=0.4:0.1:0.5;%ɨ�費ͬ��alpha����gain/loss����
G=0.3:0.1:0.4; %ɨ�費ͬ��G�����ڵ��ɾ���ϵ��
solutions = [];%�þ����д洢�õ��⣨��Ƶ�ʵȣ�

%ɨ�費ͬ��alpha����gain/loss����
for a = 1:length(alpha)
    %ɨ�費ͬ��G�����ڵ��ɾ���ϵ��
    for g=1:length(G)
        %resonant frequency
        w_0_loss = sqrt(2*G(g)*(1+alpha(a)*i)/m);
        w_0_gain = sqrt(2*G(g)*(1-alpha(a)*i)/m);
        %ɨ�貨ʸ
        w=0:0.01*pi:pi;
        %�鲿��Ϊ0��Ƶ�ʼ���
        jg=0;
        %�鲿Ϊ0��Ƶ�ʼ���
        jf=0;
        %�½���ͼ����,��ɫɢ����
        %figure;
        for h=1:length(w)
            %dispersion relation
            omega=solve(cos(w(h))+1-(2*K-x^2*(M+m*w_0_loss^2/(w_0_loss^2-x^2)))*(2*K-x^2*(M+m*w_0_gain^2/(w_0_gain^2-x^2)))/(2*K^2),x);%����һ����λ�������Ӧ�Ľ�Ƶ��
            for d=1:length(omega)
                %��¼��(ʵ�����鲿�ֿ�),��Ӧ��Ƶ���Լ����Ӧ������(alpha��G)
                solutions = [solutions;[alpha(a),G(g),omega(d),w(h)/(2*pi)]];
            end
        end
    end
end

%�ü�¼�����ݻ�ɫɢ����ͼ
for a = 1:length(alpha)
    for g=1:length(G)
        figure;
        row_colum = size(solutions); %ȡ�����������
        for r = 1:row_colum(1)
            if solutions(r,1) == alpha(a) && solutions(r,2) == G(g)
                plot(solutions(r,4),real(solutions(r,3)),'*','color',[0 0 1]);
                plot(solutions(r,4),imag(solutions(r,3)),'*','color',[1 0 0]);
                hold on;
            end
        end
        %����ͼ��ӱ���
        title(['alpha=',num2str(alpha(a)),'��G=',num2str(G(g))]);
    end
end



%����Ч��������
for a = 1:length(alpha)
    for g=1:length(G)
        figure;
        row_colum = size(solutions); %ȡ�����������
        for r = 1:row_colum(1)
            if solutions(r,1) == alpha(a) && solutions(r,2) == G(g)
                M_eff_loss = M+m*w_0_loss^2/(w_0_loss^2-solutions(r,3)^2);
                plot(solutions(r,4),real(M_eff_loss),'*','color',[0 0 1]);
                plot(solutions(r,4),imag(M_eff_loss),'*','color',[1 0 0]);
                hold on;
            end
        end
        %����ͼ��ӱ���
        title(['alpha=',num2str(alpha(a)),'��G=',num2str(G(g))]);
    end
end

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

%ȫ����ʸ��ʵ�����鲿�ֱ�ȡ������ͼ
% jg=jg+1;
% mag(jg)=imag(q);
% magf(jg)=w(h);
% jf=jf+1;
% kw(jf)=real(q);
% fr(jf)=w(h);

%��ͼ��ɫɢ����
% scatter(kw,fr,'*');
% hold on;
% scatter(mag,magf,'*','r');