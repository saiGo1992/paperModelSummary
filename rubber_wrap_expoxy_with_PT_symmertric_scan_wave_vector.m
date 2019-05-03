clear all;
clc;
syms x;
% parametric input
K=1.25;
M=2.0;
m=0.5;
alpha=0.4:0.1:0.5;%扫描不同的alpha，即gain/loss因子
G=0.3:0.1:0.4; %扫描不同的G，即内弹簧劲度系数
solutions = [];%用矩阵按行存储得到解（如频率等）

%扫描不同的alpha，即gain/loss因子
for a = 1:length(alpha)
    %扫描不同的G，即内弹簧劲度系数
    for g=1:length(G)
        %resonant frequency
        w_0_loss = sqrt(2*G(g)*(1+alpha(a)*i)/m);
        w_0_gain = sqrt(2*G(g)*(1-alpha(a)*i)/m);
        %扫描波矢
        w=0:0.01*pi:pi;
        %虚部不为0的频率集合
        jg=0;
        %虚部为0的频率集合
        jf=0;
        %新建画图窗口,画色散曲线
        %figure;
        for h=1:length(w)
            %dispersion relation
            omega=solve(cos(w(h))+1-(2*K-x^2*(M+m*w_0_loss^2/(w_0_loss^2-x^2)))*(2*K-x^2*(M+m*w_0_gain^2/(w_0_gain^2-x^2)))/(2*K^2),x);%代入一个相位，求出对应的角频率
            for d=1:length(omega)
                %记录解(实部和虚部分开),对应的频率以及相对应的条件(alpha和G)
                solutions = [solutions;[alpha(a),G(g),omega(d),w(h)/(2*pi)]];
            end
        end
    end
end

%用记录的数据画色散曲线图
for a = 1:length(alpha)
    for g=1:length(G)
        figure;
        row_colum = size(solutions); %取出解矩阵行数
        for r = 1:row_colum(1)
            if solutions(r,1) == alpha(a) && solutions(r,2) == G(g)
                plot(solutions(r,4),real(solutions(r,3)),'*','color',[0 0 1]);
                plot(solutions(r,4),imag(solutions(r,3)),'*','color',[1 0 0]);
                hold on;
            end
        end
        %给绘图添加标题
        title(['alpha=',num2str(alpha(a)),'　G=',num2str(G(g))]);
    end
end



%看等效质量曲线
for a = 1:length(alpha)
    for g=1:length(G)
        figure;
        row_colum = size(solutions); %取出解矩阵行数
        for r = 1:row_colum(1)
            if solutions(r,1) == alpha(a) && solutions(r,2) == G(g)
                M_eff_loss = M+m*w_0_loss^2/(w_0_loss^2-solutions(r,3)^2);
                plot(solutions(r,4),real(M_eff_loss),'*','color',[0 0 1]);
                plot(solutions(r,4),imag(M_eff_loss),'*','color',[1 0 0]);
                hold on;
            end
        end
        %给绘图添加标题
        title(['alpha=',num2str(alpha(a)),'　G=',num2str(G(g))]);
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