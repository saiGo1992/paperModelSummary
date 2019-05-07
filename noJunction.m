clear all;
clc;
syms x;
% parametric input
N=16;%N个超元胞

alpha=0.3;%扫描不同的alpha，即gain/loss因子
% G=0.3; %扫描不同的G，即内弹簧劲度系数
% m=M*0.13;
% k1=0.6;
% k2=0.6;
% w_0_loss = sqrt(2*G*(1.0+alpha*1i)/m);%损耗等效本征频率
% M_0_loss = M+m*(w_0_loss^2/(w_0_loss^2-x^2));
% w_0_gain = sqrt(2*G*(1.0-alpha*1i)/m);%增益等效本征频率
% M_0_gain = M+m*(w_0_gain^2/(w_0_gain^2-x^2));
mag=[];
magf=[];
kw=[];
fr=[];
w=0:0.01*pi:pi;

%位移矩阵
syms displacement_matrix;
%前半个矩阵，质量比稍小
%前半个矩阵参数
% m=M*0.13;
% k1=0.6;
% k2=0.6;
% w_0_loss = sqrt(2*G*(1.0+alpha*1i)/m);%损耗等效本征频率
% M_0_loss = M+m*(w_0_loss^2/(w_0_loss^2-x^2));
% w_0_gain = sqrt(2*G*(1.0-alpha*1i)/m);%增益等效本征频率
% M_0_gain = M+m*(w_0_gain^2/(w_0_gain^2-x^2));
for d = 1:N/2
    M=5.0;
    m = M*0.15;
    G=0.3;
    w_0_loss = sqrt(2*G*(1.0+alpha*1i)/m);%损耗等效本征频率
    M_0_loss = M+m*(w_0_loss^2/(w_0_loss^2-x^2));
    w_0_gain = sqrt(2*G*(1.0-alpha*1i)/m);%增益等效本征频率
    M_0_gain = M+m*(w_0_gain^2/(w_0_gain^2-x^2));
    %偶数位置为增益
    if mod(d,2) == 0
       displacement_matrix(d,d)=-M_0_gain*x^2;
       continue;
    end
    %奇数位置为衰减
    if mod(d,2) ~= 0
       displacement_matrix(d,d)=-M_0_loss*x^2;
       continue;
    end
end

%后半个矩阵
%后半个矩阵，质量比稍大
% m=M*0.25;
% k1=1.25;
% k2=1.25;
% w_0_loss = sqrt(2*G*(1.0+alpha*1i)/m);%损耗等效本征频率
% M_0_loss = M+m*(w_0_loss^2/(w_0_loss^2-x^2));
% w_0_gain = sqrt(2*G*(1.0-alpha*1i)/m);%增益等效本征频率
% M_0_gain = M+m*(w_0_gain^2/(w_0_gain^2-x^2));
for d = N/2+1:N
    M=4.0;
    m = M*0.3;
    G=0.55;
    w_0_loss = sqrt(2*G*(1.0+alpha*1i)/m);%损耗等效本征频率
    M_0_loss = M+m*(w_0_loss^2/(w_0_loss^2-x^2));
    w_0_gain = sqrt(2*G*(1.0-alpha*1i)/m);%增益等效本征频率
    M_0_gain = M+m*(w_0_gain^2/(w_0_gain^2-x^2));
    %偶数位置为增益
    if mod(d,2) == 0
       displacement_matrix(d,d)=-M_0_gain*x^2;
       continue;
    end
    %奇数位置为衰减
    if mod(d,2) ~= 0
       displacement_matrix(d,d)=-M_0_loss*x^2;
       continue;
    end
end


% factor_matrix = factor_matrix - displacement_matrix;
%前半部分矩阵用小的k，后半部分用大的k
k1=1.25;
k2=1.25;

%扫描相位
for kk=1:length(w)
    %构建系数矩阵
    %几个特殊位置先赋值
    factor_matrix=zeros(N);%NxN的全0矩阵
    factor_matrix(1,1)=-(k1+k2);
    factor_matrix(1,2)=k1;
    factor_matrix(1,N)=k2*exp(-1i*w(kk));
    factor_matrix(N,N)=-2*k2;
    factor_matrix(N,N-1)=k2;
    factor_matrix(N,1)=k2*exp(1i*w(kk));
    factor_matrix(N/2+1,N/2+1)=-(k1+k2);
    factor_matrix(N/2+1,N/2+2)=k2;
    factor_matrix(N/2+1,N/2)=k1;
    for index = 2:N/2
        factor_matrix(index,index)=-2*k1;
        factor_matrix(index,index-1)=k1;
        factor_matrix(index,index+1)=k1;
    end
    for index=N/2+2:N-1
        factor_matrix(index,index)=-2*k2;
        factor_matrix(index,index-1)=k2;
        factor_matrix(index,index+1)=k2;
    end
    
    
% (原来的计算方法)
    %前半部分
%     for index = 2:N-1
%         %i为奇数的情况
%         if mod(index,2)~=0
%             factor_matrix(index,index+1)=k1;%k1
%             factor_matrix(index,index)=-(k1+k2);%-(k1+k2)
%             factor_matrix(index,index-1)=k2;%k2
%             continue;
%         end
%          %i为偶数的情况
%         if mod(index,2)==0
%             factor_matrix(index,index+1)=k2;%k2
%             factor_matrix(index,index)=-(k2+k1);%-(k2+k1)
%             factor_matrix(index,index-1)=k1;%k1
%             continue;
%         end
%     end
    
    %建立起最后完整的矩阵,并算出对应的角频率
    factor_matrix= factor_matrix-displacement_matrix;
    func = det(factor_matrix);
    result = solve(func);
    for n = 1:length(result)
        if double(real(result(n)))>0
%         if double(real(result(n)))>0 && double(abs(imag(result(n))))<=0.01
            plot(w(kk)/pi,real(result(n)),'ro','MarkerFaceColor','g','MarkerSize',2);
            hold on;
        end
    end
    fprintf('%f\n',kk);
end



% (未加虚部时的计算方法)
% eigValue = eig(factor_matrix);%得出本征值向量
% 遍历本征值，得出角频率
%     for h=1:length(eigValue)
%         omega=roots([M,0.0,eigValue(h)-M*w_0_loss^2-m*w_0_loss^2,0.0,-eigValue(h)*w_0_loss^2]);%代入一个本征值，求出对应的角频率
%         for d=1:length(omega)
%              % 遍历得到的角频率
%              kw = [kw;real(omega(d))];
%              magf = [magf;imag(omega(d))];
%               plot(w(kk)/(2*pi),imag(omega(d)),'*','color',[1 0 0]);
%              if real(omega(d))>0 && abs(imag(omega(d)))<=0.01
%                 plot(w(kk)/pi,real(omega(d)),'ro','MarkerFaceColor','g','MarkerSize',2);
%              end
%              hold on;
%         end
%     end



%扫描相位
% for i=1:length(w)
%     factor_matrix(1,N)=k1*exp(-i*w(i));
%     factor_matrix(N,1)=k2*exp(i*w(i));
%     eigValue = eig(factor_matrix);%得出本征值向量
%     %遍历本征值，得出角频率
%     for h=1:length(eigValue)
%         omega=solve((k1+k2-x^2*(M+m*w_0_loss/(w_0_loss^2-x^2)))*(k1+k2-x^2*(M+m*w_0_gain/(w_0_gain^2-x^2)))-2*k1*k2*cos(w(h))-k1^2-k2^2,x);%代入一个相位，求出对应的角频率
%         for d=1:length(omega)
%              % 遍历得到的角频率
%              kw = [kw;real(omega(d))];
%              magf = [magf;imag(omega(d))];
%              plot(w(h)/(2*pi),real(omega(d)),'*','color',[0 0 1]);
%              plot(w(h)/(2*pi),imag(omega(d)),'*','color',[1 0 0]);
%              hold on;
%         end
%     end
% end
