clear all;
clc;
syms x;
% parametric input
N=16;%N����Ԫ��

alpha=0.31;%ɨ�費ͬ��alpha����gain/loss����
% G=0.3; %ɨ�費ͬ��G�����ڵ��ɾ���ϵ��
% m=M*0.13;
% k1=0.6;
% k2=0.6;
% w_0_loss = sqrt(2*G*(1.0+alpha*1i)/m);%��ĵ�Ч����Ƶ��
% M_0_loss = M+m*(w_0_loss^2/(w_0_loss^2-x^2));
% w_0_gain = sqrt(2*G*(1.0-alpha*1i)/m);%�����Ч����Ƶ��
% M_0_gain = M+m*(w_0_gain^2/(w_0_gain^2-x^2));
mag=[];
magf=[];
kw=[];
fr=[];
w=0:0.01*pi:pi;

%λ�ƾ���
syms displacement_matrix;
%ǰ���������������С
%ǰ����������
% m=M*0.13;
% k1=0.6;
% k2=0.6;
% w_0_loss = sqrt(2*G*(1.0+alpha*1i)/m);%��ĵ�Ч����Ƶ��
% M_0_loss = M+m*(w_0_loss^2/(w_0_loss^2-x^2));
% w_0_gain = sqrt(2*G*(1.0-alpha*1i)/m);%�����Ч����Ƶ��
% M_0_gain = M+m*(w_0_gain^2/(w_0_gain^2-x^2));
for d = 1:N/2
    M=7.0;
    m = M*0.15;
    G=0.5;
    w_0_loss = sqrt(2*G*(1.0+alpha*1i)/m);%��ĵ�Ч����Ƶ��
    M_0_loss = M+m*(w_0_loss^2/(w_0_loss^2-x^2));
    w_0_gain = sqrt(2*G*(1.0-alpha*1i)/m);%�����Ч����Ƶ��
    M_0_gain = M+m*(w_0_gain^2/(w_0_gain^2-x^2));
    %ż��λ��Ϊ����
    if mod(d,2) == 0
       displacement_matrix(d,d)=-M_0_gain*x^2;
       continue;
    end
    %����λ��Ϊ˥��
    if mod(d,2) ~= 0
       displacement_matrix(d,d)=-M_0_loss*x^2;
       continue;
    end
end

%��������
%���������������Դ�
% m=M*0.25;
% k1=1.25;
% k2=1.25;
% w_0_loss = sqrt(2*G*(1.0+alpha*1i)/m);%��ĵ�Ч����Ƶ��
% M_0_loss = M+m*(w_0_loss^2/(w_0_loss^2-x^2));
% w_0_gain = sqrt(2*G*(1.0-alpha*1i)/m);%�����Ч����Ƶ��
% M_0_gain = M+m*(w_0_gain^2/(w_0_gain^2-x^2));
for d = N/2+1:N
    M=6.0;
    m = M*0.3;
    G=0.9;
    w_0_loss = sqrt(2*G*(1.0+alpha*1i)/m);%��ĵ�Ч����Ƶ��
    M_0_loss = M+m*(w_0_loss^2/(w_0_loss^2-x^2));
    w_0_gain = sqrt(2*G*(1.0-alpha*1i)/m);%�����Ч����Ƶ��
    M_0_gain = M+m*(w_0_gain^2/(w_0_gain^2-x^2));
    %ż��λ��Ϊ����
    if mod(d,2) == 0
       displacement_matrix(d,d)=-M_0_gain*x^2;
       continue;
    end
    %����λ��Ϊ˥��
    if mod(d,2) ~= 0
       displacement_matrix(d,d)=-M_0_loss*x^2;
       continue;
    end
end


% factor_matrix = factor_matrix - displacement_matrix;
%ǰ�벿�־�����С��k����벿���ô��k
k1=2;
k2=2;

%ɨ����λ
for kk=1:length(w)
    %����ϵ������
    %��������λ���ȸ�ֵ
    factor_matrix=zeros(N);%NxN��ȫ0����
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
    
    
% (ԭ���ļ��㷽��)
    %ǰ�벿��
%     for index = 2:N-1
%         %iΪ���������
%         if mod(index,2)~=0
%             factor_matrix(index,index+1)=k1;%k1
%             factor_matrix(index,index)=-(k1+k2);%-(k1+k2)
%             factor_matrix(index,index-1)=k2;%k2
%             continue;
%         end
%          %iΪż�������
%         if mod(index,2)==0
%             factor_matrix(index,index+1)=k2;%k2
%             factor_matrix(index,index)=-(k2+k1);%-(k2+k1)
%             factor_matrix(index,index-1)=k1;%k1
%             continue;
%         end
%     end
    
    %��������������ľ���,�������Ӧ�Ľ�Ƶ��
    factor_matrix= factor_matrix-displacement_matrix;
    func = det(factor_matrix);
    result = solve(func);
    for n = 1:length(result)
%         if double(real(result(n)))>0 && double(abs(imag(result(n))))<=0.01
            plot(w(kk)/pi,real(result(n)),'ro','MarkerFaceColor','g','MarkerSize',2);
            hold on;
%         end
    end
    fprintf('%f\n',kk);
end



% (δ���鲿ʱ�ļ��㷽��)
% eigValue = eig(factor_matrix);%�ó�����ֵ����
% ��������ֵ���ó���Ƶ��
%     for h=1:length(eigValue)
%         omega=roots([M,0.0,eigValue(h)-M*w_0_loss^2-m*w_0_loss^2,0.0,-eigValue(h)*w_0_loss^2]);%����һ������ֵ�������Ӧ�Ľ�Ƶ��
%         for d=1:length(omega)
%              % �����õ��Ľ�Ƶ��
%              kw = [kw;real(omega(d))];
%              magf = [magf;imag(omega(d))];
%               plot(w(kk)/(2*pi),imag(omega(d)),'*','color',[1 0 0]);
%              if real(omega(d))>0 && abs(imag(omega(d)))<=0.01
%                 plot(w(kk)/pi,real(omega(d)),'ro','MarkerFaceColor','g','MarkerSize',2);
%              end
%              hold on;
%         end
%     end



%ɨ����λ
% for i=1:length(w)
%     factor_matrix(1,N)=k1*exp(-i*w(i));
%     factor_matrix(N,1)=k2*exp(i*w(i));
%     eigValue = eig(factor_matrix);%�ó�����ֵ����
%     %��������ֵ���ó���Ƶ��
%     for h=1:length(eigValue)
%         omega=solve((k1+k2-x^2*(M+m*w_0_loss/(w_0_loss^2-x^2)))*(k1+k2-x^2*(M+m*w_0_gain/(w_0_gain^2-x^2)))-2*k1*k2*cos(w(h))-k1^2-k2^2,x);%����һ����λ�������Ӧ�Ľ�Ƶ��
%         for d=1:length(omega)
%              % �����õ��Ľ�Ƶ��
%              kw = [kw;real(omega(d))];
%              magf = [magf;imag(omega(d))];
%              plot(w(h)/(2*pi),real(omega(d)),'*','color',[0 0 1]);
%              plot(w(h)/(2*pi),imag(omega(d)),'*','color',[1 0 0]);
%              hold on;
%         end
%     end
% end
