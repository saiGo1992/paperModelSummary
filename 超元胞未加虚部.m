clear all;
clc;
syms x;
% parametric input
N=40;%N����Ԫ��
factor_matrix=zeros(N);%NxN��ȫ0����
k1=1.0;
k2=1.5;
M=2.0;
m=0.5;
alpha=0;%ɨ�費ͬ��alpha����gain/loss����
G=0.25; %ɨ�費ͬ��G�����ڵ��ɾ���ϵ��
w_0_loss = sqrt(2*G*(1.0+alpha*1i)/m);%��ĵ�Ч����Ƶ��
M_0_loss = M+m*(w_0_loss^2/(w_0_loss^2-x^2));
w_0_gain = sqrt(2*G*(1.0-alpha*1i)/m);%�����Ч����Ƶ��
M_0_gain = M+m*(w_0_loss^2/(w_0_loss^2-x^2));
mag=[];
magf=[];
kw=[];
fr=[];
w=0:0.005*pi:pi;
%����ϵ������
%��������λ���ȸ�ֵ
factor_matrix(1,1)=-2*k1;%-2*k1;
factor_matrix(N/2,N/2)=-2*k1;%-2*k1;
factor_matrix(N,N)=-2*k2;%-2*k2;
factor_matrix(1,2)=k1;k1;
%λ�ƾ���
syms displacement_matrix;
for d = 1:N
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
%ɨ����λ
for kk=1:length(w)
    factor_matrix(1,N)=k1*exp(-1i*w(kk));
    factor_matrix(N,1)=k2*exp(1i*w(kk));
    %ǰ�벿��
    for index = 2:N/2
        %i=N/2���������
        if index==N/2
            factor_matrix(index,index+1)=k1;%k1
            factor_matrix(index,index-1)=k1;%k1
            continue;
        end
        %iΪ���������
        if mod(index,2)~=0
            factor_matrix(index,index+1)=k1;%k1
            factor_matrix(index,index)=-(k1+k2);%-(k1+k2)
            factor_matrix(index,index-1)=k2;%k2
            continue;
        end
         %iΪż�������
        if mod(index,2)==0
            factor_matrix(index,index+1)=k2;%k2
            factor_matrix(index,index)=-(k2+k1);%-(k2+k1)
            factor_matrix(index,index-1)=k1;%k1
            continue;
        end
    end
    %��벿��
    for index = N/2+1:N
        %i=N���������
        if index==N
            factor_matrix(index,index-1)=k2;%k2
            continue;
        end
        %iΪ���������
        if mod(index,2)~=0
            factor_matrix(index,index+1)=k2;%k2
            factor_matrix(index,index)=-(k2+k1);%-(k2+k1)
            factor_matrix(index,index-1)=k1;%k1
            continue;
        end
         %iΪż�������
        if mod(index,2)==0
            factor_matrix(index,index+1)=k1;%k1
            factor_matrix(index,index)=-(k1+k2);%-(k1+k2)
            factor_matrix(index,index-1)=k2;%k2
            continue;
        end
    end
    % (δ���鲿ʱ�ļ��㷽��)
    eigValue = eig(factor_matrix);%�ó�����ֵ����
    % ��������ֵ���ó���Ƶ��
    for h=1:length(eigValue)
        omega=roots([M,0.0,eigValue(h)-M*w_0_loss^2-m*w_0_loss^2,0.0,-eigValue(h)*w_0_loss^2]);%����һ������ֵ�������Ӧ�Ľ�Ƶ��
        for d=1:length(omega)
             % �����õ��Ľ�Ƶ��
             kw = [kw;real(omega(d))];
             magf = [magf;imag(omega(d))];
%               plot(w(kk)/(2*pi),imag(omega(d)),'*','color',[1 0 0]);
             if real(omega(d))>0 && abs(imag(omega(d)))<=0.01
                %����ڽ������򻭺���
                if real(omega(d)>0.7808 && real(omega(d))<0.866)
                    plot(w(kk)/pi,real(omega(d)),'r.','LineWidth',0.75);
                elseif real(omega(d))>0.9312 && real(omega(d))<1.118
                    plot(w(kk)/pi,real(omega(d)),'r.','LineWidth',0.75);
                elseif real(omega(d))>1.281 && real(omega(d))<1.414
                    plot(w(kk)/pi,real(omega(d)),'r.','LineWidth',0.75);
                else
                    %���໭����
                    plot(w(kk)/pi,real(omega(d)),'b.','LineWidth',0.75);
                end
             end
             hold on;
        end
    end
end




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
