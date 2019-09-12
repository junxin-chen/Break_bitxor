%ѡ�����Ĺ�������Ч�Բ���

clear all
clc
close all

load Hua_INS_K
encrypt=@(m)(Hua_2016_INS(m,'en',K));
decrypt=@(m)(Hua_2016_INS(m,'de',K));

%����ԭ��
m=imread('lenna32.bmp');
[M_plain,N_plain]=size(m);

cc=encrypt(m);
[M,N]=size(cc);

%%%���ֹ���ģ���¹�����ӵ�н��ܻ���Ȩ�ޣ����������ڽ��ܴ��ƽ��ciphertext
c0=uint8(zeros(M,N));
m0=uint8(decrypt(c0));
sum_del_m=uint8(zeros(M_plain,N_plain));
for i=1:M
    i
    for j=1:N
        for k=1:8 
            %ʵ�ʱ���У�����һ������һ�����صĻָ���������Ҫ�ķѴ�����Դ�������뱾
            %������del_c��Ȼ����del_c����Ϊc0=zeros������del_c=c,���ܵø����Ķ�Ӧ�����ģ�Ȼ����m0��һ������del_m
            %�жϵ�ǰ����λ�Ƿ�λ1����ȷ��lambdaΪ0��1�������0������Ҫִ�����²���
            lambda=bitand(cc(i,j),2^(8-k));
            if lambda==0
                continue;
            else                 
                    current_c=uint8(zeros(M,N));
                    current_c(i,j)=2^(8-k);%��ǰ�������ΰ�ÿһ��bit����
                    current_c=uint8(current_c);
                    current_m=uint8(decrypt(current_c));
                    current_del_m=bitxor(current_m,m0);
                    sum_del_m=bitxor(sum_del_m,current_del_m);                 
            end
        end
    end
end
mm_crack=bitxor(sum_del_m,m0);
dd=double(mm_crack)-double(m);

