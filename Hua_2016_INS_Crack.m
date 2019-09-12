%选择密文攻击的有效性测试

clear all
clc
close all

load Hua_INS_K
encrypt=@(m)(Hua_2016_INS(m,'en',K));
decrypt=@(m)(Hua_2016_INS(m,'de',K));

%加密原文
m=imread('lenna32.bmp');
[M_plain,N_plain]=size(m);

cc=encrypt(m);
[M,N]=size(cc);

%%%这种攻击模型下攻击者拥有解密机的权限，但不能用于解密待破解的ciphertext
c0=uint8(zeros(M,N));
m0=uint8(decrypt(c0));
sum_del_m=uint8(zeros(M_plain,N_plain));
for i=1:M
    i
    for j=1:N
        for k=1:8 
            %实际编程中，可以一个比特一个比特的恢复，而不需要耗费大量资源建立密码本
            %先生成del_c，然后由del_c，因为c0=zeros，所以del_c=c,解密得该密文对应的明文，然后由m0进一步生成del_m
            %判断当前比特位是否位1，以确定lambda为0或1，如果是0，则不需要执行以下操作
            lambda=bitand(cc(i,j),2^(8-k));
            if lambda==0
                continue;
            else                 
                    current_c=uint8(zeros(M,N));
                    current_c(i,j)=2^(8-k);%有前到后，依次把每一个bit置零
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

