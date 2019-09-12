%Information Sciences 355–356 (2016) 314–327
 
function p_enc = Diaconu_crack_func_self_permu_enc(p)

[M,N]=size(p);
p_permu=zeros(M,N);
for i=1:M
    %首先拆分为比特矩阵
    for j=1:N
        for k=1:8  %预设为256级灰度图像
            if(bitand((p(i,j)),2^(8-k)))
                p_bit(i,(j-1)*8+k)=1;
            else
               p_bit(i,(j-1)*8+k)=0; 
            end
        end
    end
    %计算循环移位的位数
    N_shifts=sum(p_bit(i,:));
    %计算移位的方向，并进行移位
    if mod(N_shifts,2)
        p_bit(i,:)=circshift(p_bit(i,:),-N_shifts,2);
    else 
        p_bit(i,:)=circshift(p_bit(i,:),N_shifts,2);
    end
    %组合成新的像素
    for j=1:N
        for k=1:8
            p_permu(i,j)=p_permu(i,j)+2^(8-k)*p_bit(i,(j-1)*8+k);
        end
    end
end

p_enc=uint8(p_permu);


