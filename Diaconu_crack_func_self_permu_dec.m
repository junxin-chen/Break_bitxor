% Information Sciences 355�C356 (2016) 314�C327
 
function p_dec = Diaconu_crack_func_self_permu_dec(p)
[M,N]=size(p);
%����������
p_permu_dec=zeros(M,N);
for i=1:M
    %���Ȳ��Ϊ���ؾ���
    for j=1:N
        for k=1:8  %Ԥ��Ϊ256���Ҷ�ͼ��
            if(bitand((p(i,j)),2^(8-k)))
                p_bit3(i,(j-1)*8+k)=1;
            else
               p_bit3(i,(j-1)*8+k)=0; 
            end
        end
    end
    %����ѭ����λ��λ��
    N_shifts=sum(p_bit3(i,:));
    %������λ�ķ��򣬲�������λ
    if mod(N_shifts,2)
        p_bit3(i,:)=circshift(p_bit3(i,:),N_shifts,2);
    else 
        p_bit3(i,:)=circshift(p_bit3(i,:),-N_shifts,2);
    end
    %��ϳ��µ�����
    for j=1:N
        for k=1:8
            p_permu_dec(i,j)=p_permu_dec(i,j)+2^(8-k)*p_bit3(i,(j-1)*8+k);
        end
    end
end
p_dec=uint8(p_permu_dec);