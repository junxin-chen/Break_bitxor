% the basic cipher (decryption part) used for illustration in the paper
% logistic map is introduced for generating the encryption elements, i.e., permutation vector and
% substitution mask
% the encryption elements are different in distinct interations
% the encryption core is repeated 10 times

function d5 = basic_dec(d1)

% the logistic map, and its initial parameters

[M,N]=size(d1);
miu0=3.99976;
x00=0.92837471;
count=10;% iteration 10 times

% initial iteration 300 times for avoiding non-chaotic points 
for j=1:300
    x00=miu0*x00*(1-x00);
end


% decryption core in a iteration
for k=1:count
    d2=reshape(d1.',1,M*N);
    
    miu=miu0+(count+1-k)*10^(-5);
    x0=x00;
    for i=1:M
    for j=1:N        
        x0=miu*x0*(1-x0);
        key(i,j)=mod(floor(x0*10^14),256);        
    end
    end
    key=uint8(key);
    key2=reshape(key.',1,M*N);   
    
    % in the decryption, subsitution and permutation should be decoded in a
    % reverse order
    % decode the subsitution firstly
    for j=1:M*N
        tmp=d2(j);
        d3(j)=bitxor(d2(j),key2(j));
    end
    d4=reshape(d3,M,N).';    
    
    % decode the permutation part
    d5=dec_bitpermu(d4,miu,x0);
    d1=d5; 
end

% decryption core of the permutation

function p_permu = dec_bitpermu(p,miu,x0)

[M,N]=size(p);

p_bit=zeros(M,8*N);
p_bit_permu=zeros(M,8*N);
p_permu=zeros(M,N);

for i=1:M
    for j=1:N
        for k=1:8  
            if(bitand((p(i,j)),2^(8-k)))
                p_bit(i,(j-1)*8+k)=1;
            else
               p_bit(i,(j-1)*8+k)=0; 
            end
        end
    end
end


% initial iteration 300 times for avoiding non-chaotic points 
for i=1:300
    x0=miu*x0*(1-x0);
end

% produce the row-shuffling vector
for i=1:M
    x0=miu*x0*(1-x0);
    hang(1,i)=x0;
end
[hang_sortint,hang_index]=sort(hang);

% produce the column-shuffling vector
for i=1:8*M
    x0=miu*x0*(1-x0);
    lie(1,i)=x0;
end
[lie_sorting,lie_index]=sort(lie);

% decryption using the column-shuffling vector firstly
for i=1:8*N
    p_bit_permu(:,lie_index(i))=p_bit(:,i);
end

% decryption using the row-shuffling vector secondly
for i=1:M
    p_bit_permu2(hang_index(i),:)=p_bit_permu(i,:);
end

% convert to the gray-level
for i=1:M
    for j=1:N
        for k=1:8
            p_permu(i,j)=p_permu(i,j)+2^(8-k)*p_bit_permu2(i,(j-1)*8+k);
        end
    end
end
p_permu=uint8(p_permu);

