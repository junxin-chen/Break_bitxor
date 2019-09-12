% the basic cipher (encryption part) used for illustration in the paper
% logistic map is introduced for generating the encryption elements, i.e., permutation vector and
% substitution mask
% the encryption elements are different in distinct interations
% the encryption core is repeated 10 times

function p5 = basic_enc(p)

[M,N]=size(p);
% the logistic map, and its initial parameters
miu0=3.99976;
x00=0.92837471;
count=10; % iteration 10 times

% initial iteration 300 times for avoiding non-chaotic points 
for j=1:300
    x00=miu0*x00*(1-x00);
end

% encryption core in a iteration
for k=1:count
    
    % new parameter in each iteration so as to produce different encryption
    % elements
    miu=miu0+k*10^(-5);
    x0=x00;    
    for i=1:M
    for j=1:N
        x0=miu*x0*(1-x0);
        key(i,j)=mod(floor(x0*10^14),256);   % produce the substituion masks
    end
    end
    key=uint8(key);
    key2=reshape(key.',1,M*N);   
        
   % encryption: bit-level permutation
    p2=enc_bitpermu(p,miu,x0);

    % encryption: substitution using bit-level XOR
    p3=reshape(p2.',1,M*N);
    for j=1:M*N
        p4(j)=bitxor(p3(j),key2(j));
    end
    p5=reshape(p4,M,N).';
    p=p5;
end


% The following is the permutation function used above 
function p_permu = enc_bitpermu(p,miu,x0)

[M,N]=size(p);

p_bit=zeros(M,8*N);
p_bit_permu=zeros(M,8*N);
p_bit_permu2=zeros(M,8*N);
p_permu=zeros(M,N);

%split the original image into binary matrix
for i=1:M
    for j=1:N
        for k=1:8  % assume the gray level is 256
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

% produce the row-shuffling matrix
for i=1:M
    x0=miu*x0*(1-x0);
    hang(1,i)=x0;
end
[hang_sortint,hang_index]=sort(hang);

% produce the column-shuffling matrix
for i=1:8*M
    x0=miu*x0*(1-x0);
    lie(1,i)=x0;
end
[lie_sorting,lie_index]=sort(lie);

% sort the binary matrix using the row-shuffling matrix
for i=1:M
    p_bit_permu(i,:)=p_bit(hang_index(i),:);
end
% sort the binary matrix using the column-shuffling matrix
for i=1:8*N
    p_bit_permu2(:,i)=p_bit_permu(:,lie_index(i));
end

% convert the binary matrix into a gray-level image
for i=1:M
    for j=1:N
        for k=1:8
            p_permu(i,j)=p_permu(i,j)+2^(8-k)*p_bit_permu2(i,(j-1)*8+k);
        end
    end
end
p_permu=uint8(p_permu);

