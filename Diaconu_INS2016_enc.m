% This is Junxin Chen's implementation of Diaconu's image cipher (encryption part) which
% is published in: Diaconu, A. V. (2016). Circular inter¨Cintra pixels bit-level permutation 
% and chaos-based image encryption. Information Sciences, 355, 314-327.

% All copyrights are reserved by Junxin Chen. E-mail:chenjx@bmie.neu.edu.cn

% All of the source codes are free to distribute, to use, and to modify
%    for research and study purposes, but absolutely NOT for commercial uses.
% If you use any of the following code in your academic publication(s), 
%    please cite the corresponding paper, as aforementioned. 
% If you have any questions, please email me and I will try to response you ASAP.
% It worthwhile to note that all following source codes are written under MATLAB R2018a.

% Interesting readers are also encouraged to refer to the original paper
% and to contact with the original authors for more details of this cipher

%
 
function p_enc = Diaconu_INS2016_enc(p)

[M,N]=size(p);
miu=3.99976;
x0=0.92837471;

p_permu=zeros(M,N);
for i=1:M
    % split the plaintext into binary planes
    for j=1:N
        for k=1:8  % 256 gra-scale image is used for illustration
            if(bitand((p(i,j)),2^(8-k)))
                p_bit(i,(j-1)*8+k)=1;
            else
               p_bit(i,(j-1)*8+k)=0; 
            end
        end
    end
    % calculate the loops of circular shift
    N_shifts=sum(p_bit(i,:));
    % calculate the direction of circular shift and then perform it
    if mod(N_shifts,2)
        p_bit(i,:)=circshift(p_bit(i,:),-N_shifts,2);
    else 
        p_bit(i,:)=circshift(p_bit(i,:),N_shifts,2);
    end
    % recombing to the permutation image
    for j=1:N
        for k=1:8
            p_permu(i,j)=p_permu(i,j)+2^(8-k)*p_bit(i,(j-1)*8+k);
        end
    end
end

p_permu=uint8(p_permu);

% generate the I_cipher_col and I_cipher_row
for i=1:300
    x0=miu*x0*(1-x0);
end
%substitution masks
for i=1:M
    for j=1:N
        x0=miu*x0*(1-x0);
        I_cipher_row(i,j)=mod(floor(x0*10^14),256);
    end
end
for i=1:M
    for j=1:N
        x0=miu*x0*(1-x0);
        I_cipher_col(i,j)=mod(floor(x0*10^14),256);
    end
end
I_cipher_row=uint8(I_cipher_row);
I_cipher_col=uint8(I_cipher_col);

% perform substitution
p_enc=bitxor(bitxor(p_permu,I_cipher_row),I_cipher_col);
p_enc=uint8(p_enc);


