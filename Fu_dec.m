% This is Junxin Chen's implementation of Fu's medical image cipher  (decryption part) which
% is published in: C. Fu, W.-H. Meng, Y.-F. Zhan, Z.-L. Zhu, F. C. Lau, K. T. Chi, and
% H.-F. Ma, ¡°An efficient and secure medical image protection scheme
% based on chaotic maps,¡± Computers in Biology and Medicine, vol. 43,
% no. 8, pp. 1000¨C1010, 2013.

% All copyrights are reserved by Junxin Chen. E-mail:chenjx@bmie.neu.edu.cn

% All of the source codes are free to distribute, to use, and to modify
%    for research and study purposes, but absolutely NOT for commercial uses.
% If you use any of the following code in your academic publication(s), 
%    please cite the corresponding paper, as aforementioned. 
% If you have any questions, please email me and I will try to response you ASAP.
% It worthwhile to note that all following source codes are written under MATLAB R2018a.

% Interesting readers are also encouraged to refer to the original paper
% and to contact with the original authors for more details of this cipher


function d5 = Fu_dec(d1)

[M,N]=size(d1);
miu0=3.99976;
x00=0.92837471;
a0=40;
b0=141;
count=3;

% initial iteration 300 times for avoiding non-chaotic points 
% the same with that in the encryption part 
for j=1:300
    x00=miu0*x00*(1-x00);
end

% the iteration of the decryption
% please refer to the encrypttion file for better understanding of the
% procedures of the decryption process
for k=1:count
    c0=123;
    d2=reshape(d1.',1,M*N);
   
    % update the control parameter for subsitution mask and permutation
    % vector
    miu=miu0+(count+1-k)*10^(-5);
    a=a0+count+1-k;
    b=b0+count+1-k;
    
    % produce the substitution mask
    x0=x00;
    for i=1:M
    for j=1:N        
        x0=miu*x0*(1-x0);
        key(i,j)=mod(floor(x0*10^14),256);        
    end
    end
    key=uint8(key);
    key2=reshape(key.',1,M*N);   
    % decryption of the substitution part
    for j=1:M*N
        tmp=d2(j);
        d3(j)=bitxor(bitxor(d2(j),key2(j)),c0);
        c0=tmp;
    end
    d4=reshape(d3,M,N).';    
    % decryption of the bit-level permutation part
    d5=Fu_bitpermu_dec(d4,a,b);
    d1=d5; 
end

% decryption of the bit-level permutation part
function p_permu = Fu_bitpermu_dec(p,a,b)

[M,N]=size(p);

for i=1:8 % assume the gray-level is 8, without loss of the generality 
    % split the original ciphertext into bit-planes
    eval (['p' num2str(8-i) '=bitand(p,2^' num2str(8-i) ');']);
    % decryption of the permutation part 
    eval (['p' num2str(8-i) '=inv_arnold_trans(p' num2str(8-i) ',a+' num2str(8-i) ',b+' num2str(8-i) ',3);']);
end
p_permu=uint8(double(p7)+double(p6)+double(p5)+double(p4)+double(p3)+double(p2)+double(p1)+double(p0));

% the decryption of cat map
function img_transed = inv_arnold_trans(img,a,b,count)

[M,N] = size(img);
img_transed = zeros(M,N);

for loop = 1:count
    for x=1:M
        for y=1:N
            x1=mod((x-1)*(a*b+1)-a*(y-1),M)+1;
            y1=mod(-b*(x-1)+(y-1),N)+1;
            img_transed(x1,y1)=img(x,y);
        end
    end
    img = img_transed;
end
img_transed = uint8(img_transed);


