% This is Junxin Chen's implementation of Fu's medical image cipher (encryption part) which
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

function p5 = Fu_enc(p)

[M,N]=size(p);
miu0=3.99976;
x00=0.92837471;
a0=40;
b0=141;
count=3;

% initial iteration 300 times for avoiding non-chaotic points 
for j=1:300
    x00=miu0*x00*(1-x00);
end

% the encryption core
for k=1:count
    
    % new parameter in each iteration so as to produce different encryption
    % elements
    miu=miu0+k*10^(-5);
    x0=x00;
    a=a0+k;
    b=b0+k;
    
    % produce the substitution mask in this iteration
    for i=1:M
    for j=1:N
        x0=miu*x0*(1-x0);
        key(i,j)=mod(floor(x0*10^14),256);        
    end
    end
    key=uint8(key);
    key2=reshape(key.',1,M*N);   
    
    
    % Fu's bit-level permutation function, given below
    p2=Fu_bitpermu_enc(p,a,b);
    

    % the substitution part using c(n)=p(n)*k(n)*c(n-1)
    c0=123; % initial value of c(n-1), i.e., the previous ciphertext for avalanche effect 
    p3=reshape(p2.',1,M*N);
    for j=1:M*N
        p4(j)=bitxor(bitxor(p3(j),key2(j)),c0);% the substitution part using c(n)=p(n)*k(n)*c(n-1)
        c0=p4(j); % renew the value of  the previous ciphertext c(n-1)
    end
    p5=reshape(p4,M,N).'; % reshape to the two-dimensional image 
    p=p5;
end


% the bit-level permutation used in Fu's cipher
function p_permu = Fu_bitpermu_enc(p,a,b)

[M,N]=size(p);
% split to 8 bit-planes
for i=1:8
    % split to 8 bit-planes
    eval (['p' num2str(8-i) '=bitand(p,2^' num2str(8-i) ');']);
    % scramble this bit-plane using cat map, the cat's parameter is
    % different for different bit-planes
    eval (['p' num2str(8-i) '=arnold_trans(p' num2str(8-i) ',a+' num2str(8-i) ',b+' num2str(8-i) ',3);']);
end
% re-combine the bit-planes into one image
p_permu=uint8(double(p7)+double(p6)+double(p5)+double(p4)+double(p3)+double(p2)+double(p1)+double(p0));

% the cat map permutation
function img_transed = arnold_trans(img,a,b,count)

[M,N] = size(img);
img_transed = zeros(M,N);

for loop = 1:count
    for x=1:M
        for y=1:N
            x1=mod((x-1)+a*(y-1),M)+1;
            y1=mod(b*(x-1)+(a*b+1)*(y-1),N)+1;
            img_transed(x1,y1)=img(x,y);
        end
    end
    img = img_transed;
end

