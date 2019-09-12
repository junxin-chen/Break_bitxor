% This is Junxin Chen's implementation of Dai's medical image cipher (encryption part) which
% is published in: Y. Dai, H. Wang, and Y. Wang, ¡°Chaotic medical image encryption algorithm based 
% on bit-plane decomposition,¡± International Journal of Pattern Recognition and Artificial 
% Intelligence, vol. 30, no. 04, p. 1657001, 2016.

% All copyrights are reserved by Junxin Chen. E-mail:chenjx@bmie.neu.edu.cn

% All of the source codes are free to distribute, to use, and to modify
%    for research and study purposes, but absolutely NOT for commercial uses.
% If you use any of the following code in your academic publication(s), 
%    please cite the corresponding paper, as aforementioned. 
% If you have any questions, please email me and I will try to response you ASAP.
% It worthwhile to note that all following source codes are written under MATLAB R2018a.

% Interesting readers are also encouraged to refer to the original paper
% and to contact with the original authors for more details of this cipher


function image_cipher=Dai_enc(p)
[M,N]=size(p);
%% parameters of Dai's cipher
% Step 6
% using logisitic map to produce the parameter of Henon map
logi_miu=3.899999;
logi_x0=0.1234567890;
for j=1:N
    logi_x0=logi_miu*logi_x0*(1-logi_x0);
end
% using logisitic map to produce the parameter of Henon map

henon_p=1+logi_x0;
henon_q=0.39245;
% henon_p=1.4;
% henon_q=0.3;
logi_x0=logi_miu*logi_x0*(1-logi_x0);
henon_x=logi_x0/2;
logi_x0=logi_miu*logi_x0*(1-logi_x0);
henon_y=logi_x0/2;

henon_x=0.7;
henon_y=0.3;
% the parameter of cat map is defined in Step 3 of Dai's paper
cat_a=1;
cat_b=1;

%% Produce the substitution masks
% Step 7
ll=M*N/2;
henon_x_seq=zeros(1,ll);
henon_y_seq=zeros(1,ll);

for i=1:ll
    henon_x_next=1-henon_p*(henon_x^2)+henon_y;
    henon_y_next=henon_q*henon_x;
    henon_x_seq(1,i)=henon_x_next;
    henon_y_seq(1,i)=henon_y_next;
    henon_x=henon_x_next;
    henon_y=henon_y_next;
end
% Steps 8-10
henon_variable_seq=zeros(1,M*N);
henon_variable_seq(1,1:M*N/4)=henon_x_seq(1,1:M*N/4);
henon_variable_seq(1,M*N/4+1:M*N/2)=henon_y_seq(1,1:M*N/4);
henon_variable_seq(1,M*N/2+1:3*M*N/4)=henon_x_seq(1,M*N/4+1:M*N/2);
henon_variable_seq(1,3*M*N/4+1:M*N)=henon_y_seq(1,M*N/4+1:M*N/2);
masks=uint8(mod(floor(henon_variable_seq*1000),256));
masks=reshape(masks,M,N);
% 

%Steps 2 split the plain image to binary planes

bit_plane_12=bitand(p,192);% the highest two planes
bit_plane_34=bitand(p,48);% the folowing higher two planes
bit_plane_5678=bitand(p,15);% the lowest four planes

% Steps 3-4 permutation with cat map
bit_plane_12_permu=arnold_trans(bit_plane_12,cat_a,cat_b,100);
bit_plane_34_permu=arnold_trans(bit_plane_34,cat_a,cat_b,50);
bit_plane_5678=bit_plane_5678;

% Steps 5 recombine to an gray-scale image
image_permu=uint8(double(bit_plane_12_permu)+double(bit_plane_34_permu)+double(bit_plane_5678));

% Step 10 subsitution with XOR
image_cipher=bitxor(image_permu,masks);



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

