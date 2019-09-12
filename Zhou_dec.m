% This is Junxin Chen's implementation of Zhou's medical image cipher (decryption part) which
% is published in: Signal Processing 97 (2014) 172¨C182 

% All copyrights are reserved by Junxin Chen. E-mail:chenjx@bmie.neu.edu.cn

% All of the source codes are free to distribute, to use, and to modify
%    for research and study purposes, but absolutely NOT for commercial uses.
% If you use any of the following code in your academic publication(s), 
%    please cite the corresponding paper, as aforementioned. 
% If you have any questions, please email me and I will try to response you ASAP.
% It worthwhile to note that all following source codes are written under MATLAB R2018a.

% Interesting readers are also encouraged to refer to the original paper
% and to contact with the original authors for more details of this cipher


function dec=Zhou_dec(ci)
[M,N]=size(ci);

x01=0.5644;
x02=0.78964;
x03=0.25874;
x04=0.63245;

dec4=Zhou_dec_1round(ci,x04);
dec4=rot90(dec4,3);
dec3=Zhou_dec_1round(dec4,x03);
dec3=rot90(dec3,3);
dec2=Zhou_dec_1round(dec3,x02);
dec2=rot90(dec2,3);
dec1=Zhou_dec_1round(dec2,x01);


dec=dec1;



