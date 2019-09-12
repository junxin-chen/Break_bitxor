% This is Junxin Chen's implementation of Zhou's medical image cipher (encryption part) which
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

function ci=Zhou_enc(m)

x01=0.5644;
x02=0.78964;
x03=0.25874;
x04=0.63245;
%
ci1=Zhou_enc_1round(m,x01);
ci2=rot90(ci1);
ci2=Zhou_enc_1round(ci2,x02);
ci3=rot90(ci2);
ci3=Zhou_enc_1round(ci3,x03);
ci4=rot90(ci3);
ci4=Zhou_enc_1round(ci4,x04);

ci=ci4;

