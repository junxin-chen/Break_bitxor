% This is the source code for a chosen-ciphertext attack which is used to
% crack Diaconu's image cipher which
% is published in: Diaconu, A. V. (2016). Circular inter¨Cintra pixels bit-level permutation 
% and chaos-based image encryption. Information Sciences, 355, 314-327.

% This attack is based on the chosen ciphertext attack (Cracker.m in this directory, 
% the simplified version) pulished in 
% 'Cryptanalysis of some image ciphers using permutation-substitution
% network and chaos' (IEEE Transactions on Circuits and Systems for Video
% Technology, vol **, no **, pp **-**, 2019).

% There are three modifications, in comparison with the basic cipher. The
% modifications have been annotated, as can be observed in the following 'the first/second/third difference ......

% All copyrights are reserved by Junxin Chen. E-mail:chenjx@bmie.neu.edu.cn
% All of the source codes are free to distribute, to use, and to modify
%    for research and study purposes, but absolutely NOT for commercial uses.
% If you use any of the following code in your academic publication(s), 
%    please cite the corresponding paper, as aforementioned. 
% If you have any questions, please email me and I will try to response you ASAP.
% It worthwhile to note that all following source codes are written under MATLAB R2018a.


clear 
clc
close all

mm=imread('lenna32.bmp');
[M,N]=size(mm);

%% function handles of the ciphers to be cryptanalyzed
encrypt=@(m)Diaconu_INS2016_enc(m);
decrypt=@(m)Diaconu_INS2016_dec(m);
%

%% the ciphertext to be recovered
%´ý»Ö¸´µÄÃÜÎÄ
cc=encrypt(mm);
%%  the attack procedures
c0=uint8(zeros(M,N));

m0=Diaconu_crack_func_self_permu_enc(decrypt(c0)); % the first difference when cracking Diaconu's cipher

sum_del_m=uint8(zeros(M,N));
for i=1:M
    i
    for j=1:N      
        for k=1:8 
            lambda=bitand(cc(i,j),2^(8-k));
            if lambda==0
                continue;
            else                 
                    current_c=uint8(zeros(M,N));
                    current_c(i,j)=2^(8-k);
                    current_c=uint8(current_c);
                    
                    current_m=Diaconu_crack_func_self_permu_enc(decrypt(current_c));% % the second difference when cracking Diaconu's cipher
                    
                    current_del_m=bitxor(current_m,m0);
                    sum_del_m=bitxor(sum_del_m,current_del_m);                 
            end
        end
    end
end
mm_crack=bitxor(sum_del_m,m0);

mm_crack=Diaconu_crack_func_self_permu_dec(mm_crack);% the third difference when cracking Diaconu's cipher

dd=double(mm_crack)-double(mm);


