% This is the source code for a chosen-ciphertext attack which is given in
% 'Cryptanalysis of some image ciphers using permutation-substitution
% network and chaos' (IEEE Transactions on Circuits and Systems for Video
% Technology, vol **, no **, pp **-**, 2019).

% All copyrights are reserved by Junxin Chen. E-mail:chenjx@bmie.neu.edu.cn
% All of the source codes are free to distribute, to use, and to modify
%    for research and study purposes, but absolutely NOT for commercial uses.
% If you use any of the following code in your academic publication(s), 
%    please cite the corresponding paper, as aforementioned. 
% If you have any questions, please email me and I will try to response you ASAP.
% It worthwhile to note that all following source codes are written under MATLAB R2018a.

%% Start 
clear 
clc
close all

%% function handles of the ciphers to be cryptanalyzed
% if you want to breack other ciphers, please change this part 

% % for the basic cipher
% encrypt=@(m)basic_enc(m);
% decrypt=@(m)basic_dec(m);

% for Dai's cipher: International Journal of Pattern Recognition and Artificial 
% Intelligence, vol. 30, no. 04, p. 1657001, 2016.
encrypt=@(m)Dai_enc(m);
decrypt=@(m)Dai_dec(m);

% for Fu's cipher: Computers in Biology and Medicine, vol. 43,
% no. 8, pp. 1000¨C1010, 2013.
% encrypt=@(m)Fu_enc(m);
% decrypt=@(m)Fu_dec(m);

% % for Hua's cipher: Information Sciences 339 (2016): 237-253.

% load Hua_INS_K
% encrypt=@(m)(Hua_2016_INS(m,'en',K));
% decrypt=@(m)(Hua_2016_INS(m,'de',K));


% % for Zhou's cipher:  Signal Processing 97 (2014) 172¨C182 
% encrypt=@(m)Zhou_enc(m);
% decrypt=@(m)Zhou_dec(m);

%% the ciphertext to be recovered
mm=imread('lenna16.bmp');  % the original file to be encrypted. 
% mm=imread('lenna32.bmp');  % the original file to be encrypted. 
% mm=imread('lenna256.bmp');  % the original file to be encrypted. 
% With respected to your computer's computational capacity, you can choose
% files with various sizes. 

cc=encrypt(mm); % this is the ciphertext to be recovered

% % mm=[45,129;199,235];      % for experimental illustration in Section III-D
% % cc=encrypt(mm);           % for experimental illustration in Section III-D


[M,N]=size(cc); % get the ciphertext's size 
[M_plain,N_plain]=size(mm); %get the plaintext's size 
% The sizes of the ciphertext and plaintext maybe different, for example
% Hua's cipher in Hua, Zhongyun, et al. "Image encryption using 2D Logistic-adjusted-Sine map." 
%              Information Sciences 339 (2016): 237-253.

L_gray=8; %256-gray-level image is used, without loss of the generality


%% Algorithm 1: Generate the atoms of the attack
c0=uint8(zeros(M,N));
atom_m0=uint8(decrypt(c0));
atom_delta_m=uint8(zeros(M,N,L_gray,M_plain,N_plain));
for i=1:M
    i
    for j=1:N
        for l=1:L_gray
            current_c=uint8(zeros(M,N));
            current_c(i,j)=2^(8-l);%current chosen-ciphertext
            current_c=uint8(current_c);
            current_m=uint8(decrypt(current_c));
            atom_delta_m(i,j,l,:,:)=bitxor(current_m,atom_m0);
        end
    end
end

%% Algorithm 2: Recovery of the plaintext
mm_crack=uint8(zeros(M_plain,N_plain));
for i=1:M
    i
    for j=1:N
        for l=1:8 
            bit_C=bitand(cc(i,j),2^(8-l));
            if bit_C==0
                continue; % if bit_C is zero, it is not necessary to conduct bit-wise xor
            else
                current_atom(:,:)=atom_delta_m(i,j,l,:,:);
                mm_crack=bitxor(mm_crack,current_atom);
            end
        end
    end
end
mm_crack=bitxor(mm_crack,atom_m0);
%% Check the difference, and verify the attack

dd=double(mm_crack)-double(mm);
max(max(dd))

%% The following is a simplified/fast algorithm
% corresponds to the aforementioend Algorithm 1 and Algorithm 2
% because if we just want to validate the attack, not all of the atoms have
% to be obtained, only those corresponds to bit_C!=0 are necessary


% c0=uint8(zeros(M,N));
% m0=uint8(decrypt(c0));
% sum_del_m=uint8(zeros(M_plain,N_plain));
% for i=1:M
%     i
%     for j=1:N
%         for k=1:8 
%             lambda=bitand(cc(i,j),2^(8-k));
%             if lambda==0
%                 continue;
%             else                 
%                     current_c=uint8(zeros(M,N));
%                     current_c(i,j)=2^(8-k);%
%                     current_c=uint8(current_c);
%                     current_m=uint8(decrypt(current_c));
%                     current_del_m=bitxor(current_m,m0);
%                     sum_del_m=bitxor(sum_del_m,current_del_m);                 
%             end
%         end
%     end
% end
% mm_crack=bitxor(sum_del_m,m0);
% dd=double(mm_crack)-double(mm);
% max(max(dd))