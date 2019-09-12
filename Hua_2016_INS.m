%%================================================================================
%This functionto do image encryption using the reference in
%         [1]. Hua, Zhongyun, et al. "Image encryption using 2D Logistic-adjusted-Sine map." 
%              Information Sciences 339 (2016): 237-253.
%All copyrights are reserved by Zhongyun Hua. E-mial:huazyum@gmail.com
%All following source code is free to distribute, to use, and to modify
%    for research and study purposes, but absolutely NOT for commercial uses.
%If you use any of the following code in your academic publication(s), 
%    please cite the corresponding paper. 
%If you have any questions, please email me and I will try to response you ASAP.
%It worthwhile to note that all following source code is written under MATLAB R2010a
%    and that files may call built-in functions from specific toolbox(es).
%%================================================================================
%%

function varargout = Hua_2016_INS(P,para,K)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the main function to implement image cipher
% P:    the input image;
% para: operation type, 'en' or 'de';
% K:    the key, if para = 'en', it can be given or can not be given; 
%       if para = 'de', it must be given;
% varargout: if K is not given, return the result and the randomly
%            generated key; if K is given, return the result.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% to get the key
if ~exist('K','var') && strcmp(para,'en')
    K = round(rand(1,232));
    OutNum = 2;
elseif ~exist('K','var')  && strcmp(para,'de')
    error('Can not dectrypted without a key');
else
    OutNum = 1;
end
tran = @(K,low,high) sum(K(low:high).*2.^(-(1:(high-low+1))));
x0 = tran(K,1,52);
y0 = tran(K,53,104);
u0 = tran(K,105,156);
w = tran(K,157,208);

y = blkproc(K(209:232),[1,12],@(x) bi2de(x));

if strcmp(para,'en')
    P = ImageExd(P);
end
[r,c] = size(P);
[X,Y,U] = KeyGene(x0,y0,u0,w,y);
S1 = ChaoticSeq(X(1),Y(1),U(1),r,c);
S2 = ChaoticSeq(X(2),Y(2),U(2),r,c);


if max(P(:))>1
    F = 256;
else
    F = 2;
end
F=256;

%To do the encryption/decryption
 C = double(P);
switch para
    case 'en'
        C = ImageShuffling(C,S1,para);
        C = ImageSub(C,S1,para);
        C = ImageShuffling(C,S2,para);
        C = ImageSub(C,S2,para);

    case 'de'
        C = ImageSub(C,S2,para);
        C = ImageShuffling(C,S2,para);
        C = ImageSub(C,S1,para);
        C = ImageShuffling(C,S1,para);
        C = C(2:end-1,2:end-1);
end

if F == 256
    C = uint8(C);
else
    C = logical(C);
end


if OutNum == 1
    varargout{1} = C;
else
    varargout{1} = C;
    varargout{2} = K;
end
end % end of the function ImageCipher

%%
function  varargout = KeyGene(x0,y0,u0,w,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is to schedule the key 
% (x0,y0,u0):    the basic component 
% (w,y): the parameter 
% varargout: two groups of initial states for the 2D-LASM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = zeros(1,2); Y = X; U = X;
for m = 1:2
    X(m) = mod(x0+w*y(m),1);
    if X(m) == 0
        X(m) = 0.4;
    end
    Y(m) = mod(y0+w*y(m),1);
    if Y(m) == 0
        Y(m) = 0.4;
    end
    U(m) = mod(u0+w*y(m),0.4) + 0.5;
end 
varargout{1} = X;
varargout{2} = Y;
varargout{3} = U;
end
%% end of the function

function C = ImageExd(P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is add random values to surroundings of the image
% P: input image
% C: output image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if max(P(:))>1
    F = 256;
else
    F = 2;
end
[r,c] = size(P);
Ro = randi(F-1,[2,c+2]);
Co = randi(F-1,[r,2]);
C = zeros(r+2,c+2);
C(2:r+1,2:c+1) = P;
C(1,:) = Ro(1,:); C(r+2,:) = Ro(2,:);
C(2:r+1,1) = Co(:,1); C(2:r+1,c+2) = Co(:,2);
end
%% end of the function

function S = ChaoticSeq(x,y,u,r,c)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is generate chaotic sequence
% (x,y,u): initial state
% (r,c): image size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = zeros(r,c);
for m = 1:r
    for n = 1:c
        x = sin(pi*u*(y+3)*x*(1-x));
        y = sin(pi*u*(x+3)*y*(1-y));
        S(m,n) = floor(mod((x+y)*2^(30),256));
    end
end
end
%% end of the function

function C = ImageShuffling(P,S,para)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is to bit manipulation confusion 
% P: input image
% S: chaotic sequence
% para: operation type, 'en' or 'de'
% C: output image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[r,c] = size(P);
b_i = ceil(log(r*c)/log(2));
if max(P(:))>1
    b_p = 8;
else
    b_p = 1;
end
C = zeros(r,c);
T = C;
switch para
    case 'en'
        for i = 1:r
            for j = 1:c
                index_bin = dec2bin(i*c+j,b_i);
                P_bin = dec2bin(P(i,j),b_p);
                
                s = dec2bin(S(i,j),8);
                temp_bin = strcat(s,index_bin);
                temp_bin = strcat(temp_bin,P_bin);
                
                T(i,j) = bin2dec(temp_bin);
            end
        end
        [T,I] = sort(T,2);
        [T,J] = sort(T,1);
        for i = 1:r
            for j = 1:c
                t = dec2bin(T(i,j),8+b_i+b_p);
                C(i,j) = bin2dec(t(1+8+b_i:end));
            end
        end 
        
    case 'de'
        C_t = C;
        for i = 1:r
            for j = 1:c
                index_bin = dec2bin(i*c+j,b_i);
                
                s = dec2bin(S(i,j),8);
                temp_bin = strcat(s,index_bin);
                
                T(i,j) = bin2dec(temp_bin);
            end
        end
        [T,I] = sort(T,2);
        [T,J] = sort(T,1);
        for i = 1:r
            for j = 1:c
                C_t(J(i,j),j) = P(i,j);
            end
        end
        for i = 1:r
            for j = 1:c
                C(i,I(i,j)) = C_t(i,j);
            end
        end
end   
end
%% end of the function

function C = ImageSub(P,S,para)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is to bit manipulation diffusion 
% P: input image
% S: chaotic sequence
% para: operation type, 'en' or 'de'
% C: output image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if max(P(:))<2
    S(:,:) = mod(S(:,:),2);
end
[r,c] = size(P);
C = zeros(r,c);
switch para
    case 'en'
        for m = 1:r
            for n = 1:c
                if m ==1 && n == 1
                    C(m,n) = bitxor(P(m,n),bitxor(P(r,c),S(m,n)));
                elseif n == 1
                    C(m,n) = bitxor(P(m,n),bitxor(C(m-1,c),S(m,n)));
                else
                    C(m,n) = bitxor(P(m,n),bitxor(C(m,n-1),S(m,n)));
                end
            end
        end
    case 'de'
        for m = r:-1:1
            for n = c:-1:1
                if m == 1&& n == 1
                    C(m,n) = bitxor(P(m,n),bitxor(C(r,c),S(m,n)));
                elseif n == 1
                    C(m,n) = bitxor(P(m,n),bitxor(P(m-1,c),S(m,n)));
                else
                    C(m,n) = bitxor(P(m,n),bitxor(P(m,n-1),S(m,n)));
                end
            end
        end
end
end
%% end of the function


