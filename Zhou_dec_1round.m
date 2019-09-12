% Zhou's decryption in Signal Processing 97 (2014) 172¨C182 
% The original codes are open accessible via:  https://www.fst.um.edu.mo/en/staff/documents/fstycz/Bao2014SP.rar


function dec=Zhou_dec_1round(ci,x01)
[M,N]=size(ci);

% use logisitc map to replace the original map
%x2=@(x) mod(r1*(1-X0)*X0+1*(4-r1)/4*(X0/0.5*double(X0<5)+ (1-X0)/(1-0.5)*double(X0>=0.5)),1);
iteration=@(x) 3.9996*x*(1-x);
% decryption of the 1st round encryption
ci1=ci;
xx=x01;
ci1_dec=zeros(M,N);
for i=1:M
    ci1_dec(i,1)=ci1(i,1);
    pre_cipher=ci1(i,1);
    for j=2:N
        tmp_cipher=ci1(i,j);
        xx=iteration(xx);
        mxx=mod(floor(xx*10^10),256);
        ci1_dec(i,j)=bitxor(bitxor(ci1(i,j),mxx),pre_cipher);
        pre_cipher=tmp_cipher;
    end
end
dec=uint8(ci1_dec(1:M,2:N));

