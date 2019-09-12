% Zhou's encryption in Signal Processing 97 (2014) 172¨C182 
% The original codes are open accessible via:  https://www.fst.um.edu.mo/en/staff/documents/fstycz/Bao2014SP.rar


function ci=Zhou_enc_1round(m,x01)

% use logisitc map to replace the original map
%x2=@(x) mod(r1*(1-X0)*X0+1*(4-r1)/4*(X0/0.5*double(X0<5)+ (1-X0)/(1-0.5)*double(X0>=0.5)),1);
iteration=@(x) 3.9996*x*(1-x);

[M,N]=size(m);
mi1=zeros(M,N+1);
mi1(:,2:N+1)=m;
mi1(:,1)=floor(256*abs(rand(M,1)));
mi1=uint8(mi1);
xx=x01;
ci1=zeros(M,N+1);
for i=1:M
    ci1(i,1)=mi1(i,1);
    for j=2:N+1
        xx=iteration(xx);
        mxx=mod(floor(xx*10^10),256);
        ci1(i,j)=bitxor(bitxor(mi1(i,j),mxx),ci1(i,j-1));
    end
end
ci=uint8(ci1);
