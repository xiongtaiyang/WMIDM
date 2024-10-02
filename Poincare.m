function [Poincare_Index]=Poincare(rr_segement)
[COL,ROW]=size(rr_segement);
 T=ROW-1;
for i=1:T
temp1=sum((rr_segement(i)-rr_segement(i+1))^2./2);
temp2=sum((rr_segement(i)+rr_segement(i+1)-2*mean(rr_segement))^2./2);
end
SD1=sqrt(temp1./T);
SD2=sqrt(temp2./T);
Poincare_Index=SD1./SD2;
end