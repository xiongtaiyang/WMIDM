function [ fractality ] = deltahurst( rr_segement )
rr_length=length(rr_segement);
scales=[];
hurst_local=[];
for i=70:200 %�����512hz���ȡ��ΧΪ4-10��
    scales(i-69)=2*1.01^i;
    %scales(i-3)=2*1.2^i;
end
coef=cwt(rr_segement,scales,'gaus3');
k=max(floor(scales))+1;
for i=k:rr_length-k  %%ȥ��RR�������ͷ��β���������С���任�߶ȵ�С��ϵ�����ݵ㣬�ų�С���任�ı߽�ЧӦ
    coef1=coef(:,i);
    hurst_local(i-1)=abs(max(coef1));
end
delta_hurst=max(hurst_local)-min(hurst_local);
fractality=[delta_hurst];
end