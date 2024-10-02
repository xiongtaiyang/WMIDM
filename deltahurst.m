function [ fractality ] = deltahurst( rr_segement )
rr_length=length(rr_segement);
scales=[];
hurst_local=[];
for i=70:200 %如果是512hz则可取范围为4-10；
    scales(i-69)=2*1.01^i;
    %scales(i-3)=2*1.2^i;
end
coef=cwt(rr_segement,scales,'gaus3');
k=max(floor(scales))+1;
for i=k:rr_length-k  %%去掉RR间隔序列头和尾部等于最大小波变换尺度的小波系数数据点，排除小波变换的边界效应
    coef1=coef(:,i);
    hurst_local(i-1)=abs(max(coef1));
end
delta_hurst=max(hurst_local)-min(hurst_local);
fractality=[delta_hurst];
end