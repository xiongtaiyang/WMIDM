function D=FractalDim(y,cellmax)
%������һά�źŵļƺз���ά��
%y ��һά�ź�
%cellmax: �����ӵ����߳� ,����ȡ 2 ��ż�����ݴ� (1,2,4,8...),ȡ�������ݳ��ȵ�ż��
%D �� y �ļƺ�ά����һ������� D>=1��,D=lim(log(N(e))/log(k/e));
%�����ж�cellmax�ĳ����Ƿ�С�����ݵĳ��ȣ��ǣ�����
if cellmax<length(y)
error('cellmax must be larger than input signal!')
end
L=length(y);% ��������ĸ���
y_min=min(y);
%��λ�������� y_min �Ƶ����� 0 ��
y_shift=y-y_min;
%�ز�����ʹ�ܵ������� cellmax+1
%vq = interp1(x,v,xq) ʹ�����Բ�ֵ����һά�������ض���ѯ��Ĳ���ֵ������ x ���������㣬
%v ������Ӧֵ v(x)������ xq ������ѯ������ꡣ
%������ж����ͬһ��������������ݼ���
%����Խ� v ���������ʽ���д��ݡ����� v ��ÿһ�ж�����һ�鲻ͬ��һά����ֵ��
x_ord=[0:L-1]./(L-1);
xx_ord=[0:cellmax]./(cellmax);
y_interp=interp1(x_ord,y_shift,xx_ord);
%���������� y��ʹ���ֵΪ 2^^c
ys_max=max(y_interp);
factory=cellmax/ys_max;
yy=abs(y_interp*factory);
t=log2(cellmax)+1;% ��������
for e=1:t
Ne=0;%�ۻ������źŵĸ��ӵ�����
cellsize=2^(e-1);% ÿ�εĸ��Ӵ�С
NumSeg(e)=cellmax/cellsize;% ���Ữ�ֳɵĶ���
for j=1:NumSeg(e) % �ɺ����һ������ͨ�����������Խ�ĸ������ۻ� N(e)
begin=cellsize*(j-1)+1;% ÿһ�ε���ʼ
tail=cellsize*j+1;%ÿ�εĽ�β
seg=[begin:tail];% ������
yy_max=max(yy(seg));
yy_min=min(yy(seg));
up=ceil(yy_max/cellsize);%Y = ceil(X) �� X ��ÿ��Ԫ���������뵽���ڻ���ڸ�Ԫ�ص���ӽ�����
down=floor(yy_min/cellsize);%��������
Ns=up-down;% ��������ռ�еĸ�����
Ne=Ne+Ns;% �ۼ�ÿһ�θ������ߵĸ�����
end
N(e)=Ne;% ��¼ÿ e �µ� N(e)
end
%�� log(N(e)) �� log(k/e) ������С���˵�һ��������� ,б�ʾ��� D
r=-diff(log2(N));% ȥ�� r ���� 2 ��С�� 1 ��Ұ������
id=find(r<=2&r>=1);% ���������ݵ�
Ne=N(id);
e=NumSeg(id);
P=polyfit(log2(e),log2(Ne),1);% һ��������Ϸ���б�ʺͽؾ�
D=P(1);
end
