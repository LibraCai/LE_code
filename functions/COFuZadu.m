function CO=COFuZadu(x,r)
% 函数名称：COFuZadu
% 函数功能：计算序列的C0复杂度
% 输入参数x，r；r为容限度；x为混沌序列
% 输出参数：C0；C0为输出的混沌序列复杂度
Y=fft(x);
Gn=mean(abs(Y).^2);
YY=zeros(1,length(x));
for i=1:length(x)
    if abs(Y(i))^2>r*Gn
        YY(i)=Y(i);
    end
end
xx=ifft(YY);
Ssum=0;
for i=1:length(x)
    Ssum=Ssum+(xx(i)-x(i))^2;
end
CO=Ssum/sum(x.^2);
