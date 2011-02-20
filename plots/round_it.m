function value=round_it(input,precision)
[m,n]=size(input)
for i=1:m
    if(min(input(i))>(1/precision))
        value(i)=round(input(i)*precision)/precision;
    else
        value(i)=input(i);
    end
end
value=transpose(value);