function difference=xph3shoot(x,freq,P,temp,xH2,xHe,tgt)



%freq=input('Freq (GHz)');
%freq=2.3;
%P=input('P (bars)');
%xH2=input('H2 mixing ratio:');
%xH2=0.81;
%xHe=input('He mixing ratio:');
%xHe=0.19;
%temp=input('Temp (K)');
%tgt=input('target opacity (db/km)');


%load entry.mat;
alpha=vvw_alpha(freq,xH2,xHe,x,P,temp)

difference(1)=abs(tgt-alpha);

if x<0
   difference(2)=1e20;
end

difference=sum(difference)  
