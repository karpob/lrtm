

%Xo=6.222e-7;
%x=Xo;
%freq=2.3;
%P=6;
%temp=213;
%xH2=0.87;
%xHe=0.13;
%tgt=1e-6;
load cass;
freq=cass(:,1);
P=cass(:,2);
temp=cass(:,3);
xH2=cass(:,4);
xHe=cass(:,5);
tgt=cass(:,6);
k=size(cass,1);

options(1)=0;
options(2)=1e-12;
options(3)=1e-12;

for i=1:k
   x(i) = fmin('xph3shoot',0,100,options,freq(i),P(i),temp(i),xH2(i),xHe(i),tgt(i));
	cd ..
	cd findalpha
	alphacheck(i)=alphashape(freq(i),xH2(i),xHe(i),x(i),P(i),temp(i));
	cd ..
	cd findph3
  
end
%fprintf(1,'Target opacity is %4.2e \n',tgt(i))
%fprintf(1,'Calculated opacity is %4.2e \n',alphacheck(i))
[tgt alphacheck']