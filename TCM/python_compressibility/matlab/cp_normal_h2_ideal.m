function out=cp_normal_h2_ideal(T)
    R=8.314472/2.01594;
	ni=2.5 ;  
	ti=0.0;
    vi=[1.616,-0.4117,-0.792,0.758,1.217] ;
    ui=[531,751,1989,2484,6859];
	term1=0;
	term2=0;
	
	term1=ni*T.^ti;
	for i=1:length(ui)
		term2=term2+vi(i).*((ui(i)./T).^2).*(exp(-ui(i)./T))./((1-exp(-ui(i)./T)).^2);     
    end
    out=R.*(term1+term2)./T;