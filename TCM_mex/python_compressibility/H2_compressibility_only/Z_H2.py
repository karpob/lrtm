def Z=Z_H2(p_bars,T):

	p=0.1*p_bars #p bars=> MPa
        from pylab import zeros
        a=zeros(9)
        b=zeros(9)
        c=zeros(9)
	a=[0.05888460, -0.06136111,-0.002650473,0.002731125,0.001802374,-0.001150707,0.9588528e-4, -0.1109040e-6,0.1264403e-9];
	b=[1.325,1.87,2.5,2.8,2.938,3.14,3.37,3.75,4.0];
	c=[1.0,1.0,2.0,2.0,2.42,2.63,3.0,4.0,5.0]
                                                                   
	sumz=0;
	for i in range(0,8)
    		sumz=sumz+a(i)*(100/T)^b(i)*(p/1)^c(i);

	Z=1+sumz;

