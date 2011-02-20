load catalog.dat

load aug1
load aug3
load aug5

load dec2
load dec3
load dec6

load sep1
load sep2

load mar1
load mar3
load mar5

load marb1
load marb3
load marb5

load apr1
load apr3

load aprb1
load aprb3

subplot(2,2,1)
 m_alphaaug1=alphashape( freqaug1,0.82,0.092,0.082,1.00,298);
goodaug1=(( m_alphaaug1-alphaaug1).^2).*((sigmaaug1).^(-2));
plot(freqaug1, m_alphaaug1,'h')
title('aug1 298K 1 bar 8%')
hold on
errorbar(freqaug1,alphaaug1,sigmaaug1,'r')
hold off

subplot(2,2,2)

 m_alphaaug3=alphashape( freqaug3,0.82,0.092,0.082,3.04,298);
goodaug3=(( m_alphaaug3-alphaaug3).^2).*((sigmaaug3).^(-2));
plot(freqaug3, m_alphaaug3,'h')
title('aug3 298K 3.04 bar 8%')
hold on
errorbar(freqaug3,alphaaug3,sigmaaug3,'r')
hold off

subplot(2,2,3)
 m_alphaaug5=alphashape( freqaug5,0.82,0.092,0.082,5.08,298);
goodaug5=(( m_alphaaug5-alphaaug5).^2).*((sigmaaug5).^(-2));
plot(freqaug5, m_alphaaug5,'h')
title('aug5 298K 5.08 bar 8%')
hold on
errorbar(freqaug5,alphaaug5,sigmaaug5,'r')
hold off


figure

subplot(2,2,1)
 m_alphadec2=alphashape( freqdec2,0.82,0.092,0.082,2.02,292);
gooddec2=(( m_alphadec2-alphadec2).^2).*((sigmadec2).^(-2));
plot(freqdec2, m_alphadec2,'h')
title('dec2 292K 2.02 bar 8%')
hold on
errorbar(freqdec2,alphadec2,sigmadec2,'r')
hold off

subplot(2,2,2)
 m_alphadec3=alphashape( freqdec3,0.82,0.092,0.082,3.04,292);
gooddec3=(( m_alphadec3-alphadec3).^2).*((sigmadec3).^(-2));
plot(freqdec3, m_alphadec3,'h')
title('dec3 292K 3.04 bar 8%')
hold on
errorbar(freqdec3,alphadec3,sigmadec3,'r')
hold off

subplot(2,2,3)
 m_alphadec6=alphashape( freqdec6,0.82,0.092,0.082,6.14,292);
gooddec6=(( m_alphadec6-alphadec6).^2).*((sigmadec6).^(-2));
plot(freqdec6, m_alphadec6,'h')
title('dec6 292K 6.14 bar 8%')
hold on
errorbar(freqdec6,alphadec6,sigmadec6,'r')
hold off

figure

subplot(2,2,1)
 m_alphasep1=alphashape( freqsep1,0.82,0.092,0.082,1.01,213);
goodsep1=(( m_alphasep1-alphasep1).^2).*((sigmasep1).^(-2));
plot(freqsep1, m_alphasep1,'h')
title('sep1 213K 1.01 bar 8%')
hold on
errorbar(freqsep1,alphasep1,sigmasep1,'r')
hold off

subplot(2,2,2)
 m_alphasep2=alphashape( freqsep2,0.82,0.092,0.082,2.09,213);
goodsep2=(( m_alphasep2-alphasep2).^2).*((sigmasep2).^(-2));
plot(freqsep2, m_alphasep2,'h')
title('sep2 213K 2.09 bar 8%')
hold on
errorbar(freqsep2,alphasep2,sigmasep2,'r')
hold off


figure
subplot(2,2,1)
 m_alphamar1=alphashape( freqmar1,0.82,0.092,0.082,1.10,213);
goodmar1=(( m_alphamar1-alphamar1).^2).*((sigmamar1).^(-2));
plot(freqmar1, m_alphamar1,'h')
title('mar 213K 1.1 bar 8%')
hold on
errorbar(freqmar1,alphamar1,sigmamar1,'r')
hold off

subplot(2,2,2)
 m_alphamar3=alphashape( freqmar3,0.82,0.092,0.082,2.94,213);
goodmar3=(( m_alphamar3-alphamar3).^2).*((sigmamar3).^(-2));
plot(freqmar3, m_alphamar3,'h')
title('mar 213K 2.94 bar 8%')
hold on
errorbar(freqmar3,alphamar3,sigmamar3,'r')
hold off

subplot(2,2,3)
 m_alphamar5=alphashape( freqmar5,0.82,0.092,0.082,5.15,213);
goodmar5=(( m_alphamar5-alphamar5).^2).*((sigmamar5).^(-2));
plot(freqmar5, m_alphamar5,'h')
title('mar 213K 5.15 bar 8%')
hold on
errorbar(freqmar5,alphamar5,sigmamar5,'r')
hold off

figure
subplot(2,2,1)
 m_alphamarb1=alphashape( freqmarb1,0.88,0.098,0.022,1.12,213);
goodmarb1=(( m_alphamarb1-alphamarb1).^2).*((sigmamarb1).^(-2));
plot(freqmarb1, m_alphamarb1,'h')
title('marb 213K 1.12 bar 2%')
hold on
errorbar(freqmarb1,alphamarb1,sigmamarb1,'r')
hold off

subplot(2,2,2)

 m_alphamarb3=alphashape( freqmarb3,0.88,0.098,0.022,2.87,213);
goodmarb3=(( m_alphamarb3-alphamarb3).^2).*((sigmamarb3).^(-2));
plot(freqmarb3, m_alphamarb3,'h')
title('marb 213K 2.87 bar 2%')
hold on
errorbar(freqmarb3,alphamarb3,sigmamarb3,'r')
hold off

subplot(2,2,3)
 m_alphamarb5=alphashape( freqmarb5,0.88,0.098,0.022,5.36,213);
goodmarb5=(( m_alphamarb5-alphamarb5).^2).*((sigmamarb5).^(-2));
plot(freqmarb5, m_alphamarb5,'h')
title('marb 213K 5.37 bar 2%')
hold on
errorbar(freqmarb5,alphamarb5,sigmamarb5,'r')
hold off


figure

subplot(2,2,1)
 m_alphaapr1=alphashape( freqapr1,0.88,0.098,0.022,1.01,175);
goodapr1=(( m_alphaapr1-alphaapr1).^2).*((sigmaapr1).^(-2));
plot(freqapr1, m_alphaapr1,'h')
title('apr 175K 1.01 bar 2%')
hold on
errorbar(freqapr1,alphaapr1,sigmaapr1,'r')
hold off

subplot(2,2,2)
 m_alphaapr3=alphashape( freqapr3,0.88,0.098,0.022,3.15,175);
goodapr3=(( m_alphaapr3-alphaapr3).^2).*((sigmaapr3).^(-2));
plot(freqapr3, m_alphaapr3,'h')
title('apr 175K 3.15 bar 2%')
hold on
errorbar(freqapr3,alphaapr3,sigmaapr3,'r')
hold off

subplot(2,2,3)
 m_alphaaprb1=alphashape( freqaprb1,0.88,0.098,0.022,1.01,175);
goodaprb1=(( m_alphaaprb1-alphaaprb1).^2).*((sigmaaprb1).^(-2));
plot(freqaprb1, m_alphaaprb1,'h')
title('aprb1 175K 1.01 bar 2%')
hold on
errorbar(freqaprb1,alphaaprb1,sigmaaprb1,'r')
hold off

subplot(2,2,4)
 m_alphaaprb3=alphashape( freqaprb1,0.88,0.092,0.022,3.08,175);
goodaprb3=(( m_alphaaprb3-alphaaprb3).^2).*((sigmaaprb3).^(-2));
plot(freqaprb3, m_alphaaprb3,'h')
title('apr1 175K 3.08 bar 2%')
hold on
errorbar(freqaprb3,alphaaprb3,sigmaaprb3,'r')
hold off


goodtotal1=[goodaug1' goodaug3' goodaug5' gooddec2' gooddec3' gooddec6' goodsep1' goodsep2']'
goodtotal2=[goodmar1' goodmar3' goodmar5' goodmarb1' goodmarb3' goodmarb5' goodapr1' goodapr3' goodaprb1' goodaprb3']'