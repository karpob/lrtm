%Plot weighting functions
clear all;
hold off;
load weighting_functions_goodman.mat

recordlength=486;
k=(0:recordlength-1);

figure(1)
axes('YTickLabel',{'0.1','1','10','100','1000'},...
    'YScale','log',...
    'YMinorTick','on',...
    'YDir','reverse');
hold('on');
%title('Nadir Viewing')
semilogy(weighting_function_a_nadir(1:486,6),tcme(recordlength-k,1),'m')
semilogy(weighting_function_a_nadir(1:486,5),tcme(recordlength-k,1),'c')
semilogy(weighting_function_a_nadir(1:486,4),tcme(recordlength-k,1),'g')
semilogy(weighting_function_a_nadir(1:486,3),tcme(recordlength-k,1),'r')
semilogy(weighting_function_a_nadir(1:486,2),tcme(recordlength-k,1),'b')
semilogy(weighting_function_a_nadir(1:486,1),tcme(recordlength-k,1),'k')
xlabel('Weighting Function')
ylabel('Pressure (bars)')
legend('1.3 cm','3.125 cm','6.25 cm','12.5 cm','25 cm','50 cm','Location','SouthEast')
title('Nadir Viewing (H_2O Goodman)')
hold off;


 figure(2)
 axes('YTicklabel',{'0.1','1','10','100','1000'},...
     'Yscale','log',...
     'YminorTick','on',...
     'Ydir','reverse')
 box('on')
 hold('on')
 
 semilogy(weighting_function_a_limb(1:486,6),tcme(recordlength-k,1),'m')
 semilogy(weighting_function_a_limb(1:486,5),tcme(recordlength-k,1),'c')
 semilogy(weighting_function_a_limb(1:486,4),tcme(recordlength-k,1),'g')
 semilogy(weighting_function_a_limb(1:486,3),tcme(recordlength-k,1),'r')
 semilogy(weighting_function_a_limb(1:486,2),tcme(recordlength-k,1),'b')
 semilogy(weighting_function_a_limb(1:486,1),tcme(recordlength-k,1),'k')
 xlabel('Weighting Function')
 ylabel('Pressure (bars)')
 legend('1.3 cm','3.125 cm','6.25 cm','12.5 cm','25 cm','50 cm','Location','SouthEast')
 title('Limb Viewing (H_2O Goodman)')
 hold off;
 
 figure(3)
axes('YTickLabel',{'0.1','1','10','100','1000'},...
    'YScale','log',...
    'YMinorTick','on',...
    'YDir','reverse');
hold('on');
%title('Nadir Viewing')
semilogy(weighting_function_a_nadir(1:486,6),tcme(recordlength-k,1),'m')
semilogy(weighting_function_a_nadir(1:486,5),tcme(recordlength-k,1),'c')
semilogy(weighting_function_a_nadir(1:486,4),tcme(recordlength-k,1),'g')
semilogy(weighting_function_a_nadir(1:486,3),tcme(recordlength-k,1),'r')
semilogy(weighting_function_a_nadir(1:486,2),tcme(recordlength-k,1),'b')
semilogy(weighting_function_a_nadir(1:486,1),tcme(recordlength-k,1),'k')
xlabel('Weighting Function')
ylabel('Pressure (bars)')
legend('1.3 cm','3.125 cm','6.25 cm','12.5 cm','25 cm','50 cm','Location','SouthEast')
title('Nadir Viewing (H_2O Goodman)')
hold off;
 
 