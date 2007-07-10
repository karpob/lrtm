
clear;
f_in_GHz=1:0.1:1000;
T_in_K=293.15

for k=1:length(f_in_GHz)
    c_dielec_constant_room_temp(k)=get_complex_dielectric_constant_water(f_in_GHz(k),T_in_K);
    alpha_20(k)=rayleigh_absorption(f_in_GHz(k),1,1,c_dielec_constant_room_temp(k));
end

T_in_K=273.15
for k=1:length(f_in_GHz)
    c_dielec_constant_near_zero_C(k)=get_complex_dielectric_constant_water(f_in_GHz(k),T_in_K);
    alpha_zero(k)=rayleigh_absorption(f_in_GHz(k),1,1,c_dielec_constant_near_zero_C(k));
end

figure(1)
set(gcf,'PaperSize',[5.75,3])
loglog(f_in_GHz,real(c_dielec_constant_room_temp),'k')
hold on;
loglog(f_in_GHz,real(c_dielec_constant_near_zero_C),'--b')
title('Real Part of Dielectric Constant (\epsilon\prime)','fontsize',16)
xlabel('Frequency (GHz)','fontsize',16)
ylabel('\epsilon\prime','fontsize',16)
set(gca,'XTickLabel',{'1';'10';'100';'1000'})
set(gca,'FontSize',16)
legend('Temperature 20^{\circ}C','Temperature 0^{\circ}C')
savefig('Real_cmyk','pdf','eps','-cmyk','-r600')
savefig('Real','pdf','eps','-r600')



figure(2)
set(gcf,'PaperSize',[5.75,3])
loglog(f_in_GHz,imag(c_dielec_constant_room_temp),'k')
hold on;
loglog(f_in_GHz,imag(c_dielec_constant_near_zero_C),'--b')
title('Imaginary Part of Dielectric Constant (\epsilon\prime\prime)','fontsize',16)
xlabel('Frequency (GHz)','fontsize',16)
ylabel('\epsilon\prime\prime','fontsize',16)
set(gca,'XTickLabel',{'1';'10';'100';'1000'})
set(gca,'FontSize',16)
legend('Temperature 20^{\circ}C','Temperature 0^{\circ}C')
savefig('Imag_cmyk','pdf','eps','-cmyk','-r600')
savefig('Imag','pdf','eps','-r600')

figure(3)
set(gcf,'PaperSize',[5.75,3])
loglog(f_in_GHz,alpha_20,'k')
hold on
loglog(f_in_GHz,alpha_zero,'--b')
title('Imaginary Part of Dielectric Constant (\epsilon\prime\prime)','fontsize',16)
xlabel('Frequency (GHz)','fontsize',16)
ylabel('Absorption Coefficient \alpha (cm^{-1})','fontsize',16)
set(gca,'XTickLabel',{'1';'10';'100';'1000'})
set(gca,'FontSize',16)
legend('Temperature 20 ^{\circ}C','Temperature 0 ^{\circ}C','Location','SouthEast')
savefig('alpha_cmyk','pdf','eps','-cmyk','-r600')
savefig('alpha','pdf','eps','-r600')
