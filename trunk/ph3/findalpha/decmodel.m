fmodel=1:0.1:30;
hold on
m_alphadec2=alphashape(fmodel,0.82,0.092,0.082,2.02,292);
semilogy(1:0.1:30,m_alphadec2,':k')
hold on
yerrbar(freqdec2,alphadec2,sigmadec2,1,1)

m_alphadec3=alphashape(fmodel,0.82,0.092,0.082,3.04,292);
semilogy(fmodel,m_alphadec3,'k')
hold on
yerrbar(freqdec3,alphadec3,sigmadec3,1,1)

m_alphadec6=alphashape(fmodel,0.82,0.092,0.082,6.14,292);
semilogy(fmodel, m_alphadec6,'--k')
hold on
yerrbar(freqdec6,alphadec6,sigmadec6,1,1)



ylabel('\fontsize{18}Opacity (dB/km)')
xlabel('\fontsize{18}Frequency (GHz)')