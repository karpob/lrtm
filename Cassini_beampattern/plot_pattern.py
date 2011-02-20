from pylab import *
from scipy.io import loadmat
rc('text',usetex='True')
figure(1,figsize=(8,8))

ax_x_position=0.1
ax_y_position=0.23
ax_width=0.8
ax_height=0.68
ax=axes([ax_x_position,ax_y_position,ax_width,ax_height]);
beam_data=loadmat('beam_pattern.mat');
ColorMap=contourf(beam_data['x'],beam_data['y'],beam_data['data'],1024,antialiased=False)
xlim([-1.01,1.01])
ylim([-1.01,1.01])
xticks((-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1),('-1','-0.75','-0.5','-0.25','0','0.25','0.5','0.75','1'),fontsize=16)
yticks((-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1),('-1','-0.75','-0.5','-0.25','0','0.25','0.5','0.75','1'),fontsize=14)
xlabel('Azimuth ($^{\circ}$)',fontsize=16)
ylabel('Elevation ($^{\circ}$)',fontsize=16)
title('Cassini Radar Antenna Pattern (Beam 3)',fontsize=20)

colorbar_x_position=ax_x_position;
colorbar_y_position=0.05;
colorbar_width=ax_width;
colorbar_height=0.03
colorbar_axes=axes([colorbar_x_position,colorbar_y_position,colorbar_width,colorbar_height])
colorbar(mappable=ColorMap,cax=colorbar_axes,orientation='horizontal',format='%.1f')
xticks(fontsize=16)
text(0.3,1.27,'Normalized Antenna Gain (Linear)',fontsize=16)
#show()

savefig('antenna_pattern.png',dpi=600)

figure(2)
plot(beam_data['x'],beam_data['data'][599][:])
xlim([-1,1])
plot(beam_data['y'],beam_data['data'][:,201],'r')
xlim([-1,1])
legend(('Azimuth','Elevation'),'best')
ylabel('Gain (Linear/Normalized)',fontsize=16)
xlabel('Angle ($^{\circ}$)',fontsize=16)
xticks((-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1),('-1','-0.75','-0.5','-0.25','0','0.25','0.5','0.75','1'),fontsize=16)
yticks((0,0.25,0.5,0.75,1.0),('0','0.25','0.5','0.75','1'))
savefig('az.eps',dpi=600)

