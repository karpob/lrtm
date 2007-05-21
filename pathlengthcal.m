function pathlength=pathlengthcal(thetain,delta_hght,orig_radius,n1,n2)
h=delta_hght;
a=orig_radius;

cos_thetatrans=n1*(a+h)*cos(degtorad(thetain));
pathlength=h/cos_thetatrans;
