function deg=radtodeg(rad)

deg=180*rad/pi;

if deg<0
   m=-1
else
   m=1
end
deg=abs(deg)

while deg>360
   deg=deg-360
end

deg=m*deg