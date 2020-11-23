function [x] = infi(V,Vhalf,k)
x = 1./(1+exp((Vhalf - V)./k));
end