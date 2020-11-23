function [V,m,h,n,t] = HHmodel(I,tspan,start,over,v, mi, hi, ni)
dt = 0.001;
loop  = ceil(tspan/dt);

Ek = -12;
Ena = 120;
El = 10.6;
gk_max = 36;
gna_max = 120;
gl = 0.3;
Iapp = I;

t = (1:loop)*dt;
V = zeros(loop,1);
m = zeros(loop,1);
h = zeros(loop,1);
n = zeros(loop,1);
gna = zeros(loop,1);
gk = zeros(loop,1);
ina = zeros(loop,1);
ik = zeros(loop,1);
iall = zeros(loop,1);
ic = zeros(loop,1);

V(1)=v;
m(1)=mi;
h(1)=hi;
n(1)=ni;

for i=1:loop-1 
    if (t(i)>start)&&(t(i)<over)
        Icur = Iapp;
    else
        Icur = 0;
    end
    V(i+1) = V(i) + dt*(gna_max*m(i)^3*h(i)*(Ena-V(i)) + gk_max*n(i)^4*(Ek-V(i)) + gl*(El-V(i)) + Icur);
    m(i+1) = m(i) + dt*(alphaM(V(i))*(1-m(i)) - betaM(V(i))*m(i));
    h(i+1) = h(i) + dt*(alphaH(V(i))*(1-h(i)) - betaH(V(i))*h(i));
    n(i+1) = n(i) + dt*(alphaN(V(i))*(1-n(i)) - betaN(V(i))*n(i));
    gna(i) = gna_max*m(i)^3*h(i);
    gk(i) = gk_max*n(i)^4;
    ina(i) = -gna_max*m(i)^3*h(i)*(Ena-V(i));
    ik(i) = -gk_max*n(i)^4*(Ek-V(i));
    iall(i) = -(gna_max*m(i)^3*h(i)*(Ena-V(i)) + gk_max*n(i)^4*(Ek-V(i)) + gl*(El-V(i)));
    ic(i) = Icur;
end

figure
plot(t,V);
xlabel('Time');
ylabel('Membrane Potential');

figure
plot(t,m,t,h,t,n);
xlabel('Time');
ylabel('Gating variables');
legend('m','h','n')

figure
plot(t,gk,t,gna);
xlabel('Time');
ylabel('Conductances');
legend('gk','gna')

figure
plot(t,ik,t,ina,t,iall);
xlabel('Time');
ylabel('Currents');
legend('ik','ina','iall')

figure
plot(t,ic);
xlabel('Time');
ylabel('Appliedd current');
end
%%
function aN = alphaN(V)
aN = (0.1-0.01*V) / (exp(1-V/10) -1);
end

function bN = betaN(V)
bN = 0.125*exp(-V/80);
end

function aM = alphaM(V)
aM = (2.5-0.1*V) / (exp(2.5-0.1*V) -1);
end

function bM = betaM(V)
bM = 4*exp(-V/18);
end

function aH = alphaH(V)
aH = 0.07*exp(-V/20);
end

function bH = betaH(V)
bH = 1./(exp(3.0-0.1*V)+1);
end