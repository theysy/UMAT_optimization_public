function [a] = smoothstep(a0, a1, t0, t1, dt)
t=linspace(t0,t1,1/dt);
a=zeros(size(t));
for i=1:1/dt
    xi=(t(i)-t0)/(t1-t0);
    a(i)=a0+(a1-a0)*xi^3*(10-15*xi+6*xi^2);
end
plot(t,a)
end

