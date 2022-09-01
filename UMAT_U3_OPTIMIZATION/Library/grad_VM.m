function [dfds] = grad_VM(stress)
%% deviators
ndim3=max(size(stress));
if ndim3==3
    ndim1=2;
    ndim6=ndim3+1;
    dev_sig=zeros([ndim6,1]);
    dev_sig(1)=stress(1);
    dev_sig(2)=stress(2);
    dev_sig(3)=0;
    dev_sig(4)=stress(3);
else
    ndim1=3;
    ndim6=ndim3;
    dev_sig=zeros([ndim6,1]);
    dev_sig(:)=stress(:);
end
sig_kk=0;
delta=zeros([ndim6,1]);
for i=1:3
    sig_kk=sig_kk+dev_sig(i)/3;
    delta(i)=1;
end
dev_sig=dev_sig-sig_kk.*delta;

%% Von-Mises equivalent stress
psi=0;
for i=1:3
    psi= psi+dev_sig(i)*dev_sig(i);
end
for i=4:ndim6
    psi= psi+2*dev_sig(i)*dev_sig(i);
end
    phi=sqrt((3/2)*psi);
%% Gradient
if ndim3==3
    dev_sig0(1)=dev_sig(1);
    dev_sig0(2)=dev_sig(2);
    dev_sig0(3)=dev_sig(4);
else
    dev_sig0=dev_sig;
end
dfds=(3/2)*dev_sig0/phi;
