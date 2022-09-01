function [tens3] = rotmat(ang, vec, itens)
% This function roates 2nd-order tensor.
% itens=1: Strain | itens=2: Stress
th=ang*pi/180;
r=[cos(th) sin(th) 0; -sin(th) cos(th) 0; 0 0 1];
ntens=max(size(vec));
tens1= trans_voigt(vec,ntens,1,itens); % Vec -> Tensor
tens2=r*tens1;
tens2=tens2*r';
tens3= trans_voigt(tens2,ntens,2,itens); % Tensor -> Vec
end

