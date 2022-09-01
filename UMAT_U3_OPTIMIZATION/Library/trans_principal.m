function [tens4] = trans_principal(tens1, iopt)
% This function transforms stress tensor into principal stress tensors.
% iopt=1: Tensor -> Principal tensor
% iopt=2: Correction: Tensor -> Principal tensor -> Tensor
% itens=1: Strain | itens=2: Stress
itens=2;
ntens=max(size(tens1));
if ntens==3
    ndi=2;
else
    ndi=3;
end
if iopt==1 % Tensor -> Principal tensor
    tens2= trans_voigt(tens1,ntens,1,itens); % Vec -> Tensor
    tens3= eig(tens2);
    tens4= sort(tens3, 'descend');
else % Correction: Tensor -> Principal tensor -> Tensor
    s=zeros([ntens,1]);
    tens32=zeros([ntens,1]);
    s(1)=1;
    tens2= trans_voigt(tens1,ntens,1,itens); % Vec -> Tensor
    [pvec, tens3]=eig(tens2);
    [tens31, indx]= sort(diag(tens3), 'descend');
    tens32(1:ndi)=tens31(1:ndi);
    r=[pvec(:,indx(1)), pvec(:,indx(2)), pvec(:,indx(3))];
    tens32=tens32.*s;
    tens33= trans_voigt(tens32,ntens,1,itens); % Vec -> Tensor
    tens33=r*tens33;
    tens33=tens33*r';
    tens4= trans_voigt(tens33,ntens,2,itens); % Tensor -> Vec
end
end

