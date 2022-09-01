function [tens2] = trans_voigt(tens1,ntens,iopt,itens)
% This function converts tensor to vector based on Voigt notation
% iopt=1: Vector -> 2nd-order Tensor
% iopt=2: 2nd-order Tensor -> Vector
% itens=1: Strain | itens=2: Stress
if itens==1
    val1=0.5;
    val2=1/val1;
else
    val1=1;
    val2=1;
end
% Tens=[s11 s22 s33 s12 s13 s23]

if iopt ==1
    tens2=zeros(3);
    if ntens==3
        tens2(1,1)=tens1(1);
        tens2(2,2)=tens1(2);
        tens2(1,2)=val1*tens1(3);
        tens2(2,1)=val1*tens1(3);
    else
        tens2(1,1)=tens1(1);
        tens2(2,2)=tens1(2);
        tens2(3,3)=tens1(3);
        tens2(1,2)=val1*tens1(4);
        tens2(2,1)=val1*tens1(4);
        tens2(1,3)=val1*tens1(5);
        tens2(3,1)=val1*tens1(5);
        tens2(2,3)=val1*tens1(6);
        tens2(3,2)=val1*tens1(6);
    end
elseif iopt==2
    tens2=zeros([1,ntens]);
    if ntens==3
        tens2(1)=tens1(1,1);
        tens2(2)=tens1(2,2);
        tens2(3)=val2*tens1(1,2);
    else
        tens2(1)=tens1(1,1);
        tens2(2)=tens1(2,2);
        tens2(3)=tens1(3,3);
        tens2(4)=val2*tens1(1,2);
        tens2(5)=val2*tens1(1,3);
        tens2(6)=val2*tens1(2,3);
    end
end

