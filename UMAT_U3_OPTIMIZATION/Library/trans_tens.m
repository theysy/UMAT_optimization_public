function [tens] = trans_tens(flag,ntens,tens0)
%   This function converts 1st-order tensor into 2nd-order and vice versa.
%   flag=0: Vector -> Matrix
%   flag=1: Matrix -> Matrix
%   Sequence of a tensor components:
%   S=[s11 s22 s33 s12 s13 s23];

if flag==0 % vec -> tens
    tens=zeros(3);
    if ntens==3
        tens(1,1)=tens0(1);
        tens(2,2)=tens0(2);
        tens(1,2)=tens0(3);
        tens(2,1)=tens0(3);
    else
        tens(1,1)=tens0(1);
        tens(2,2)=tens0(2);
        tens(3,3)=tens0(3);
        tens(1,2)=tens0(4);
        tens(1,3)=tens0(5);
        tens(2,3)=tens0(6);
        tens(2,1)=tens0(4);
        tens(3,1)=tens0(5);
        tens(3,2)=tens0(6);
    end
else % tens -> vec
    tens=zeros([ntens,1]);
    if ntens==3
        tens(1)=tens0(1,1);
        tens(2)=tens0(2,2);
        tens(3)=tens0(1,2);
    else
        tens(1)=tens0(1,1);
        tens(2)=tens0(2,2);
        tens(3)=tens0(3,3);
        tens(4)=tens0(1,2);
        tens(5)=tens0(1,3);
        tens(6)=tens0(2,3);
    end 
end

