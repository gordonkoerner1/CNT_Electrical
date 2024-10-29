function [Aelect]=ConnectionElect_ElectrodesL(leftContacts,nodeCount)
%erased unused quantities E, Ac, Ic, LL, G, EI
   
    ii_local=zeros(4,1); jj_local=zeros(4,1);
    Ni=(leftContacts(:));
    Nj=(nodeCount+1)*ones(size(leftContacts,2),1);
       
    ii_local=[Ni Nj Ni Nj]';
    jj_local=[Ni Ni Nj Nj]';
   
   ConnectionResistance=1/10000; %This should actually be the inverse of resistance (conductance)
                                    % was 1000
   
   Kg=zeros(4,size(leftContacts,2));
   
   Kg(1,:)= ConnectionResistance;
   Kg(2,:)= -ConnectionResistance;
   Kg(3,:)= -ConnectionResistance;
   Kg(4,:)= ConnectionResistance;
   
Aelect=sparse(jj_local(:), ii_local(:), Kg(:), nodeCount+2, nodeCount+2);
%Aelect=fsparse(jj_local(:), ii_local(:), Kg(:),[nodeCount+2, nodeCount+2]);
