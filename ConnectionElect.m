function [Aelect] = ConnectionElect(closeNodes,nodeCount)
% erased unused quantities E, Ac, Ic, LL, G, EI
   
    sizeClose = size(closeNodes,1);

    ii_local = zeros(4,1); jj_local = zeros(4,1);
   
    Ni = (closeNodes(:,1));
    Nj = (closeNodes(:,2));
       
    ii_local = [Ni Nj Ni Nj]';
    jj_local = [Ni Ni Nj Nj]';
   
    ConnectionResistance = (1/50000)*ones(1,sizeClose); % This should actually be the inverse of resistance (conductance)
                        % was 1000
   
    Kg = zeros(4,sizeClose);
   
    Kg(1,:) = ConnectionResistance;
    Kg(2,:) = -ConnectionResistance;
    Kg(3,:) = -ConnectionResistance;
    Kg(4,:) = ConnectionResistance;
   
    Aelect = sparse(jj_local(:), ii_local(:), Kg(:), nodeCount+2, nodeCount+2);
    %Aelect = fsparse(jj_local(:), ii_local(:), Kg(:),[nodeCount+2, nodeCount+2]);
