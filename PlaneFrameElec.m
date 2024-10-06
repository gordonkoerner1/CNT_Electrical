function K_Elect = PlaneFrameElec(A,L,elementNodes,GDof,conductivity,nodeCount,numberBeams)
% Heat Transfer Conduction Matrix
    %actives=setdiff(1:nodeCount-numberBeams,tops)';
    actives = (1:nodeCount-numberBeams);
    %actives=setdiff(1:nodeCount-numberBeams,crossedElements);
   
    ii_local_heat = zeros(4,1); jj_local_heat = zeros(4,1);
    Ni = (elementNodes(actives,1));
    Nj = (elementNodes(actives,2));
   
    ii_local_heat = [Ni Nj Ni Nj]';
    jj_local_heat = [Ni Ni Nj Nj]';
   
    g = conductivity*A./L; %Conductivity - the inverse of resistivity
    g = conductivity.*L(actives); %Conductivity - the inverse of resistivity
    Gg(1,:) = g';
    Gg(2,:) = -g';
    Gg(3,:) = -g';
    Gg(4,:) = g';
    
    K_Elect = sparse(jj_local_heat(:), ii_local_heat(:), Gg(:), nodeCount+2, nodeCount+2);