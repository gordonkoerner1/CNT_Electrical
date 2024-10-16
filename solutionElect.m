function [V] = solutionElect(GDof,prescribedDof_Elect,K_Elect,V,nodeCount,numberBeams)

    current = zeros(nodeCount+2,1);
%   activeDof=setdiff([1:nodeCount+2], prescribedDof_Elect');
    activeDof = (1:nodeCount);
 
    current(activeDof) = K_Elect(activeDof,prescribedDof_Elect)*V(prescribedDof_Elect)';
    
    S = (K_Elect(activeDof,activeDof) + K_Elect(activeDof,activeDof)').*0.5;
    V(activeDof) = S\(-current(activeDof));
