function [displacements,Fpush] = solution_inactivesAN(GDof,K,numberBeams,E,A,ang,phi,rate,oldForce,inactiveCNTs,t)
    
    % set up displacement & force of each node
    U = zeros(GDof,1);   
    F = zeros(GDof,1);

    % set up forces along each DOF for each bottom node
    Fpush = zeros(numberBeams*6,1);

    % define active cnts (removes delaminated cnts)
    activeCNTs = setdiff(1:numberBeams,inactiveCNTs);
 
    % pdof = 1:numberBeams
    pdof = activeCNTs; % only populates displacement of active CNTs

    % calculate displacements of active cnts
        % reference bottom node DOFs of active cnts
    U(GDof-6*numberBeams+(6*(pdof-1))+1) = rate(pdof).*sin(ang(pdof)).*cos(phi(pdof)); % Push nodes up
    U(GDof-6*numberBeams+(6*(pdof-1))+2) = rate(pdof).*sin(ang(pdof)).*sin(phi(pdof));
    U(GDof-6*numberBeams+(6*(pdof-1))+3) = rate(pdof).*cos(ang(pdof));
    U(GDof-6*numberBeams+(6*(pdof-1))+4) = 0; % No change in beam angles
    U(GDof-6*numberBeams+(6*(pdof-1))+5) = 0;
    U(GDof-6*numberBeams+(6*(pdof-1))+6) = 0;

    % define bottom node DOFs of inactive cnts
    inactiveDofs = GDof - 6*numberBeams + (6*(inactiveCNTs - 1)) + (1:3); % this keeps the dead CNTs from spinning
    
    % define bottom node DOFs of active cnts
    bottomDof = setdiff((GDof-6*numberBeams+1):GDof,inactiveDofs); % DOFs pushed upwards
            % allDof = [(GDof-6*numberBeams+1):GDof] --> both active & inactive bottom node DOFs
    prescribedDof = bottomDof'; % DOFs pushed upwards

    % define all upper node DOFs of active cnts, leaving the inactive bottom DOFs
        % to solve for bottom inactive nodes later
    activeDof = setdiff((1:GDof)',prescribedDof); % forces = 0 for these DOF
    
    % track calculation time: start
    % tic()

%% Displacement-based growth rather than force-based growth

    % calculate forces exerted at bottom nodes (for all DOFs)
    F(activeDof) = -K(activeDof,prescribedDof)*U(prescribedDof);

    % calculate displacements of upper nodes based on forces at bottom nodes (for all DOFs)
   
    S = tril(K(activeDof,activeDof),-1)+tril(K(activeDof,activeDof))';
    U(activeDof) = S\(F(activeDof));

    % displacements along each DOF for all nodes (active & inactive)
    displacements = U;
   
%%% BEAM ANGLE DISPLACEMENTS found to be HUGE for delaminated cnts   
%     if t == 13
%         error('delam check')
%             % cnt 171 delaminates at t=13s (stop at t=12s before reposition)
%             % cnt 73 delaminates at t=17s (stop at t=16s before reposition)
%     end
     
    % Fpush(prescribedDof) = K(prescribedDof,:)*U;

    % forces at all bottom nodes (active & inactive)
    Fpush = K((GDof-6*numberBeams+1):GDof,:)*U;

    % track calculation time: stop
    % decomp = toc()
    end