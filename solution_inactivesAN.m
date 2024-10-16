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
    mat = K(activeDof,activeDof);
    b = F(activeDof);
    
    % initialization of a matlab MUMPS structure
    id = initmumps;
    id.SYM = 1;
    
    % here JOB = -1, the call to MUMPS will initialize C 
    % and fortran MUMPS structure
    id = dmumps(id);
    % load a sparse matrix
    % load lhr01;
    % mat = Problem.A;
    % JOB = 6 means analysis+facto+solve
    
    %prob = UFget(373);
    %mat = prob.A;
    id.JOB = 6;
    %%%%%%% BEGIN OPTIONAL PART TO ILLUSTRATE THE USE OF MAXIMUM TRANSVERSAL
    id.ICNTL(7) = 2;
    id.ICNTL(6) = 7;
    id.ICNTL(8) = 8;
    id.ICNTL(14) = 30;
    % set the right hand side
    id.RHS = b;
    %call to mumps
    id = dmumps(id,mat);
    % % we see that there is a memory problem in INFOG(1) and INFOG(2)
    % id.INFOG(1)
    % id.INFOG(2)
    % % we activate the numerical maximun transversal 
    % fprintf('total number of nonzeros in factors %d\n', id.INFOG(10));
    % 
    % %%%%%%% END OPTIONAL PART %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % if(norm(mat*id.SOL - ones(size(mat,1),1),'inf') > sqrt(eps))
    % 	disp('WARNING : precision may not be OK');
    % else
    % 	disp('SOLUTION OK');
    % end
    % norm(mat*id.SOL - ones(size(mat,1),1),'inf')
    % destroy mumps instance
    U(activeDof) = id.SOL;
    id.JOB = -2;
    id = dmumps(id);
   
    %%% BEAM ANGLE DISPLACEMENTS found to be HUGE for delaminated cnts   
    %     if t == 13
    %         error('delam check')
    %             % cnt 171 delaminates at t=13s (stop at t=12s before reposition)
    %             % cnt 73 delaminates at t=17s (stop at t=16s before reposition)
    %     end
     
    % displacements along each DOF for all nodes (active & inactive)
    displacements = U;

    % forces at all bottom nodes (active & inactive)
    Fpush = K((GDof-6*numberBeams+1):GDof,:)*U;

    % track calculation time: stop
    % decomp = toc()
    end