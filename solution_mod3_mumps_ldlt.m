function displacements=solution_mod3_mumps_ldlt(GDof,K,numberBeams,E,A,ang,phi)
  
    
    U=zeros(GDof,1);
    
    ii=zeros(numberBeams,1);
    jj=zeros(numberBeams,1);
    kkforce=zeros(numberBeams,1);
    
    force=sparse(GDof,1);
    
  for ff=1:numberBeams
    
    %forces(GDof-12*numberBeams+6*(ff-1)+1)=E*A(ff)*sin(ang(ff))*cos(phi(ff));
    ii(ff)=GDof-12*numberBeams+6*(ff-1)+1;
    jj(ff)=1;
    kkforce(ff)=E*A(ff)*sin(ang(ff))*cos(phi(ff));
    
    %forces(GDof-12*numberBeams+6*(ff-1)+2)=E*A(ff)*sin(ang(ff))*sin(phi(ff));
    ii(numberBeams+ff)=GDof-12*numberBeams+6*(ff-1)+2;
    jj(numberBeams+ff)=1;
    kkforce(numberBeams+ff)=E*A(ff)*sin(ang(ff))*sin(phi(ff));
    
    %forces(GDof-12*numberBeams+6*(ff-1)+3)=E*A(ff)*cos(ang(ff));
    ii(2*numberBeams+ff)=GDof-12*numberBeams+6*(ff-1)+3;
    jj(2*numberBeams+ff)=1;
    kkforce(2*numberBeams+ff)=E*A(ff)*cos(ang(ff));
    
  end
    force=sparse(ii',jj',kkforce',GDof,1);
    
    pdof=1:numberBeams;
    xdof=GDof-6*(pdof-1)-5; %% horizontal displacement DOF for each Beam
    ydof=GDof-6*(pdof-1)-4; %% vertical displacement DOF for each beam
    zdof=GDof-6*(pdof-1)-3;
    theta1dof=GDof-6*(pdof-1)-2;
    theta2dof=GDof-6*(pdof-1)-1;
    theta3dof=GDof-6*(pdof-1)-0;
    alldof=[xdof ydof zdof theta1dof theta2dof theta3dof];
    prescribedDof=[alldof]';
    
    activeDof=setdiff([1:GDof]',[prescribedDof]);
    
    %U=K(activeDof,activeDof)\(force(activeDof));
    
    % Set A and b
    %A = (K(activeDof,activeDof)+K(activeDof,activeDof)').*0.5;
    mat = K(activeDof,activeDof);
    b = force(activeDof);
    
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
    displacements = id.SOL;
    id.JOB = -2;
    id = dmumps(id);
