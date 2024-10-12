function displacements=solution_mod3_mumps(GDof,K,numberBeams,E,A,ang,phi)
  

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
   A = K(activeDof,activeDof);
   b = force(activeDof);

   % Initialize MUMPS
   id = initmumps;
   
   % Initialize MUMPS structures
   id = dmumps(id);
   
   % Set the analysis + factorization + solve option in MUMPS (JOB = 6)
   id.JOB = 6;
   
   % Disable output of MUMPS messages (optional)
   id.ICNTL(6) = 0;
   
   % Set the right-hand side
   id.RHS = b;
   
   % Call MUMPS to factorize and solve the system
   id = dmumps(id, A);
   
   % Check MUMPS output for success or issues
   % disp('INFOG(1) and INFOG(2) give information on any potential issues:');
   % disp(['INFOG(1): ', num2str(id.INFOG(1))]);
   % disp(['INFOG(2): ', num2str(id.INFOG(2))]);
   
   % Optionally activate numerical maximum transversal to improve stability
   % id.ICNTL(6) = 6;
   % id = dmumps(id, A);
   
   % Check the solution by comparing A*SOL with the right-hand side b
   % residual_norm = norm(A*id.SOL - b, 'inf');
   % disp(['Residual norm (should be small): ', num2str(residual_norm)]);
   
   % Check if the solution is close enough
   % if residual_norm < 1e-6
   %     disp('Solution is accurate!');
   % else
   %     disp('Solution may need further checks.');
   % end
   
   displacements=zeros(GDof,1);
   displacements(activeDof)=id.SOL;

   % Destroy the MUMPS instance to free resources
   id.JOB = -2;
   id = dmumps(id);
   
   % Final output to indicate the test completed
   % fprintf('Test completed successfully.\n');

