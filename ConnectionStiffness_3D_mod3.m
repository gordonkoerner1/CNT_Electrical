function [Ac,vdwk]=ConnectionStiffness_3D_mod3(closeNodes,GDof,beamType,nodeCoordinates)

    vdk=1000;
    
    
     NumberBars=size(closeNodes,1);
     xb=nodeCoordinates(closeNodes(:,2),1)-nodeCoordinates(closeNodes(:,1),1);
     yb=nodeCoordinates(closeNodes(:,2),2)-nodeCoordinates(closeNodes(:,1),2);
     zb=nodeCoordinates(closeNodes(:,2),3)-nodeCoordinates(closeNodes(:,1),3);
     ndf=3;%number of DOFs per node in the bar element
     ntt=2*NumberBars*ndf;%total number of DOFs
 
          Ab=sparse(1,1,0,GDof,GDof);
             ll=(xb.*xb + yb.*yb + zb.*zb).^0.5;
             CXx=xb./ll;
             CYx=yb./ll;
             CZx=zb./ll;
    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   alpha=1;
     if beamType==0
         epsilon=0;
     else
         epsilon=5000; %12.*E.*Ic./(G.*(Ac./alpha));
     end
 
     vdwk=1000;
 
 
     ii_local=zeros(36,1); jj_local=zeros(36,1);
     
     Ni=6*(closeNodes(:,1)-1);
     Nj=6*(closeNodes(:,2)-1);
         
 ii_local=[Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6]';
  
 jj_local=[Ni+1 Ni+1 Ni+1 Ni+1 Ni+1 Ni+1 Ni+1 Ni+1 Ni+1 Ni+1 Ni+1 Ni+1 Ni+2 Ni+2 Ni+2 Ni+2 Ni+2 Ni+2 Ni+2 Ni+2 Ni+2 Ni+2 Ni+2 Ni+2 Ni+3 Ni+3 Ni+3 Ni+3 Ni+3 Ni+3 Ni+3 Ni+3 Ni+3 Ni+3 Ni+3 Ni+3 Ni+4 Ni+4 Ni+4 Ni+4 Ni+4 Ni+4 Ni+4 Ni+4 Ni+4 Ni+4 Ni+4 Ni+4 Ni+5 Ni+5 Ni+5 Ni+5 Ni+5 Ni+5 Ni+5 Ni+5 Ni+5 Ni+5 Ni+5 Ni+5 Ni+6 Ni+6 Ni+6 Ni+6 Ni+6 Ni+6 Ni+6 Ni+6 Ni+6 Ni+6 Ni+6 Ni+6 Nj+1 Nj+1 Nj+1 Nj+1 Nj+1 Nj+1 Nj+1 Nj+1 Nj+1 Nj+1 Nj+1 Nj+1 Nj+2 Nj+2 Nj+2 Nj+2 Nj+2 Nj+2 Nj+2 Nj+2 Nj+2 Nj+2 Nj+2 Nj+2 Nj+3 Nj+3 Nj+3 Nj+3 Nj+3 Nj+3 Nj+3 Nj+3 Nj+3 Nj+3 Nj+3 Nj+3 Nj+4 Nj+4 Nj+4 Nj+4 Nj+4 Nj+4 Nj+4 Nj+4 Nj+4 Nj+4 Nj+4 Nj+4 Nj+5 Nj+5 Nj+5 Nj+5 Nj+5 Nj+5 Nj+5 Nj+5 Nj+5 Nj+5 Nj+5 Nj+5 Nj+6 Nj+6 Nj+6 Nj+6 Nj+6 Nj+6 Nj+6 Nj+6 Nj+6 Nj+6 Nj+6 Nj+6]';
 
       Kg=zeros(144,size(closeNodes,1));
    
   
       Kg(1,:)=vdwk*(  CXx.^2  )';
       Kg(2,:)=vdwk*(  CXx.*CYx  )';
       Kg(3,:)=vdwk*(  CXx.*CZx  )';
       Kg(7,:)=vdwk*(  -CXx.^2  )';
       Kg(8,:)=vdwk*(  -CXx.*CYx  )';
       Kg(9,:)=vdwk*(  -CXx.*CZx  )';
       
       Kg(13,:)=vdwk*(  CXx.*CYx  )';
       Kg(14,:)=vdwk*(  CYx.^2  )';
       Kg(15,:)=vdwk*(  CYx.*CZx  )';
       Kg(19,:)=vdwk*(  -CXx.*CYx  )';
       Kg(20,:)=vdwk*(  -CYx.^2  )';
       Kg(21,:)=vdwk*(  -CYx.*CZx  )';
       
       Kg(25,:)=vdwk*(  CXx.*CZx  )';
       Kg(26,:)=vdwk*(  CYx.*CZx  )';
       Kg(27,:)=vdwk*(  CZx.^2  )';
       Kg(31,:)=vdwk*(  -CXx.*CZx  )';
       Kg(32,:)=vdwk*(  -CYx.*CZx  )';
       Kg(33,:)=vdwk*(  -CZx.^2  )';
       
       Kg(73,:)=vdwk*(  -CXx.^2  )';
       Kg(74,:)=vdwk*(  -CXx.*CYx  )';
       Kg(75,:)=vdwk*(  -CXx.*CZx  )';
       Kg(79,:)=vdwk*(  CXx.^2  )';
       Kg(80,:)=vdwk*(  CXx.*CYx  )';
       Kg(81,:)=vdwk*(  CXx.*CZx  )';
       
       Kg(85,:)=vdwk*(  -CXx.*CYx  )';
       Kg(86,:)=vdwk*(  -CYx.^2  )';
       Kg(87,:)=vdwk*(  -CYx.*CZx  )';
       Kg(91,:)=vdwk*(  CXx.*CYx  )';
       Kg(92,:)=vdwk*(  CYx.^2  )';
       Kg(93,:)=vdwk*(  CYx.*CZx  )';
       
       Kg(97,:)=vdwk*(  -CXx.*CZx  )';
       Kg(98,:)=vdwk*(  -CYx.*CZx  )';
       Kg(99,:)=vdwk*(  -CZx.^2  )';
       Kg(103,:)=vdwk*(  CXx.*CZx  )';
       Kg(104,:)=vdwk*(  CYx.*CZx  )';
       Kg(105,:)=vdwk*(  CZx.^2  )';
   
       
   
    %Ac=sparse(jj_local(:), ii_local(:), Kg(:), GDof, GDof);
 %%AA=sparse(ii_local(:), jj_local(:), Kg(:), GDof, GDof);
    Ac=fsparse(jj_local(:), ii_local(:), Kg(:),[GDof,GDof]);
 
