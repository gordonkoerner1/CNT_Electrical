function [K,LL] = SpaceBeamElemStiffness_mod4(E,A,GDof,Ic,G,elementNodes,element,nodeCoordinates,nodeCount,PeriodicBoundary,h_span,v_span) 
    %SpaceBeamElementStiffness This function returns the element stiffness
    %matrix for a space beam element with modulus of elasticity
    %E,cross-sectional area A,length L.
    %The size of the element stiffness matrix is 12x12
    
    
    % % %      xa=zeros(element,1);
    % % %      ya=zeros(element,1);
    % % %      za=zeros(element,1);
            
        xa=nodeCoordinates(elementNodes(:,2),1)-nodeCoordinates(elementNodes(:,1),1); % x component of the length of the element where elementNodes is the connectivity matrix   
        ya=nodeCoordinates(elementNodes(:,2),2)-nodeCoordinates(elementNodes(:,1),2); % y component of the length of the element 
        za=nodeCoordinates(elementNodes(:,2),3)-nodeCoordinates(elementNodes(:,1),3); % z component of the length of the element
        
              if PeriodicBoundary>0;
                 crossedl=find(xa>h_span/2);% moves across left boundary
                 xa(crossedl)=xa(crossedl)-h_span;     
                 crossedr=find( xa <(-h_span/2 ) );
                 xa(crossedr)=xa(crossedr)+h_span;
                  
                 crossed2=find(ya>v_span/2);% moves across left boundary
                 ya(crossed2)=ya(crossed2)-v_span;     
                 crossedr2=find( ya <(-v_span/2 ) );
                 ya(crossedr2)=ya(crossedr2)+v_span;
                
             end
         
    % % %       NDF=6; % number of DOFs per node
    % % %       NTT=nodeCount*NDF;% total number of DOFs
          %K=zeros(NTT,NTT);
    % % %       tt=zeros(3,3);
    % % %       localx=zeros(1,3);
    % % %       localy=zeros(1,3);
    % % %       localz=zeros(1,3);
    
        LL=(xa.*xa+ ya.*ya + za.*za).^0.5;
    if LL==0
        error('LL is 0')
    end
        l=xa./LL;
        m=ya./LL;
    % % %     n=za./LL; 
    
    % if find(LL==0) > 0
    %     error('checkhere')
    % end
    
        D=sqrt(l.*l+m.*m);
    % % %     Lvx=-m./D;
    % % %     Lvy=l./D;
    % % %     Lvz=0./D;
    % % %     %Lvector2=[Lvx Lvy Lvx];
    % % %     
    % % %     Lv3x=-(l.*n./D);
    % % %     Lv3y=-(m.*n./D);
        
    CXx = (xa)./LL;
    CYx = (ya)./LL;
    CZx = (za)./LL;
    D = sqrt(CXx.*CXx + CYx.*CYx);
    CXy = -CYx./D;
    CYy = CXx./D;
    CZy = 0;
    CXz = -CXx.*CZx./D;
    CYz = -CYx.*CZx./D;
    CZz = D;
    
    
    
    EI=E*Ic(1); %% ONLY VALID WHEN ALL CNTS ARE IDENTICAL
    J=Ic(1);
    A=A(1);
    %K = spalloc(GDof,GDof,120*element);
    %K=sparse(GDof,GDof);
    
    Ni=6*(elementNodes(:,1)-1); %This represents a list of the first degree of freedom of the first node number of each frame element
    Nj=6*(elementNodes(:,2)-1); %This represents a list of the first degree of freedom of the second node number of each frame element
    
    ii_local=[
            Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6 Ni+1 Ni+2 Ni+3 Ni+4 Ni+5 Ni+6 Nj+1 Nj+2 Nj+3 Nj+4 Nj+5 Nj+6]';
     
    jj_local=[
        Ni+1 Ni+1 Ni+1 Ni+1 Ni+1 Ni+1 Ni+1 Ni+1 Ni+1 Ni+1 Ni+1 Ni+1 Ni+2 Ni+2 Ni+2 Ni+2 Ni+2 Ni+2 Ni+2 Ni+2 Ni+2 Ni+2 Ni+2 Ni+2 Ni+3 Ni+3 Ni+3 Ni+3 Ni+3 Ni+3 Ni+3 Ni+3 Ni+3 Ni+3 Ni+3 Ni+3 Ni+4 Ni+4 Ni+4 Ni+4 Ni+4 Ni+4 Ni+4 Ni+4 Ni+4 Ni+4 Ni+4 Ni+4 Ni+5 Ni+5 Ni+5 Ni+5 Ni+5 Ni+5 Ni+5 Ni+5 Ni+5 Ni+5 Ni+5 Ni+5 Ni+6 Ni+6 Ni+6 Ni+6 Ni+6 Ni+6 Ni+6 Ni+6 Ni+6 Ni+6 Ni+6 Ni+6 Nj+1 Nj+1 Nj+1 Nj+1 Nj+1 Nj+1 Nj+1 Nj+1 Nj+1 Nj+1 Nj+1 Nj+1 Nj+2 Nj+2 Nj+2 Nj+2 Nj+2 Nj+2 Nj+2 Nj+2 Nj+2 Nj+2 Nj+2 Nj+2 Nj+3 Nj+3 Nj+3 Nj+3 Nj+3 Nj+3 Nj+3 Nj+3 Nj+3 Nj+3 Nj+3 Nj+3 Nj+4 Nj+4 Nj+4 Nj+4 Nj+4 Nj+4 Nj+4 Nj+4 Nj+4 Nj+4 Nj+4 Nj+4 Nj+5 Nj+5 Nj+5 Nj+5 Nj+5 Nj+5 Nj+5 Nj+5 Nj+5 Nj+5 Nj+5 Nj+5 Nj+6 Nj+6 Nj+6 Nj+6 Nj+6 Nj+6 Nj+6 Nj+6 Nj+6 Nj+6 Nj+6 Nj+6]';
    
        
    y=zeros(144,element);
    
    y(1,:)=( (12*CXy.^2.*EI)./LL.^3 + (12.*CXz.^2.*EI)./LL.^3 + (A.*CXx.^2.*E)./LL )';
    y(13,:)=(  (12*CXy.*CYy.* EI)./LL.^3 + (12.*CXz.*CYz.*EI)./LL.^3 + (A.*CXx.*CYx.*E)./LL  )';
    y(25,:)=(  (12*CXy.*CZy.*EI)./LL.^3 + (12.*CXz.*CZz.*EI)./LL.^3 + (A.*CXx.*CZx.*E)./LL  )';
    %y(37,:)=0;
    y(49,:)=(  -((6*CXz.*CYy.*EI)./LL.^2) + (6*CXy.*CYz.*EI)./LL.^2  )';
    y(61,:)=(  -((6*CXz.*CZy.*EI)./LL.^2) + (6*CXy.*CZz.*EI)./LL.^2 )';
    y(73,:)=(  -((12*CXy.^2.*EI)./LL.^3) - (12*CXz.^2.*EI)./LL.^3 - (A.*CXx.^2.*E)./LL )';
    y(85,:)=(  -((12*CXy.*CYy.*EI)./LL.^3) - (12.*CXz.*CYz.*EI)./LL.^3 - (A.*CXx.*CYx.*E)./LL )';
    y(97,:)=(  -((12*CXy.*CZy.*EI)./LL.^3) - (12*CXz.*CZz.*EI)./LL.^3 - (A.*CXx.*CZx.*E)./LL )';
    %y(109,:)=0;
    y(121,:)=(  -((6*CXz.*CYy.*EI)./LL.^2) + (6*CXy.*CYz.*EI)./LL.^2 )';
    y(133,:)=(  -((6*CXz.*CZy.*EI)./LL.^2) + (6*CXy.*CZz.*EI)./LL.^2 )';
    
    y(2,:)=(  (12*CXy.*CYy.*EI)./LL.^3 + (12.*CXz.*CYz.*EI)./LL.^3 + (A.*CXx.*CYx.*E)./LL )';
    y(14,:)=(  (12*CYy.^2.*EI)./LL.^3 + (12.*CYz.^2.*EI)./LL.^3 + (A.*CYx.^2.*E)./LL )';
    y(26,:)=(  (12*CYy.*CZy.*EI)./LL.^3 + (12.*CYz.*CZz.*EI)./LL.^3 + (A.*CYx.*CZx.*E)./LL )';
    y(38,:)=(  (6*CXz.*CYy.*EI)./LL.^2 - (6*CXy.*CYz.*EI)./LL.^2 )';
    %y(50,:)=0;
    y(62,:)=( -((6.*CYz.*CZy.*EI)./LL.^2) + (6.*CYy.*CZz.*EI)./LL.^2 )';
    y(74,:)=( -((12.*CXy.*CYy.*EI)./LL.^3) - (12.*CXz.*CYz.*EI)./LL.^3 - (A.*CXx.*CYx.*E)./LL )';
    y(86,:)=( -((12.*CYy.^2.*EI)./LL.^3) - (12.*CYz.^2.*EI)./LL.^3 - (A.*CYx.^2.*E)./LL )';
    y(98,:)=( -((12.*CYy.*CZy.*EI)./LL.^3) - (12.*CYz.*CZz.*EI)./LL.^3 - (A.*CYx.*CZx.*E)./LL )';
    y(110,:)=( (6.*CXz.*CYy.*EI)./LL.^2 - (6.*CXy.*CYz.*EI)./LL.^2 )';
    %y(122,:)=0;
    y(134,:)=( -((6.*CYz.*CZy.*EI)./LL.^2) + (6.*CYy.*CZz.*EI)./LL.^2 )';
    
    
    y(3,:)=(  (12.*CXy.*CZy.*EI)./LL.^3 + (12.*CXz.*CZz.*EI)./LL.^3 + (A.*CXx.*CZx.*E)./LL)';
    y(15,:)=( (12.*CYy.*CZy.*EI)./LL.^3 + (12.*CYz.*CZz.*EI)./LL.^3 + (A.*CYx.*CZx.*E)./LL)';
    y(27,:)=( (12.*CZy.^2.*EI)./LL.^3 + (12.*CZz.^2.*EI)./LL.^3 + (A.*CZx.^2.*E)./LL)';
    y(39,:)=( (6.*CXz.*CZy.*EI)./LL.^2 - (6.*CXy.*CZz.*EI)./LL.^2)';
    y(51,:)=( (6.*CYz.*CZy.*EI)./LL.^2 - (6.*CYy.*CZz.*EI)./LL.^2)';
    %y(63,:)=0;
    y(75,:)=( -((12.*CXy.*CZy.*EI)./LL.^3) - (12.*CXz.*CZz.*EI)./LL.^3 - (A.*CXx.*CZx.*E)./LL)';
    y(87,:)=( -((12.*CYy.*CZy.*EI)./LL.^3) - (12.*CYz.*CZz.*EI)./LL.^3 - (A.*CYx.*CZx.*E)./LL)';
    y(99,:)=( -((12.*CZy.^2.*EI)./LL.^3) - (12.*CZz.^2.*EI)./LL.^3 - (A.*CZx.^2.*E)./LL)';
    y(111,:)=( (6.*CXz.*CZy.*EI)./LL.^2 - (6.*CXy.*CZz.*EI)./LL.^2)';
    y(123,:)=( (6.*CYz.*CZy.*EI)./LL.^2 - (6.*CYy.*CZz.*EI)./LL.^2)';
    %y(135,:)=0;
    
    
    %y(4,:)=0;
    y(16,:)=( (6.*CXz.*CYy.*EI)./LL.^2 - (6.*CXy.*CYz.*EI)./LL.^2)';
    y(28,:)=( (6.*CXz.*CZy.*EI)./LL.^2 - (6.*CXy.*CZz.*EI)./LL.^2)';
    y(40,:)=( (CXx.^2.*G.*J)./LL + (4.*CXy.^2.*EI)./LL + (4.*CXz.^2.*EI)./LL)';
    y(52,:)=( (CXx.*CYx.*G.*J)./LL + (4.*CXy.*CYy.*EI)./LL + (4.*CXz.*CYz.*EI)./LL)';
    y(64,:)=( (CXx.*CZx.*G.*J)./LL + (4.*CXy.*CZy.*EI)./LL + (4.*CXz.*CZz.*EI)./LL)';
    %y(76,:)=0;
    y(88,:)=( -((6.*CXz.*CYy.*EI)./LL.^2) + (6.*CXy.*CYz.*EI)./LL.^2)';
    y(100,:)=( -((6.*CXz.*CZy.*EI)./LL.^2) + (6.*CXy.*CZz.*EI)./LL.^2)';
    y(112,:)=( -((CXx.^2.*G.*J)./LL) + (2.*CXy.^2.*EI)./LL + (2.*CXz.^2.*EI)./LL)';
    y(124,:)=( -((CXx.*CYx.*G.*J)./LL) + (2.*CXy.*CYy.*EI)./LL + (2.*CXz.*CYz.*EI)./LL)';
    y(136,:)=( -((CXx.*CZx.*G.*J)./LL) + (2.*CXy.*CZy.*EI)./LL + (2.*CXz.*CZz.*EI)./LL)';
    
    y(5,:)=( -((6.*CXz.*CYy.*EI)./LL.^2) + (6.*CXy.*CYz.*EI)./LL.^2 )';
    %y(17,:)=0;
    y(29,:)=( (6.*CYz.*CZy.*EI)./LL.^2 - (6.*CYy.*CZz.*EI)./LL.^2 )';
    y(41,:)=( (CXx.*CYx.*G.*J)./LL + (4.*CXy.*CYy.*EI)./LL + (4.*CXz.*CYz.*EI)./LL )';
    y(53,:)=( (CYx.^2.*G.*J)./LL + (4.*CYy.^2.*EI)./LL + (4.*CYz.^2.*EI)./LL )';
    y(65,:)=( (CYx.*CZx.*G.*J)./LL + (4.*CYy.*CZy.*EI)./LL + (4.*CYz.*CZz.*EI)./LL )';
    y(77,:)=( (6.*CXz.*CYy.*EI)./LL.^2 - (6.*CXy.*CYz.*EI)./LL.^2 )';
    %y(89,:)=0;
    y(101,:)=( -((6.*CYz.*CZy.*EI)./LL.^2) + (6.*CYy.*CZz.*EI)./LL.^2 )';
    y(113,:)=( -((CXx.*CYx.*G.*J)./LL) + (2.*CXy.*CYy.*EI)./LL + (2.*CXz.*CYz.*EI)./LL )';
    y(125,:)=( -((CYx.^2.*G.*J)./LL) + (2.*CYy.^2.*EI)./LL + (2.*CYz.^2.*EI)./LL )';
    y(137,:)=( -((CYx.*CZx.*G.*J)./LL) + (2.*CYy.*CZy.*EI)./LL + (2.*CYz.*CZz.*EI)./LL )';
    
    
    y(6,:)=(  -((6.*CXz.*CZy.*EI)./LL.^2) + (6.*CXy.*CZz.*EI)./LL.^2 )';
    y(18,:)=(  -((6.*CYz.*CZy.*EI)./LL.^2) + (6.*CYy.*CZz.*EI)./LL.^2 )';
    %y(30,:)=0;
    y(42,:)=(  (CXx.*CZx.*G.*J)./LL + (4.*CXy.*CZy.*EI)./LL + (4.*CXz.*CZz.*EI)./LL )';
    y(54,:)=(  (CYx.*CZx.*G.*J)./LL + (4.*CYy.*CZy.*EI)./LL + (4.*CYz.*CZz.*EI)./LL )';
    y(66,:)=(  (CZx.^2.*G.*J)./LL + (4.*CZy.^2.*EI)./LL + (4.*CZz.^2.*EI)./LL )';
    y(78,:)=(  (6.*CXz.*CZy.*EI)./LL.^2 - (6.*CXy.*CZz.*EI)./LL.^2 )';
    y(90,:)=(  (6.*CYz.*CZy.*EI)./LL.^2 - (6.*CYy.*CZz.*EI)./LL.^2 )';
    %y(102,:)=0;
    y(114,:)=(  -((CXx.*CZx.*G.*J)./LL) + (2.*CXy.*CZy.*EI)./LL + (2.*CXz.*CZz.*EI)./LL )';
    y(126,:)=(  -((CYx.*CZx.*G.*J)./LL) + (2.*CYy.*CZy.*EI)./LL + (2.*CYz.*CZz.*EI)./LL )';
    y(138,:)=(  -((CZx.^2.*G.*J)./LL) + (2.*CZy.^2.*EI)./LL + (2.*CZz.^2.*EI)./LL )';
    
    y(7,:)=(  -((12.*CXy.^2.*EI)./LL.^3) - (12.*CXz.^2.*EI)./LL.^3 - (A.*CXx.^2.*E)./LL  )';
    y(19,:)=(  -((12.*CXy.*CYy.*EI)./LL.^3) - (12.*CXz.*CYz.*EI)./LL.^3 - (A.*CXx.*CYx.*E)./LL)';
    y(31,:)=(  -((12.*CXy.*CZy.*EI)./LL.^3) - (12.*CXz.*CZz.*EI)./LL.^3 - (A.*CXx.*CZx.*E)./LL)';
    %y(43,:)=0;
    y(55,:)=(  (6.*CXz.*CYy.*EI)./LL.^2 - (6.*CXy.*CYz.*EI)./LL.^2)';
    y(67,:)=(  (6.*CXz.*CZy.*EI)./LL.^2 - (6.*CXy.*CZz.*EI)./LL.^2)';
    y(79,:)=(  (12.*CXy.^2.*EI)./LL.^3 + (12.*CXz.^2.*EI)./LL.^3 + (A.*CXx.^2.*E)./LL)';
    y(91,:)=(  (12.*CXy.*CYy.*EI)./LL.^3 + (12.*CXz.*CYz.*EI)./LL.^3 + (A.*CXx.*CYx.*E)./LL)';
    y(103,:)=(  (12.*CXy.*CZy.*EI)./LL.^3 + (12.*CXz.*CZz.*EI)./LL.^3 + (A.*CXx.*CZx.*E)./LL)';
    %y(115,:)=0;
    y(127,:)=(  (6.*CXz.*CYy.*EI)./LL.^2 - (6.*CXy.*CYz.*EI)./LL.^2)';
    y(139,:)=(  (6.*CXz.*CZy.*EI)./LL.^2 - (6.*CXy.*CZz.*EI)./LL.^2)';
    
    
    y(8,:)=(  -((12.*CXy.*CYy.*EI)./LL.^3) - (12.*CXz.*CYz.*EI)./LL.^3 - (A.*CXx.*CYx.*E)./LL  )';
    y(20,:)=(  -((12.*CYy.^2.*EI)./LL.^3)- (12.*CYz.^2.*EI)./LL.^3 - (A.*CYx.^2.*E)./LL  )';
    y(32,:)=(  -((12.*CYy.*CZy.*EI)./LL.^3) - (12.*CYz.*CZz.*EI)./LL.^3 - (A.*CYx.*CZx.*E)./LL  )';
    y(44,:)=(  -((6.*CXz.*CYy.*EI)./LL.^2) + (6.*CXy.*CYz.*EI)./LL.^2  )';
    %y(56,:)=(  0  )';
    y(68,:)=(  (6.*CYz.*CZy.*EI)./LL.^2 - (6.*CYy.*CZz.*EI)./LL.^2  )';
    y(80,:)=(  (12.*CXy.*CYy.*EI)./LL.^3 + (12.*CXz.*CYz.*EI)./LL.^3 + (A.*CXx.*CYx.*E)./LL  )';
    y(92,:)=(  (12.*CYy.^2.*EI)./LL.^3 + (12.*CYz.^2.*EI)./LL.^3 + (A.*CYx.^2.*E)./LL  )';
    y(104,:)=(  (12.*CYy.*CZy.*EI)./LL.^3 + (12.*CYz.*CZz.*EI)./LL.^3 + (A.*CYx.*CZx.*E)./LL  )';
    y(116,:)=(  -((6.*CXz.*CYy.*EI)./LL.^2) + (6.*CXy.*CYz.*EI)./LL.^2  )';
    %y(128,:)=(  0  )';
    y(140,:)=(  (6.*CYz.*CZy.*EI)./LL.^2 - (6.*CYy.*CZz.*EI)./LL.^2  )';
    
    y(9,:)=(  -((12.*CXy.*CZy.*EI)./LL.^3) - (12.*CXz.*CZz.*EI)./LL.^3 - (A.*CXx.*CZx.*E)./LL  )';
    y(21,:)=(  -((12.*CYy.*CZy.*EI)./LL.^3) - (12.*CYz.*CZz.*EI)./LL.^3 - (A.*CYx.*CZx.*E)./LL  )';
    y(33,:)=(  -((12.*CZy.^2.*EI)./LL.^3) - (12.*CZz.^2.*EI)./LL.^3 - (A.*CZx.^2.*E)./LL  )';
    y(45,:)=(  -((6.*CXz.*CZy.*EI)./LL.^2) + (6.*CXy.*CZz.*EI)./LL.^2  )';
    y(57,:)=(  -((6.*CYz.*CZy.*EI)./LL.^2) + (6.*CYy.*CZz.*EI)./LL.^2  )';
    %y(69,:)=(  0  )';
    y(81,:)=(  (12.*CXy.*CZy.*EI)./LL.^3 + (12.*CXz.*CZz.*EI)./LL.^3 + (A.*CXx.*CZx.*E)./LL  )';
    y(93,:)=(  (12.*CYy.*CZy.*EI)./LL.^3 + (12.*CYz.*CZz.*EI)./LL.^3 + (A.*CYx.*CZx.*E)./LL  )';
    y(105,:)=(  (12.*CZy.^2.*EI)./LL.^3 +(12.*CZz.^2.*EI)./LL.^3 + (A.*CZx.^2.*E)./LL  )';
    y(117,:)=(  -((6.*CXz.*CZy.*EI)./LL.^2) + (6.*CXy.*CZz.*EI)./LL.^2  )';
    y(129,:)=(  -((6.*CYz.*CZy.*EI)./LL.^2) + (6.*CYy.*CZz.*EI)./LL.^2  )';
    %y(141,:)=(  0  )';
    
    
    %y(10,:)=(  0  )';
    y(22,:)=(  (6.*CXz.*CYy.*EI)./LL.^2 - (6.*CXy.*CYz.*EI)./LL.^2  )';
    y(34,:)=(  (6.*CXz.*CZy.*EI)./LL.^2 - (6.*CXy.*CZz.*EI)./LL.^2  )';
    y(46,:)=(  (2.*CXy.^2.*EI)./LL + (2.*CXz.^2.*EI)./LL - (CXx.^2.*G.*J)./LL  )';
    y(58,:)=(  (2.*CXy.*CYy.*EI)./LL + (2.*CXz.*CYz.*EI)./LL - (CXx.*CYx.*G.*J)./LL  )';
    y(70,:)=(  (2.*CXy.*CZy.*EI)./LL + (2.*CXz.*CZz.*EI)./LL - (CXx.*CZx.*G.*J)./LL  )';
    %y(82,:)=(  0  )';
    y(94,:)=(  -((6.*CXz.*CYy.*EI)./LL.^2) + (6.*CXy.*CYz.*EI)./LL.^2  )';
    y(106,:)=(  -((6.*CXz.*CZy.*EI)./LL.^2) + (6.*CXy.*CZz.*EI)./LL.^2  )';
    y(118,:)=(  (4.*CXy.^2.*EI)./LL + (4.*CXz.^2.*EI)./LL + (CXx.^2.*G.*J)./LL  )';
    y(130,:)=(  (4.*CXy.*CYy.*EI)./LL + (4.*CXz.*CYz.*EI)./LL + (CXx.*CYx.*G.*J)./LL  )';
    y(142,:)=(  (4.*CXy.*CZy.*EI)./LL + (4.*CXz.*CZz.*EI)./LL + (CXx.*CZx.*G.*J)./LL  )';
    
    y(11,:)=(  -((6.*CXz.*CYy.*EI)./LL.^2) + (6.*CXy.*CYz.*EI)./LL.^2  )';
    %y(23,:)=(  0  )';
    y(35,:)=(  (6.*CYz.*CZy.*EI)./LL.^2 - (6.*CYy.*CZz.*EI)./LL.^2  )';
    y(47,:)=(  (2.*CXy.*CYy.*EI)./LL + (2.*CXz.*CYz.*EI)./LL - (CXx.*CYx.*G.*J)./LL  )';
    y(59,:)=(  (2.*CYy.^2.*EI)./LL + (2.*CYz.^2.*EI)./LL- (CYx.^2.*G.*J)./LL  )';
    y(71,:)=(  (2.*CYy.*CZy.*EI)./LL + (2.*CYz.*CZz.*EI)./LL - (CYx.*CZx.*G.*J)./LL  )';
    y(83,:)=(  (6.*CXz.*CYy.*EI)./LL.^2 - (6.*CXy.*CYz.*EI)./LL.^2  )';
    %y(95,:)=(  0  )';
    y(107,:)=(  -((6.*CYz.*CZy.*EI)./LL.^2) + (6.*CYy.*CZz.*EI)./LL.^2  )';
    y(119,:)=(  (4.*CXy.*CYy.*EI)./LL + (4.*CXz.*CYz.*EI)./LL + (CXx.*CYx.*G.*J)./LL  )';
    y(131,:)=(  (4.*CYy.^2.*EI)./LL + (4.*CYz.^2.*EI)./LL + (CYx.^2.*G.*J)./LL  )';
    y(143,:)=(  (4.*CYy.*CZy.*EI)./LL + (4.*CYz.*CZz.*EI)./LL + (CYx.*CZx.*G.*J)./LL  )';
    
    y(12,:)=(  -((6.*CXz.*CZy.*EI)./LL.^2) + (6.*CXy.*CZz.*EI)./LL.^2  )';
    y(24,:)=(  -((6.*CYz.*CZy.*EI)./LL.^2) + (6.*CYy.*CZz.*EI)./LL.^2  )';
    %y(36,:)=(  0  )';
    y(48,:)=(  (2.*CXy.*CZy.*EI)./LL + (2.*CXz.*CZz.*EI)./LL - (CXx.*CZx.*G.*J)./LL  )';
    y(60,:)=(  (2.*CYy.*CZy.*EI)./LL + (2.*CYz.*CZz.*EI)./LL - (CYx.*CZx.*G.*J)./LL  )';
    y(72,:)=(  (2.*CZy.^2.*EI)./LL + (2.*CZz.^2.*EI)./LL - (CZx.^2.*G.*J)./LL  )';
    y(84,:)=(  (6.*CXz.*CZy.*EI)./LL.^2 - (6.*CXy.*CZz.*EI)./LL.^2  )';
    y(96,:)=(  (6.*CYz.*CZy.*EI)./LL.^2 - (6.*CYy.*CZz.*EI)./LL.^2  )';
    %y(108,:)=(  0  )';
    y(120,:)=(  (4.*CXy.*CZy.*EI)./LL + (4.*CXz.*CZz.*EI)./LL + (CXx.*CZx.*G.*J)./LL  )';
    y(132,:)=(  (4.*CYy.*CZy.*EI)./LL + (4.*CYz.*CZz.*EI)./LL + (CYx.*CZx.*G.*J)./LL  )';
    y(144,:)=(  (4.*CZy.^2.*EI)./LL + (4.*CZz.^2.*EI)./LL + (CZx.^2.*G.*J)./LL  )';
    
    
    
    %K=sparse(jj_local(:), ii_local(:), y(:),GDof,GDof);
    %K=sparse(jj_local, ii_local, y);
    
    K=fsparse(jj_local, ii_local, y,[GDof,GDof]);
    
    % % % K(jj_local(:), ii_local(:))=y(:);
    end
    