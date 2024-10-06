E=0.7e12;               %%% CNT elastic modulus
warning('off','all')

rout=5e-9;              %%% CNT outer radius
rate_stdev=10;			%%% CNT population growth rate std. deviation [%]
ang_stdev=5;            %%% CNT population orientation angle std. deviation [%]
avgRate=60e-9;          %%% Average growth per time step (meters)

Pcritical = 200*1e6 %Delamination stress, Pa

ri=0.8*rout; A=pi*(rout^2-ri^2); Ic=pi/2*(rout^4-ri^4); EA=E*A; EI=E*Ic; G=5e11;

title1=char('_Dens03_diam10_rate5_225CNT_');
title=strcat(int2str(Pcritical/1e6),'MPa',title1);   %%% File name for images

fname = 'C:\Users\maschmannm\Documents\MATLAB\Image Folder\Delam\100CNTStudy\Dense03\';
fname2 = 'C:\Users\maschmannm\Documents\MATLAB\Image Folder\Delam\100CNTStudy\Dense03\';
plotint=500;            %%% Plotting incriment

% % % CNT_num = density_factor*1e2;
% % % cnt_num = sqrt(1/CNT_num);
% % % maxNumCompThreads('automatic')

steps = 3000;% Number of time steps for full CNT forest consideration
numberBeamsx=ceil(15);% Number of CNTs to be modeled lengthwise
numberBeamsy=ceil(15);% Number of CNTs to be modeled breadthwise
cnt_num = numberBeamsx*numberBeamsy;

% % h_span=1/sqrt(density_factor)*numberBeamsx*1e-7;% Span of the substrate lengthwise(meters)
% % v_span=1/sqrt(density_factor)*numberBeamsy*1e-7;% Span of the substrate breadthwise(meters)

segmentTime=zeros(steps,1);

numberBeams=numberBeamsx*numberBeamsy;% Total number of Beams on the Substrate
ContinuousPlot=0; 
% 0=plotting off.  1=plotting on
PeriodicBoundary=1;% 0=off, 1=on
beamType=0;% 0=Euler Beam ; 1=Timoshenko Beam
totalCompress=0;
compressiveLoad=0;
element=numberBeams;
nodeCount=2*numberBeams;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%/ nucleating CNTs with distributed properties according to inputs%%%%/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[elementNodes,ang,rate,nodeCoordinates,nodeCount,element,nucleationSite,growthNode,ro,rinn,phi]=nucleate3D_Curve(numberBeamsx,numberBeamsy,h_span,v_span,avgRate,rout,ri,rate_stdev,ang_stdev);

A=pi*(ro.^2-ri.^2)';
Ic=pi/2*(ro.^4-ri.^4)'; 
EA=E.*A; EI=E.*Ic;
L=rate'; plotcount=0;closeNodes=[0 0]; closeNodesOLD=[0 0];

removedNodes1=[]; removedNodes2=[];
totalNodes=numberBeams*(steps+2); totalDOFs=numberBeams*(steps+2)*6;
inactiveNodes=[];
inactive_CNT_Time_Matrix=zeros(1,2);
basePress=ones(numberBeams,1);
oldDelam=[];
inactiveElements=[];

for t=1:steps %% Growing CNT forest for a quantity of "steps" time steps
 
    tic
    t                   %%% Display time step to screen
    GDof=6*nodeCount;  
   
    if closeNodes(1)==0;
        closeNodes=[0 0]; %%[nodeCount,0] closeNodes==[nodeCount,0]
        closeNodesOLD=[0 0];
        oldCounter=0;
    else
        closeNodesOLD=closeNodes; %storing old closeNodes
        oldCounter=oldCounter+1;
    end
     
   
    if t>1
            e=element-numberBeams+1:element; L(e)=(rate(e-element+numberBeams));   A(e)=(A(e-element+numberBeams));  Ic(e)=Ic(e-element+numberBeams)';        
    end
    
    UTrack=zeros(GDof,50); TotalForceTrack=zeros(GDof,50);
    BeamForceTrack=zeros(GDof,50);
    PForceTrack=zeros(GDof,50); GeometricForce=zeros(element,6);
        
    entercloseNodes=0;
         [closeNodes]=FindCloseNodes_Range(nodeCoordinates,nodeCount);

         if oldCounter>1
            %closeNodesCombined=vertcat(closeNodes,closeNodesOLD);
            [closeNodesCombined]=[closeNodes;closeNodesOLD];
            closeNodes=unique(closeNodesCombined,'rows');
         end
         
         size(closeNodes);
    leaveCloseNodes=1;
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Calculate the delaminated CNTs %%%%%%%%%%%%%%%%%%%%
   delamCNTs = find(basePress>Pcritical);%%% This will be a list of CNTs
   
   newDelam = setdiff(delamCNTs,oldDelam);  %%% List of new delaminated CNTs
   
   
   
   numNew = length(newDelam);
   inactiveCNTs = [[newDelam];[oldDelam]];
   numberInactive = length(inactiveCNTs) %Display
   
    if newDelam
        for i = 1:numNew 
            inactive_CNT_Time_Matrix(length(oldDelam) + i, 1) = newDelam(i);
            inactive_CNT_Time_Matrix(length(oldDelam) + i, 2) = t;
            inactive_CNT_Time_Matrix(length(oldDelam) + i, 3) = basePress(newDelam(i));
            inactive_CNT_Time_Matrix(length(oldDelam) + i, 4) = max(nodeCoordinates(:,3));
            
        end
    end
   
if numberInactive > 0
    
   for i=1:size(inactive_CNT_Time_Matrix,1)
       start=numberBeams*(inactive_CNT_Time_Matrix(i,2)-1)+inactive_CNT_Time_Matrix(i,1);
       inactives=(start:numberBeams:element)'; %All inactive nodes
       %inactiveElements=vertcat(inactives,inactiveElements);
       [inactiveElements]=[inactives;inactiveElements];
   end
       inactiveNodePairs=elementNodes(inactiveElements,:);
       inactiveNodes=inactiveNodePairs(:);   
                  
       [val,ia,ib] = intersect(closeNodes(:,1),inactiveNodes);
       [val,iia,iib] = intersect(closeNodes(:,2),inactiveNodes);
       %inactiveConnect=vertcat(ia,iia);
       [inactiveConnect]=[ia;iia];
       
       if length(inactiveConnect)>0
          closeNodes(inactiveConnect,:)=[];
       end

end
   oldDelam = unique(inactiveCNTs);
  %basePressDelam = basePress(oldDelam)
   trackDelam(t)=length(oldDelam);
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

buildK=0;
K= SpaceBeamElemStiffness_mod4(E,A,GDof,Ic,G,elementNodes,element,nodeCoordinates,nodeCount,PeriodicBoundary,h_span,v_span)  ;
buildK=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%/ ATTRACIVE FORCES BY VAN DER WAALS ATTRACTION %%%%%%%%%%%%%%%%%%%%%%%%%%/    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/
    

if size(closeNodes,1)>1  
[Ac,vdwk]=ConnectionStiffness_3D_mod3(rout,closeNodes,GDof,beamType,nodeCoordinates);
%Ac = distributed(Ac);


K=K+Ac; 

end

entersolution=0

% % % if t < 11
% % %     U = solution_mod3 (GDof,K,numberBeams,E,A,ang,phi);
% % % end
% % % if t > 10
% % %     U = solution_mod5 (GDof,K,numberBeams,E,A,ang,phi,Uo);
% % % end


    
    [U] = solution_mod4(GDof,K,numberBeams,E,A,ang,phi);
    %U=solution_mod3(GDof,KK,numberBeams,E,A,ang,phi);


    
    %Uo = U;
    
    
    %alpha = 0.1;
    %L = ichol(K,struct('type','ict','droptol',1e-2,'diagcomp',alpha));
   % M1 = L; M2 = L'; L_T=M2; U=U0;
    
    %for ii=1:NL
    %  [U(:)] = pcg(K,F(:,ii),tol_pcg,50000,M1,M2,U0); % Fastest
    %end
    
    
    
    %U = U0;
    %[L_T,~,s] = chol(K,'lower','vector');
    %U(s,:) = L_T'\(L_T\F(s,:)); 
    
    
    
leavesolution=1



   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   %%%%%%%%%%%%% Calculate Force of Bottom Elements %%%%%%%%%%%%%%%%%%%% 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if t>5
        it=zeros(1,numberBeams);
        it=element-2*numberBeams+1: element-numberBeams;
        evalNodes=it;
        %bottomDofs=GDof-6*numberBeams-(3*(numberBeams-1)):GDof-6*numberBeams;
       
      
        bottomDofs=GDof-12*numberBeams-(6*(numberBeams-1)):GDof-12*numberBeams;
                 
        %Kk = StiffnessForceElements(EE,AAA,II,Ldef,C,S,GG,evalNodes,elementNodes,GDof,beamType,totalDOFs);%stiffnessMatrix of just bottom elements
        Kk = SpaceBeamElemStiffness_Force(E,A,GDof,Ic,G,evalNodes,elementNodes,element,nodeCoordinates,nodeCount,PeriodicBoundary,h_span,v_span); 
        bbaseForce=Kk*U;
        baseForce=bbaseForce(6*(evalNodes-1)+3);  %Vector of force
        basePress=baseForce./A(1:numberBeams);    %Vector of pressure
        botK=zeros(numberBeams,12); botU=zeros(numberBeams,12);

        %plot(baseForce,'.')
        % baseForce(timestep,:)=-K(vdofOrig,:)*U-R(vdofOrig); %vertical force at base
       % maxBase(timestep)=max(baseForce(timestep,:)); minBase(timestep)=min(baseForce(timestep,:));
   end
    
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[nodeCoordinates,nodeCount,element,elementNodes]=reposition(nodeCoordinates,U,nodeCount,numberBeams,growthNode,element,nucleationSite,elementNodes,PeriodicBoundary,h_span,v_span);
segmentTime(t)=toc;
segmentTime(t)

closeNodesOLD=closeNodes;




    if rem(t,plotint)==0 
    %save([fullfile(fname,title),'.mat'])
    %CNTPlotFast(fname,avgRate, steps, numberBeams,ContinuousPlot,nodeCount,nodeCoordinates,t,title,totalCompress, compressiveLoad)
    
    CNTPlotBW(fname,avgRate, steps, h_span,v_span, numberBeams,ContinuousPlot,nodeCount,nodeCoordinates,t,title,totalCompress, compressiveLoad,numberBeamsx,numberBeamsy,inactiveNodes)
    saveme=strcat(fname,title,'.mat');
 %save(saveme,'ang_stdev','rate_stdev','avgRate','rate','inactive_CNT_Time_Matrix','ro','h_span','v_span','numberBeamsx','numberBeamsy','Pcritical');
%save(saveme,'K');
    end

if t>100 && sum(inactive_CNT_Time_Matrix(:,1))>0 && mod(t,10)==0

  trackmatrix= zeros(length(inactive_CNT_Time_Matrix),6);
  trackmatrix(:,1) = inactive_CNT_Time_Matrix(:,1);
  trackmatrix(:,2) = inactive_CNT_Time_Matrix(:,2);
  trackmatrix(:,3) = ro(1);
  trackmatrix(:,4) = ri(1);
    
  trackmatrix(:,6) = rate(inactive_CNT_Time_Matrix(:,1));
  trackmatrix(:,5) = inactive_CNT_Time_Matrix(:,3);
  trackmatrix(:,7) = inactive_CNT_Time_Matrix(:,4);
  csvwrite(strcat(fname2,'Test_', int2str(number),'dlam_','_Pcrit_',num2str(Pcritical/1e6),'_CNTdens=',num2str(numberBeamsx*numberBeamsy),'dia=',num2str(rout/1e-9),'rate=',num2str(rate_stdev),'.csv'),trackmatrix);
end

end




