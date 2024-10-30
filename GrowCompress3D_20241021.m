function GrowCompress3D_20241021(steps,mu_ro,sigma_ro,shape_ri,scale_ri,avgRate,rate_stdev,ang_stdev,Pcritical,density_factor)
format long;format compact
addpath('/usr/src/MUMPS_5.7.3/MATLAB/')
addpath('/usr/src/MUMPS_5.7.3/MATLAB/stenglib/Fast')
warning('off','all')

E = 0.8e12; G = 5e11; % [Pa] CNT elastic modulus  

beamType = 0; % 0 = Euler Beam; 1 = Timoshenko Beam
compressiveLoad = 0;

% fname = '/home/maschmann/Documents/MATLAB/images/junk/'; %% File Path for Saving Results
fname = '/mnt/pixstor/data/glkxr9/';

% fname2 = strcat(fname,'TrackMatrixData\');
% fname3 = strcat(fname,'ResistanceData\');
% fname4 = strcat(fname,'ConductanceData\');

startTime2 = datetime('now', 'Format', 'yyyy_MM_dd_HH_mm_ss');
fname5 = strcat(fname,char(startTime2),'.csv');

% Load set parameter matrix
%paramat = load(strcat(fname,'SimuRun_Parameter_File_sorted.csv'));
% Pcritical = paramat(run,6)*1e6;                 % Pa
% density_factor = paramat(run,5);                % 1e10 CNT/cm2
% rate_stdev = paramat(run,1);                    % unitless
% h_span = paramat(run,2)*1e-6;                   % m
% v_span = paramat(run,3)*1e-6;                   % m
% electrodeWidth = paramat(run,4)*1e-6;           % m

% steps = 1000;
plotint = 50;

% Pcritical = 1.6e8;                 % Pa
% density_factor = 0.7;              % 1e10 CNT/cm2
% rate_stdev = 0.2;                  % This is now for a log-normal distribution
h_span = 5*1e-6;                   % m
v_span = 15*1e-6;                  % m
electrodeWidth = 0.5*1e-6;         % m

numberBeamsx = ceil(h_span*sqrt(density_factor)/1e-7);  % Span of the substrate lengthwise(meters)
numberBeamsy = ceil(v_span*sqrt(density_factor)/1e-7);  % Span of the substrate breadthwise(meters)
numberBeams = numberBeamsx*numberBeamsy; % Total number of Beams on the Substrate
            
ContinuousPlot = 0; % 0 = plotting off ; 1 = plotting on

element = numberBeams;
nodeCount = 2*numberBeams;

[ro,ri]=sample_radius(numberBeams,mu_ro,sigma_ro,shape_ri,scale_ri); % CNT radius [m]+

% avgRate = 60e-9;                  %nm/step
% ang_stdev = 5;                    % Degrees

totalCompress = 5;
PeriodicBoundary = 0;             % 0 = off ; 1 = on

%%%%%%%%%%%%%%%%%%%%%%%% Electrical setup info %%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vleft = 0; Vright = 0.1; %Voltages
conductivity = 3.3e8; % Electrical conductivity of 1 CNT ohm^-1 / meter
%electrodeWidth = 0.5e-6; % contact distance b/t electrodes

nuctitle = strcat(num2str(Pcritical/1e6),'MPA_dens',num2str(density_factor),'_rate',num2str(rate_stdev),'_hspan',num2str(h_span),'_vspan',num2str(v_span),'_elecW',num2str(electrodeWidth));  % Naming nucleation file

%%%%%%%%%%%%%%% Load previously saved nucleation data file %%%%%%%%%%%%%%%%
%load(strcat(fname,'NucleationFile_300MPA_Dens5_diam5_rate10_20x20CNT_30nmClose_distDiam','.mat'))
%load(strcat(fname,'NucleationFile__testfail_1_','.mat'))
% Pcritical = 100*1e6*6  
% title1 = char('_testfail_1_'); % File naming
% title = strcat(int2str(Pcritical/1e6),'MPA',title1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Nucleating CNTs with distributed properties according to inputs %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [elementNodes,ang,rate,nodeCoordinates,nodeCount,element,nucleationSite,growthNode,ro,rinn,phi]=nucleate3D_Curve(numberBeamsx,numberBeamsy,h_span,v_span,avgRate,rout,ri,rate_stdev,ang_stdev);
[elementNodes,ang,rate,nodeCoordinates,nodeCount,element,nucleationSite,growthNode,phi] = nucleate3D_Jitter(numberBeamsx,numberBeamsy,h_span,v_span,avgRate,rate_stdev,ang_stdev);


%load(strcat(fname,'NucInitData\',num2str(run),'_NucInitData_',nuctitle,'.mat'));

%%%%%%%%%%%%%%%%%%%%%% Save current nucleation data %%%%%%%%%%%%%%%%%%%%%%%%
%save(strcat(fname,'NucInitData\',num2str(run),'_NucInitData_',nuctitle,'.mat'));
A = pi*(ro.^2 - ri.^2);
Ic = pi/2*(ro.^4 - ri.^4)'; I = Ic;
EA = E.*A; EI = E.*Ic;
L = rate'; plotcount = 0; closeNodes = []; closeNodesOLD = []; oldCounter = 0;
removedNodes1 = []; removedNodes2 = [];
totalNodes = numberBeams*(steps+2); totalDOFs = numberBeams*(steps+2)*6;
inactiveNodes = [];
inactive_CNT_Time_Matrix = zeros(1,2);
basePress = ones(numberBeams,1);
oldDelam = [];
inactiveElements = [];
pushForce = zeros(6*numberBeams,1); oldForce = [];
check = 0;

startTime = datetime

for t = 1:steps % Growing CNT forest for a quantity of "steps" time steps
    looper=tic;
    % display current time step
    t
    % track calculation time: start
    tic;
    % global degrees of freedom (all system nodes)
    GDof = 6*nodeCount;  
     
    %     if closeNodes(1) == 0
    %         closeNodes = []; % closeNodes == [nodeCount,0]
    %         closeNodesOLD = [];
    %         oldCounter = 0;
    if closeNodes
        closeNodesOLD = closeNodes; % storing old closeNodes
        oldCounter = oldCounter + 1;
    end
   
    elapsedTime = toc(); % 'Section 2: Handling Close Nodes'
    fprintf('Section 2 runtime: %.5f seconds\n', elapsedTime);

    tic;
    if t > 1
        e = element - numberBeams + 1:element; L(e) = (rate(e-element+numberBeams)); A(e) = (A(e-element+numberBeams)); Ic(e) = Ic(e-element+numberBeams)';        
    end
    elapsedTime = toc; % Section 3: Updating Element Properties
    fprintf('Section 3 runtime: %.5f seconds\n', elapsedTime);

    tic;
    UTrack = zeros(GDof,50); TotalForceTrack = zeros(GDof,50);
    BeamForceTrack = zeros(GDof,50);
    PForceTrack = zeros(GDof,50); GeometricForce = zeros(element,6);
       
    entercloseNodes = 0;
    % tic()    
%     fastPath = 0
%
%     if fastPath == 0; [closeNodes,vox,uniquevox] = FindCloseNodes_RangeDec22_3D(nodeCoordinates,nodeCount); end
%
%     if fastPath == 1; [closeNodes] = FindCloseNodes_Voxel_Par(nodeCoordinates,nodeCount); end
%     if rem(t,20)==0
%     % determine nodes close enough for vdw forces
% %     if t < 1000
% tic()
% if fastPath
%         % [closeNodes,vox,uniquevox] = FindCloseNodes_RangeDec22_3D(nodeCoordinates,nodeCount);
%         single=toc()
%         tic();
%         [closeNodes] = FindCloseNodes_Voxel_Par(nodeCoordinates,nodeCount);
%         double = toc();
%     else
%          [closeNodes] = FindCloseNodes_Voxel_Par(nodeCoordinates,nodeCount);
%     end
    leaveCloseNodes = 1;
    % track calculation time: stop
    % close = toc()
    % remove repeated closeNodes rows (must also be from bottom up)

    closeNodes = [];

    % use average ro value for gap 
    r_mu = 12.34;
    gap = 2*r_mu+20e-9;


    [Idx,D] = knnsearch(nodeCoordinates,nodeCoordinates,"K",3);
    [C1,C2] = find(D>0 & D<gap);
    closeNodes = [Idx(C1,1), Idx(C1,2)];
    sizeClose=size(closeNodes)

%%%%%
    if size(closeNodes,1) == 0
        closeNodes = [0,0];
    end
    closeNodes = unique(closeNodes,'rows');

    if oldCounter > 1
        closeNodes = [closeNodes;closeNodesOLD];
    end
    closeNodes = unique(closeNodes,'rows','stable');

    elapsedTime = toc(); % Section 4: Finding Close Nodes
    fprintf('Section 4 runtime: %.5f seconds\n', elapsedTime);
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Calculate the delaminated CNTs %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic;
    delamCNTs = find(basePress > Pcritical); % This will be a list of CNTs
    newDelam = setdiff(delamCNTs,oldDelam); % List of new delaminated CNTs
    numNew = length(newDelam);
    inactiveCNTs = [[newDelam];[oldDelam]];
    numberInactive = length(inactiveCNTs); % Display

    if newDelam
        for i = 1:numNew
            inactive_CNT_Time_Matrix(length(oldDelam) + i, 1) = newDelam(i);
            inactive_CNT_Time_Matrix(length(oldDelam) + i, 2) = t;
            inactive_CNT_Time_Matrix(length(oldDelam) + i, 3) = basePress(newDelam(i));
            inactive_CNT_Time_Matrix(length(oldDelam) + i, 4) = mean(nodeCoordinates(:,3));
            inactive_CNT_Time_Matrix(length(oldDelam) + i, 5) = size(closeNodes,1);
        end
    end
    %%% Remove delaminated closeNodes (no force between active & delaminated CNTs)
    inactiveElements=[];
    if numberInactive > 0
        % for each delaminated CNT
% % % % %         for i = 1:size(inactive_CNT_Time_Matrix,1)
% % % % %             start = numberBeams*(inactive_CNT_Time_Matrix(i,2)-1) + inactive_CNT_Time_Matrix(i,1); % All free nodes for the inactiveCNTs
% % % % %             inactives = (start:numberBeams:element)'; % All inactive nodes
% % % % %             inactiveElements = [inactives;inactiveElements]; %
% % % % %         end
% % % % %         inactiveElements=zeros(size(inactive_CNT_Time_Matrix,1),t);
% % % % %         inactiveElements();
        counter=0;
        for i = 1:size(inactive_CNT_Time_Matrix,1)
            counter = counter+1;
            starts(counter) = numberBeams*(inactive_CNT_Time_Matrix(i,2)-1) + inactive_CNT_Time_Matrix(i,1); % All free nodes for the inactiveCNTs
        end
        num = 1:size(starts,2);
        inactive=zeros(size(starts,2),t);
        for num = 1:size(starts,2);
        p(num) = size(starts(num):numberBeams:element , 2);
        inactive(num,1:p(num))=starts(num):numberBeams:element;
        end
        inactiveElements=nonzeros(inactive);
   
        inactiveNodePairs = elementNodes(inactiveElements,:); % Nodal pair for inactive element
        inactiveNodes = inactiveNodePairs(:); % List of inactive nodes
        [val,ia,ib] = intersect(closeNodes(:,1),inactiveNodes);
        [val,iia,iib] = intersect(closeNodes(:,2),inactiveNodes);
        [inactiveConnect] = [ia;iia];
        % Clear vdwf for delaminated nodes that are close
        if length(inactiveConnect) > 0
            closeNodes(inactiveConnect,:) = [];
        end
    end
    oldDelam = unique(inactiveCNTs);
    trackDelam(t) = length(oldDelam);
    elapsedTime = toc();  % Section 5: Handling Delaminated CNTs
    fprintf('Section 5 runtime: %.5f seconds\n', elapsedTime);
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic;
    buildK = 0;
    [K,LL] = SpaceBeamElemStiffness_mod4(E,A,GDof,Ic,G,elementNodes,element,nodeCoordinates,nodeCount,PeriodicBoundary,h_span,v_span);
    buildK = 1;
    K_Elect = PlaneFrameElec(A,L,elementNodes,GDof,conductivity,nodeCount,numberBeams);
    elapsedTime = toc();  % Section 6: Updating Stiffness Matrix
    fprintf('Section 6 runtime: %.5f seconds\n', elapsedTime);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% ATTRACTIVE FORCES BY VAN DER WAALS ATTRACTION %%%%%%%%%%%%%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic;
    if closeNodes
        [Ac,vdwk] = ConnectionStiffness_3D_mod3(closeNodes,GDof,beamType,nodeCoordinates);
        K = K + Ac;
    end

    if closeNodes
    [Aelect] = ConnectionElect(closeNodes,nodeCount); %Contact Resistance
    K_Elect = K_Elect + Aelect;
    end
    elapsedTime = toc(); % Section 7: Calculating van der Waals Forces
    fprintf('Section 7 runtime: %.5f seconds\n', elapsedTime);

    tic;
    entersolution = 0;
    % U is displacement vector. oldForce is force vector pushing up to
    % displace the bottom nodes
    % calculate all displacements (per DOF) & bottom node forces (per DOF)

    [U,Fpush] = solution_inactivesAN(GDof,K,numberBeams,E,A,ang,phi,rate,oldForce,inactiveCNTs,t);

    forceTrack(t,:) = Fpush;
    % determine new force based on previous and additional force
    pushForce = pushForce + Fpush;
    % update active CNTs by removing delaminated cnts (based on inactive cnt list)
    activeCNTs = sort(setdiff(1:numberBeams,inactive_CNT_Time_Matrix(:,1)));
    % calculate the vertical force exerted on active cnts
    vertForce(activeCNTs) = pushForce(6*(activeCNTs-1)+3)./cos(ang(activeCNTs))';
    % update base pressure on active cnts from new vertical force
    basePress(activeCNTs) = vertForce(activeCNTs)./A(activeCNTs)';
    %basePress = pushForce(2:6:6*numberBeams)./A(1:numberBeams); The highlighted text appears to be y-force instead of z-force. Commented the entire line
    leavesolution = 1;
    elapsedTime = toc(); % Section 8: Solving Displacements
    fprintf('Section 8 runtime: %.5f seconds\n', elapsedTime);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Electrical SOLUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic;
    i_touch = (nodeCount-numberBeams+activeCNTs); centerContacts=[];rightContacts=[];leftContacts=[];
    contacts = size(i_touch,2);
    miny = min(nodeCoordinates(:,2));
    V(nodeCount+1) = Vleft;
    %nodeCoordinates(nodeCount+1,:)=[electrodeWidth/2 miny];
    V(nodeCount+2) = Vright;
    %nodeCoordinates(nodeCount+2,:)=[span-electrodeWidth/2 miny];
    %V(nodeCount+3)=Vright;
   
    nodecountleft=0; nodecountright=0; nodecountcenter=0; electContact=0;% leftContacts(1)=numberBeams/2+1; rightContacts(1)=numberBeams-5;
    % GG EEEEE GGGGG EEEEE GGG
    for i = 1:contacts
        if nodeCoordinates(i_touch(i),1) < electrodeWidth
            nodecountleft = nodecountleft + 1;
            leftContacts(nodecountleft) = i_touch(i);
        end
       
        if nodeCoordinates(i_touch(i),1) > (h_span - electrodeWidth)
            nodecountright = nodecountright + 1;
            rightContacts(nodecountright) = i_touch(i);
        end
    end

    nodecountleft;
    nodecountright;

    % if nodecountleft > 0 && nodecountcenter > 0 || nodecountright > 0 && nodecountcenter > 0
    if nodecountleft > 0 && nodecountright > 0
        %electContact = 1
        % For periodic boundary conditions, need to remove elements that cross
        % across boundary
        % crossedElements = [crossedl; crossedr];
       
        %%% First solving for the case of contact resistance
        [Aleft] = ConnectionElect_ElectrodesL(leftContacts,nodeCount);
        [Aright] = ConnectionElect_ElectrodesR(rightContacts,nodeCount);
        K_Elect = K_Elect + Aleft + Aright;
        contactsCount(t) = size(rightContacts,2) + size(leftContacts,2);
        voltageNodes = [nodeCount+1 nodeCount+2]; % nodeCount+3];
       
        prescribedDof_Elect = [voltageNodes]; % [leftContacts rightContacts centerContacts voltageNodes]';
        % activeDofElect = setdiff(1:nodeCount);%+3,prescribedDof_Elect);
        V = solutionElect(GDof,prescribedDof_Elect,K_Elect,V,nodeCount,numberBeams);
               currents = K_Elect*V';

        if currents >0
            electContact = 1
        end
        totalcurrent(t) = sum(abs(currents(leftContacts)));
        totalCurrentRight(t) = sum(abs(currents(rightContacts)));
        Resistance(t) = ((Vright-Vleft)/currents(nodeCount+2));
        res = Resistance(t);
    end

    %      %i=zeros(GDof/3,1);
    %      %m=min(nodeCoordinates(:,2)); %Nodes contacting bottom Surf
    %      i_touch=find(nodeCoordinates(:,2) < (min(nodeCoordinates(:,2))+2*10^-9))
    %      contacts=size(i_touch,1);
    %      V(nodeCount+1) = Vleft;
    %      V(nodeCount+2) = Vright;
    %
    %      nodecountleft=0; nodecountright=0; electContact=0;% leftContacts(1)=numberBeams/2+1; rightContacts(1)=numberBeams-5;
    %      for i=1:contacts
    %          if nodeCoordinates(i_touch(i),1)<electrodeWidth
    %              nodecountleft=nodecountleft+1;
    %              leftContacts(nodecountleft)=i_touch(i);
    %          end
    %          if nodeCoordinates(i_touch(i),1)> (h_span-electrodeWidth)
    %              nodecountright=nodecountright+1;
    %              rightContacts(nodecountright)=i_touch(i);
    %          end
    %      end
    %          
    %      if nodecountleft>0 && nodecountright>0
    %      electContact=1;
    %      V=zeros(nodeCount,1);  
    %  
    %      prescribedDof_Elect=[leftContacts rightContacts]';
    %      [V]=solutionElect(GDof,prescribedDof_Elect,K_Elect,V,nodeCount,numberBeams)
    %      %V=V+Vactive;
    %      
    %      currents=K_Elect*V;
    %      totalcurrent(t)=sum(currents(rightContacts));
    %      Resistance(t)=(abs(Vright-Vleft))/totalcurrent(compress)
    %      %currentmax=max(current)
    %      end

    elapsedTime = toc(); % Section 9: Solving Electical
    fprintf('Section 9 runtime: %.5f seconds\n', elapsedTime);

    %%%%%%%%%%%%%%%%%%%%%%%%% REPOSITION NODES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate new node positions based on displacements & add new nodes at bottom/surface
    tic;
    [nodeCoordinates,nodeCount,element,elementNodes] = reposition_inactiveAN(nodeCoordinates,U,nodeCount,numberBeams,growthNode,element,nucleationSite,elementNodes,PeriodicBoundary,h_span,v_span,inactiveCNTs,rate,ang,phi,LL,t);
    closeNodesOLD = closeNodes;
    elapsedTime = toc(); % Section 10: Repositioning Nodes
    fprintf('Section 10 runtime: %.5f seconds\n', elapsedTime);

    tic;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% RESISTANCE PLOT  %%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  if electContact>0 && rem(compress,5)==0;
    %        plot(compressDisp*(1:compress)*10^6,abs(Resistance),'o', 'Color','red')
    %        axis([0 10 1 1000]);
    %        xlabel('Compression (\mum)');
    %        ylabel('Resistance (Ohms)');
    %        name=char('Resistance');
    %        plotname=strcat(title,name,num2str(t));
    %        hold on
    %       % saveas(gcf,fullfile(fname,plotname),'png');
    %         %close()
    %  end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%     if electContact > 0 && rem(t,3000) == 0
%        %CNTPlotFastCompress(fname,avgRate, steps, span, numberBeams,ContinuousPlot,nodeCount,nodeCoordinates,t,title,totalCompress, compressiveLoad,tgrowth)
%        CNTPlot_Elect(fname, h_span, numberBeams,nodeCount,nodeCoordinates,t,title, V, elementNodes,Resistance)
%        %CNTPlotCompress(fname,avgRate, steps, span, numberBeams,ContinuousPlot,nodeCount,nodeCoordinates,t,title,totalCompress, compressiveLoad,bottomNode)
%     end 
    if rem(t,plotint) == 0
%         if numberInactive
            CNTPlotBWFast(fname,avgRate,steps,h_span,v_span,numberBeams,ContinuousPlot,nodeCount,nodeCoordinates,t,totalCompress,compressiveLoad,numberBeamsx,numberBeamsy,inactiveNodes)
%         end
        %saveme = strcat(fname,title,'plot','.png');
    end
    if t > 1 && sum(inactive_CNT_Time_Matrix(:,1)) > 0
        trackmatrix = [];
        % trackmatrix= zeros(length(inactive_CNT_Time_Matrix),6);
        trackmatrix(:,1) = inactive_CNT_Time_Matrix(:,1); % CNT that delaminated
        trackmatrix(:,2) = inactive_CNT_Time_Matrix(:,2); % time of delamination
        trackmatrix(:,3) = ro(1);
        trackmatrix(:,4) = ri(1);
        trackmatrix(:,5) = inactive_CNT_Time_Matrix(:,3); % delamination critical pressure     
        trackmatrix(:,6) = rate(inactive_CNT_Time_Matrix(:,1));
        trackmatrix(:,7) = inactive_CNT_Time_Matrix(:,4); % forest height at time step
        trackmatrix(:,8) = inactive_CNT_Time_Matrix(:,5); % number of close nodes at time step
        % csvwrite(strcat(fname2,num2str(run),'dlam_','Pcrit_',num2str(Pcritical/1e6),'_CNTdens=',num2str(numberBeamsx*numberBeamsy)',num2str(rout/1e-9),'rate=',num2str(rate_stdev),'imaging','.csv'),trackmatrix);
    end

    if nodecountleft > 0 && nodecountright > 0
        resdata(t) = Resistance(t);
        condata = 1./resdata;
        % csvwrite(strcat(fname3,num2str(run),'ResistanceData)_imaging','.csv'),resdata');
        % csvwrite(strcat(fname4,num2str(run),'ConductanceData_imaging','.csv'),condata');
    end

    numberActive=numberBeams-numberInactive

    height(t)=mean(nodeCoordinates(:,3))*1e6;
    density(t)=double(numberActive)/double(h_span*v_span*1e12);
    
    

    
    elapsedTime = toc(); % Section 11: Remainder
    fprintf('Section 11 runtime: %.5f seconds\n', elapsedTime);

    elapsedTime = toc(looper); % Entire step
    fprintf('Total step runtime: %.5f seconds\n', elapsedTime);
end

% Write the units
%writematrix(['steps (int), ','E (10^11 Pa), ','G (10^11 Pa), ','avgRate (nm/step), ','rate_stdev (nm/step), ','ang_stdev (deg), ','Pcritical (MPa), ','ri (nm), ','ro (nm), '],fname5)

% Prepare the first line with actual variable values
firstLine = [steps, E*1e-11, G*1e-11, avgRate*1e9, rate_stdev, ang_stdev, Pcritical*1e-6];

% Write the first line
writematrix(firstLine, fname5,'WriteMode','append','Delimiter', ',');

% Write the units
%writematrix('HEIGHT (nm)',fname5,'WriteMode','append','Delimiter',',');
% Append the data
writematrix(height, fname5,'WriteMode','append','Delimiter',',');

% Write the units
%writematrix('RESISTANCE (Ohms? - need to check)',fname5,'WriteMode','append','Delimiter',',');
% Append the data
writematrix(Resistance,fname5,'WriteMode','append','Delimiter',',');

% Write the units
%writematrix('DENSITY (#active CNTs / sq micron)',fname5,'WriteMode','append','Delimiter',',')
% Append the data
writematrix(density,fname5,'WriteMode','append','Delimiter',',');

stopTime = datetime
ending = toc()
% plot CNT growth (*check plotting interval, plotint)
% CNTPlotBW(fname,avgRate,steps,h_span,v_span,numberBeams,ContinuousPlot,nodeCount,nodeCoordinates,t,title,totalCompress,compressiveLoad,numberBeamsx,numberBeamsy,inactiveNodes)
   
GDof = 6*size(nodeCoordinates,1);
UU = zeros(GDof,1); forceCompress = zeros(GDof,1); U = zeros(GDof,1); compressiveLoad = zeros(totalCompress,1); force = zeros(GDof,1);
reactionForce = zeros(GDof,1);

for e = element - numberBeams + 1:element
    L(e) = (rate(e-element+numberBeams));
    A(e) = (A(e-element+numberBeams));
    I(e) = I(e-element+numberBeams)';
end

topNode = max(nodeCoordinates(:,3)); % Top-most node in forest at initial time step
closeNodesOld = closeNodes; % final update on all close nodes
Fbot = zeros(GDof,1); Ftop = zeros(GDof,1);
