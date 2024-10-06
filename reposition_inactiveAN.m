function [nodeCoordinates,nodeCount,element,elementNodes] = reposition_inactiveAN(nodeCoordinates,U,nodeCount,numberBeams,growthNode,element,nucleationSite,elementNodes,PeriodicBoundary,h_span,v_span,inactiveCNTs,rate,ang,phi,LL,t)
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% REPOSITIONING NODES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % define active cnts
    activeCNTs = 1:numberBeams;
    
    % update active cnts
    if inactiveCNTs
        activeCNTs = setdiff(1:numberBeams,inactiveCNTs); % may be over written if inactiveCNTs are present
    end

%% Shift all upper nodes to their new positions based on {U}. Not affected by delaminated CNTs
    
    % reference to all upper nodes (no bottom nodes)
    ii = 1:nodeCount - numberBeams;
    
    % displace all upper nodes
    nodeCoordinates(ii,1) = nodeCoordinates(ii,1) + U(6*(ii-1)+1); % displaced nodal locations in x
    nodeCoordinates(ii,2) = nodeCoordinates(ii,2) + U(6*(ii-1)+2); % displaced nodal locations in y
    nodeCoordinates(ii,3) = nodeCoordinates(ii,3) + U(6*(ii-1)+3); % displaced nodal locations in z
       
%% Moving the bottom nodes of inactive CNTs - they are positioned in a straight line relative to the last CNT
    
    if inactiveCNTs  
        % reference to bottom nodes of inactive cnts
        nn = nodeCount - numberBeams + inactiveCNTs;
     
        % check for identical updated locations
        [Bb,~] = size(nodeCoordinates);
        [B,~] = size(unique(nodeCoordinates,'rows'));
        if Bb > B
            xk1 = find(setdiff(Bb,B))
            error('1')
        end 

        % displace bottom nodes of inactive cnts
        nodeCoordinates(nn,1) = nodeCoordinates(nn,1) + U(6*(nn-1)+1); % displaced nodal locations in x
        nodeCoordinates(nn,2) = nodeCoordinates(nn,2) + U(6*(nn-1)+2); % displaced nodal locations in y
        nodeCoordinates(nn,3) = nodeCoordinates(nn,3) + U(6*(nn-1)+3); % displaced nodal locations in z 

        % check for identical updated locations
        [Bb,~] = size(nodeCoordinates);
        [B,~] = size(unique(nodeCoordinates,'rows'));
        if Bb > B
            xk2 = find(setdiff(Bb,B))
            error('2')
        end 
    end

%% Add new bottom nodes and elements (at the surface)

    % reference to CURRENT bottom nodes of active cnts
    jj = nodeCount - numberBeams + activeCNTs;
    
    % displace CURRENT bottom nodes of active cnts to second row (growthNode locations)
    nodeCoordinates(jj,1) = growthNode(jj-nodeCount+numberBeams,1); % displaced nodal locations in x
    nodeCoordinates(jj,2) = growthNode(jj-nodeCount+numberBeams,2); % displaced nodal locations in y
    nodeCoordinates(jj,3) = growthNode(jj-nodeCount+numberBeams,3); % displaced nodal locations in z

    % check for identical updated locations
    [Bb,~] = size(nodeCoordinates);
    [B,~] = size(unique(nodeCoordinates,'rows'));
    if Bb > B
        xk3 = find(setdiff(Bb,B))
        error('3')
    end 

    % define nodes for each element (row = nodes of one element)
    elementNodes(element,1) = element;
    elementNodes(element,2) = nodeCount;
    
    % track current nodal count before new nodes added
    currentCount = nodeCount;
    
    % reference to NEW active bottom nodes
    kk = (currentCount + 1:currentCount + numberBeams)';

    % set locations of NEW bottom nodes at nucleation points
    nodeCoordinates(kk,1) = nucleationSite(kk-currentCount,1); % x-location of new nodes
    nodeCoordinates(kk,2) = nucleationSite(kk-currentCount,2); % y-location of new nodes
    nodeCoordinates(kk,3) = nucleationSite(kk-currentCount,3); % z-location of new nodes (= 0)

    % check for identical updated locations
    [Bb,~] = size(nodeCoordinates);
    [B,~] = size(unique(nodeCoordinates,'rows'));
%     if Bb > B
%         % check which position (row) has been duplicated: C = values, IA = indices
%         [C,IA] = setdiff(nodeCoordinates,unique(nodeCoordinates,'rows'),'rows','legacy')
%         error('4')
%     end 

    %%%%% overwrite locations of new inactive bottom nodes (extrapolate locations instead of using nucleation points)
        if inactiveCNTs
            % is error that causes NaN present
            if find(LL==00)
                error('5')
            end

            % reference to bottom nodes of inactive cnts
            gg = nodeCount + inactiveCNTs;
         
            % set locations of new inactive bottom nodes (vs. by displacement determined from force calculation)
            % (a) displace new nodes by same magnitude and direction of previous nodes
%             nodeCoordinates(gg,1) = nodeCoordinates(gg-2*numberBeams,1) - nodeCoordinates(gg-numberBeams,1); 
%             nodeCoordinates(gg,2) = nodeCoordinates(gg-2*numberBeams,2) - nodeCoordinates(gg-numberBeams,2); 
%             nodeCoordinates(gg,3) = nodeCoordinates(gg-2*numberBeams,3) - nodeCoordinates(gg-numberBeams,3);
            % (b) extrapolate from previous two nodes for displacement of new nodes
            nodeCoordinates(gg,1) = nodeCoordinates(gg-numberBeams,1) + (nodeCoordinates(gg-numberBeams,1) - nodeCoordinates(gg-2*numberBeams,1)); 
            nodeCoordinates(gg,2) = nodeCoordinates(gg-numberBeams,2) + (nodeCoordinates(gg-numberBeams,2) - nodeCoordinates(gg-2*numberBeams,2)); 
            nodeCoordinates(gg,3) = nodeCoordinates(gg-numberBeams,3) + (nodeCoordinates(gg-numberBeams,3) - nodeCoordinates(gg-2*numberBeams,3)); 
        
            % check for identical updated locations
            [Bb,~] = size(nodeCoordinates);
            [B,~] = size(unique(nodeCoordinates,'rows'));
            if Bb > B
                xk5 = find(setdiff(Bb,B))
                error('6')
            end          
        end
    
    % reference to new bottom row of nodes
    ii = nodeCount + 1:nodeCount + numberBeams;
    % reference to new bottom row of elements
    jj = element + 1:element + numberBeams;
    % update element number & nodal pairs
    elementNodes(jj,1) = jj';
    elementNodes(jj,2) = ii';
    
    % increase current node count by new row count of nodes
    nodeCount = nodeCount + numberBeams;
    % increase current element count by new row count of elements
    element = element + numberBeams;

    % periodic boundary condition check & adjustment
    if PeriodicBoundary > 0
         crossedxl = find(nodeCoordinates(:,1) < 0); % moves across left boundary
         nodeCoordinates(crossedxl,1) = nodeCoordinates(crossedxl,1) + h_span; 
         
         crossedxr = find(nodeCoordinates(:,1) > h_span); % moves across left boundary
         nodeCoordinates(crossedxr,1) = nodeCoordinates(crossedxr,1) - h_span;

         crossedy2 = find(nodeCoordinates(:,2) < 0); % moves across left boundary
         nodeCoordinates(crossedy2,2) = nodeCoordinates(crossedy2,2) + v_span;     
         
         crossedry2 = find(nodeCoordinates(:,2) > v_span); % moves across left boundary
         nodeCoordinates(crossedry2,2) = nodeCoordinates(crossedry2,2) - v_span;  
    end
end