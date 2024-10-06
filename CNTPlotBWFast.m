function CNTPlotBWFast(fname,avgRate, steps, h_span,v_span, numberBeams,ContinuousPlot,nodeCount,nodeCoordinates,t,totalCompress, compressiveLoad,numberBeamsx,numberBeamsy,inactiveNodes)

hold off
allnodes = 1:nodeCount;
activeNodes = setdiff(allnodes,inactiveNodes');

% Adjust grayscale to offest whitest pixels (front = white, back = black)
scale = max(nodeCoordinates(activeNodes,1)) - min(nodeCoordinates(activeNodes,1)); % distance between furthest front and backmost nodes
normalizedX = (min(nodeCoordinates(activeNodes,1)) + nodeCoordinates(activeNodes,1)) ./ scale; % normalized distance of nodes in x from 0 to 1

%array = (nodeCoordinates(activeNodes,1)/ (h_span)).^3;
array = normalizedX ./ (scale).^2; % array of normalized x values for activeNodes and divide by factor of scale: factor determines how quickly grayscale changes in x
color = [0 0 0] + array./1.3; % division by 1.1 offsets grayscale so that whitest pixels are not pure white/blend with background

scatter3(nodeCoordinates(activeNodes,1)*1e6,nodeCoordinates(activeNodes,2)*1e6, nodeCoordinates(activeNodes,3)*1e6,0.25,color, 'filled')
hold on
%scatter3(nodeCoordinates(activeNodes,1)*1e6,(nodeCoordinates(activeNodes,2)+v_span)*1e6, nodeCoordinates(activeNodes,3)*1e6,0.5,color, 'filled')
% *** pick scale
axis([-1*h_span*1e6 h_span*1e6 -3-v_span*1e6 3+v_span*1e6 0 50])            
axis off
            daspect([1 1 1])
            figureHandle = gcf;
            set(gca,'FontSize',22)
            hold on
            %view(96,15)
            view(96,5)
           
% Add artificial points between nodeCoordinates to make plots more full
l = find(activeNodes<nodeCount-numberBeams);
p = activeNodes(l);
for i = 0.2:0.2:0.8
    fakeCoordinates(p,:) = nodeCoordinates(p,:)+i*(nodeCoordinates(p,:) - nodeCoordinates((p+numberBeams),:));
    arrayf = (min(nodeCoordinates(activeNodes,1))+(fakeCoordinates(p,1))/(h_span)).^3;
    colorf = [0 0 0] + arrayf/1.1;
    scatter3(fakeCoordinates(p,1)*1e6,fakeCoordinates(p,2)*1e6, fakeCoordinates(p,3)*1e6,0.25, colorf, 'filled')
    hold on
end
hold off

%scatter3(fakeCoordinates(p,1)*1e6,(fakeCoordinates(p,2)+v_span)*1e6, fakeCoordinates(p,3)*1e6,0.5,colorf, 'filled')
% ll = find(activeNodes<nodeCount-numberBeams);
% pp=activeNodes(ll);
% fakeCoordinates(pp,:)=(nodeCoordinates(pp,:)+ nodeCoordinates((pp+numberBeams),:))*1/3;
% arrayf = fakeCoordinates(pp,1)/(1.05*h_span);
% colorff = [0 0 0] + arrayf;
% scatter3(fakeCoordinates(pp,1)*1e6,fakeCoordinates(pp,2)*1e6, fakeCoordinates(pp,3)*1e6,0.5,colorff, 'filled')
% hold on
       
plotname = strcat(fname,num2str(t));  
print(gcf,plotname,'-dpng','-r300');

       %pause(5)
       %saveas(gcf,fullfile(fname,plotname),'jpg');
       %saveas(gcf,fullfile(fname,plotname),'fig');
            close()


