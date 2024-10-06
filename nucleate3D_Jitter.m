function [elementNodes,ang,rate,nodeCoordinates,nodeCount,element,nucleationSite,growthNode,phi]=nucleate3D_Jitter(numberBeamsx,numberBeamsy,h_span,v_span,avgRate,rate_stdev,ang_stdev)

for o=1:numberBeamsx;
    origin1(o)=h_span/(numberBeamsx)*(o-1/2);
end

for p=1:numberBeamsy
    origin2(p)=v_span/(numberBeamsy)*(p-1/2);
end
    [xx,yy]=meshgrid(origin1,origin2);
    size(xx);
    size(yy);
    numberBeams=numberBeamsx*numberBeamsy;

    jitterx=-35e-9+70e-9*rand(numberBeamsx,numberBeamsy)';
    jittery=-35e-9+70e-9*rand(numberBeamsx,numberBeamsy)';
    
    jitterx=-15e-9+7.5e-9*rand(numberBeamsx,numberBeamsy)';
    jittery=-15e-9+7.5e-9*rand(numberBeamsx,numberBeamsy)';
    
    size(jitterx);
    size(jittery);
   [x]=xx+jitterx;
   [y]=yy+jittery;

   %% NORMAL DISTRIBUTION
   
% % % for ii=1:numberBeams;
% % %     mu_ang = 0;
% % %     Sigma_ang = (ang_stdev*3.1415/180)^2; %Converted to radians
% % %     R_ang = chol(Sigma_ang);
% % %     ang(ii) = repmat(mu_ang,1) + randn(1)*R_ang;
% % %     ro(ii)=rout;
% % %     rinn(ii)=ri;
% % %     phi(ii)=rand*pi;
% % %     mu_rate=avgRate;
% % %     Sigma_rate = [(rate_stdev/100*avgRate)^2;];%changed
% % %     R_rate = chol(Sigma_rate);
% % %     rate(ii) = repmat(mu_rate,1) + randn(1)*R_rate;  
% % % end


for ii=1:numberBeams;
    mu_ang = 0;
    Sigma_ang = (ang_stdev*3.1415/180)^2; %Converted to radians
    R_ang = chol(Sigma_ang);
    ang(ii) = repmat(mu_ang,1) + randn(1)*R_ang;
    phi(ii)=rand*pi;
    mu_rate=2;
    Sigma_rate = rate_stdev;%changed
    rate(ii) = logninv(rand,mu_rate,Sigma_rate)*10e-9;  
end

nodeCount=0;
element=0;
for num = numberBeams+1:2*numberBeams %% Setting up CNT bases
    nodeCount=nodeCount+1;
    nodeCoordinates(num,1)=x(num-numberBeams);
    nodeCoordinates(num,2)=y(num-numberBeams); 
    nodeCoordinates(num,3)=0;
    nucleationSite(num-numberBeams,1)=nodeCoordinates(num,1);%%xxx(num)
    nucleationSite(num-numberBeams,2)=nodeCoordinates(num,2);%%yyy(num)
    nucleationSite(num-numberBeams,3)=nodeCoordinates(num,3);%%zzz(num)
end

%%Making the free tips of the CNTs the initial nodes
for num=1:numberBeams;%%Setting position of CNT free ends
        element=element+1; 
        nodeCount=nodeCount+1;
        nodeCoordinates(num,1)=nodeCoordinates(num+numberBeams,1)+sin(ang(num))*cos(phi(num))*rate(num);%changed '*' to '+' before sin
        nodeCoordinates(num,2)=nodeCoordinates(num+numberBeams,2)+sin(ang(num))*sin(phi(num))*rate(num);%changed '*' to '+' before sin
        nodeCoordinates(num,3)=nodeCoordinates(num+numberBeams,3)+cos(ang(num))*rate(num);%changed '*' to '+' before cos
        growthNode(num,1)=nodeCoordinates(num,1);
        growthNode(num,2)=nodeCoordinates(num,2);
        growthNode(num,3)=nodeCoordinates(num,3);
        elementNodes(element,1)=element;
        elementNodes(element,2)=nodeCount;
end   



