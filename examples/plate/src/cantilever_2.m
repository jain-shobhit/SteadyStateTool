function [nodes,elements,load,freedofs,prop]=cantilever_2(L,H,n,m,ll,w,t)

% L = lenght
% H = height
% n = number of divisions in L direction
% m = number of divisions in H direction

% the plate lies in the yz plane
% x direction is the out-of-plane direction

seeds=1:(n+1)*(m+1);
seeds=reshape(seeds,n+1,m+1)';

count=1;
for j=1:m+1
    for i=1:n+1
        nodes(count,:)=[count ((i-1)/n)*L H-(j-1)/m*H 0];
        count=count+1;
    end
end

% rotating the plate to xz plane
% x=nodes(:,2); y=nodes(:,3); z=nodes(:,4);


% add curvature

if w~=0
R = 1/2*(w^2+(H/2)^2)/w;
theta0 = asin(H/2/R);

for i = 1:size(nodes,1)
    th = asin((nodes(i,3)-H/2)/R);
    nodes(i,4) = R*cos(th)-cos(theta0);
end
end


count=1;
for j=1:m
    for i=1:n
    
      elements(count:count+1,:)=[count   seeds(j,i)   seeds(j+1,i+1)   seeds(j,i+1);
                                 count+1 seeds(j,i) seeds(j+1,i)     seeds(j+1,i+1)];
  
      count=count+2;
       
   end
end

elements = [elements'; ones(1,size(elements,1))]';
% ----- Boundary conditions...


% for axial compression, restrain only x 

if  strcmp (ll,'C') || strcmp (ll ,'S') ;
    % side W
    bnodes1=1:n+1:m*(n+1)+1;
	bounddof1=bnodes1*6-3; % bnodes1*6-1]; % bnodes1*6-3];
	% side E
	bnodes2=n+1:n+1:m*(n+1)+n+1;
	bounddof2=bnodes2*6-3; % bnodes2*6-1]; % bnodes2*6-3];
	% side N
	bnodes3=1:n+1;
	bounddof3= bnodes3*6-3; % bnodes3*6]; % bnodes3*6]; % bnodes3*6-3];
	% side S
	bnodes4=m*(n+1)+1:m*(n+1)+n+1;
	bounddof4= bnodes4*6-3; % bnodes4*6]; % bnodes4*6]; % bnodes4*6-3];

% top corner

elseif strcmp(ll, 'B')
    % side W
    bnodes1=1:n+1:m*(n+1)+1;
% 	bounddof1=[bnodes1*6-5  bnodes1*6-4 bnodes1*6-3];
    bounddof1=[];
	% side E
	bnodes2=n+1:n+1:m*(n+1)+n+1;
% 	bounddof2=[bnodes2*6-5  bnodes2*6-4 bnodes2*6-3];
    bounddof2=[];
	% side N
	bnodes3=1:n+1;
	bounddof3=[ bnodes3*6-5 bnodes3*6-4 bnodes3*6-3];
	% side S
	bnodes4=m*(n+1)+1:m*(n+1)+n+1;
	bounddof4=[ bnodes4*6-5 bnodes4*6-4 bnodes4*6-3]; 
%     bounddof4= [];
    elseif strcmp(ll , 'DC')
    % side W
    bnodes1=1:n+1:m*(n+1)+1;
% 	bounddof1 = [bnodes1*6-5 bnodes1*6-4 bnodes1*6-3 bnodes1*6-2 bnodes1*6-1 bnodes1*6]; % bnodes1*6-3];
	bounddof1 = [];
    % side E
	bnodes2=n+1:n+1:m*(n+1)+n+1;
	bounddof2 = []; 
%     bounddof2 = [bnodes2*6-5 bnodes2*6-4 bnodes2*6-3 bnodes2*6-2 bnodes2*6-1 bnodes2*6]; % bnodes2*6-3];
% 	% side N
 	bnodes3=1:n+1;
 	% bounddof3 = [bnodes3*6-5 bnodes3*6-4 bnodes3*6-3]; % bnodes3*6-3];
%     bounddof3 = [];
    bounddof3 = [bnodes3*6-5 bnodes3*6-4 bnodes3*6-3 bnodes3*6-2 bnodes3*6-1 bnodes3*6]; % bnodes1*6-3];

% 	% side S
 	bnodes4=m*(n+1)+1:m*(n+1)+n+1;
 	bounddof4=[ bnodes4*6-5 bnodes4*6-4 bnodes4*6-3 bnodes4*6-2 bnodes4*6-1 bnodes4*6 ]; %[ bnodes4*6-5 bnodes4*6-4 bnodes4*6-3];
%     bounddof4 = [];
    

end

% extrabounddof = 0; %(n+1)*6-4;

%bounddof=unique([2 3 extrabounddof,bounddof1, bounddof2, bounddof3, bounddof4]);

%  for rigib body rotation
bounddof=unique([bounddof1, bounddof2 , bounddof3, bounddof4]);
ndof=6*size(nodes,1);
alldof=1:ndof;

freedofs=setdiff(alldof,bounddof);

% ----- External load...

load=sparse(ndof,1);

% ----- compression...

if strcmp(ll,'C')
	loadeddof1=bnodes1*6-5;
	loadeddof2=bnodes2*6-5;
	
	load(loadeddof1)=1/m;
	load(loadeddof1(1))=0.5/m; load(loadeddof1(end))=0.5/m;
	load(loadeddof2)=-1/m;
	load(loadeddof2(1))=-0.5/m; load(loadeddof2(end))=-0.5/m;

	cc = m*(n+1)+1;
	% 
	m0=L/12/m/m;
	% 
	load(4) = -m0;
	load(6*cc-2) = m0;
	load((n+1)*6-2)=m0;
	load(6*cc+6*n-2)=-m0;
end

% ----- central bending load...

if strcmp(ll ,'B')|| strcmp(ll, 'DC')
    %loadednode=(n+1)*m/2+1+n/2; % center node
    % loadednode = (n+1)*(m/2+1);  % center node at the tip
    loadednode = n+1; % tip corner
    load(6*loadednode-3)=-1;
end

% ----- shear load...

if strcmp(ll ,'S')
    loadeddof1 = bnodes1*6-4;
    loadeddof2 = bnodes2*6-4;
    loadeddof3 = bnodes3*6-5;
    loadeddof4 = bnodes4*6-5;
    load(loadeddof1(2:end-1)) = 1/m;
    load(loadeddof1([1 end])) = 0.5/m;
    load(loadeddof2(2:end-1)) = -1/m;
    load(loadeddof2([1 end])) = -0.5/m;
    load(loadeddof3(2:end-1)) = -1/n*L/H;
    load(loadeddof3([1 end])) = -0.5/n*L/H;
    load(loadeddof4(2:end-1)) = 1/n*L/H;
    load(loadeddof4([1 end])) = 0.5/n*L/H;
end
%% Units in mm - Force is in Newton
E = 70e3;       % MPa
nu = 0.33;      
rho = 2700e-12; % Megagram/mm^3
prop=[1 E nu t rho];
%% SI units 
% E = 70e9;       % Pa
% nu = 0.33;      
% t = 8e-4;        % m
% rho = 2700; % kg/m^3
% prop=[1 E nu t rho];
