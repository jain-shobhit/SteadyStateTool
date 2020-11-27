function [model,Fext,constr_nodes,bounddofs,ndof,freedofs,nodes,elements] = generate_model(mshfile,p,BC,units)
%Generate Model. 
%Generate model based on file "mshfile" generated using GMSH and pressure
%in z direction over nbc nodes. 
disp('Generating the model...')

% cd GMSH %% change Directory

% mshfile= 'Mesh1.msh'; % define variable for file name

[dmn_ids,ebc_ids,nbc_ids]=n_pid(mshfile);

%function reads the mesh file and returns 
%the different physical domain ids -dmn_ids , (these are numbered from 1-100)- domain ID corresponds to the property ID of these elements 
%the different natural BC ids - ebc_ids,(these are numbered from 101-200)
%the different essential BC ids - nbc_ids (these are numbered from 201-300)


%% generating nodes Matrix

%nodes contains matrix: #node_id #coor_x #coor_y #coor_z
[nodeid,node]=readnodes(mshfile);
nodes = [nodeid,node];
all_nodes = []; 
for i=1:size(dmn_ids,1)
all_nodes = [all_nodes; readnodeset(mshfile,dmn_ids(i))]; % read node ids of triangular shell elements
end
all_nodes = unique(all_nodes);
red_nodes = setdiff(all_nodes,nodeid');

%% generating element connectivity Matrix
%element contains matrix: #element_id #node1 #node2 #node3 #property number
% [conn,ele_id,ele_typ]=readelements(mshfile,prop_id(1,1));
% elements = [ele_id,conn,prop_id(1,2)*ones(size(ele_id,1))];
elements=[];
for i=1:size(dmn_ids,1)

[conn,ele_id,~]=readelements(mshfile,dmn_ids(i)); %return elements belonging to particular physical ID. ele_id, eletyp are row matrices
elements = [elements;[ele_id',conn,dmn_ids(i)*ones(size(ele_id,2),1)]];

end

nelm = length(elements(:,1));
el_id = elements(:,1);  % el_id for mapping
elements(:,1) = 1:nelm; % renumber elements
%% shell properties of the different wing sections
% prop = [prod_id E nu thickness density] all in mm
%  prop = [1 6.9e7 0.33 2 2.7e-6;
%          2 6.9e7 0.33 2 2.7e-6;
%          3 6.9e7 0.33 2 2.7e-6;
%          4 6.9e7 0.33 0.5 2.7e-6];
if strcmp(units, 'SI')
E = 70e9; nu = 0.33; t = 1.5e-3; rho = 2700; %% SI Units
elseif strcmp(units, 'mm')
E = 70e3; nu = 0.33; t = 1.5; rho = 2700e-12; %% mm units
end
prop=[1 E nu t rho;   % thickness for skin elements
      2 E nu t rho;  % thickness for stiffener elements
      3 E nu t rho];  % thickness for stiffener elements




%% creates a structure with all the useful quantities for the triangular 3 node shell element     
model = model_database_2(nodes,elements,prop);
%% nodes that are constrained
constr_nodes = readnodeset(mshfile,ebc_ids);
%% corresponding tied dofs
if strcmp(BC, 'SS')
    bounddofs = [constr_nodes*6-5; constr_nodes*6-4; constr_nodes*6-3; red_nodes*6-5; red_nodes*6-4; red_nodes*6-3; red_nodes*6-2; red_nodes*6-1; red_nodes*6]; % contrain translation in 3 directions
elseif strcmp(BC, 'C')
    bounddofs = [constr_nodes*6-5; constr_nodes*6-4; constr_nodes*6-3; constr_nodes*6-2;constr_nodes*6-1;constr_nodes*6; red_nodes*6-5; red_nodes*6-4; red_nodes*6-3; red_nodes*6-2; red_nodes*6-1; red_nodes*6]; % contrain translation in 3 directions
end
ndof=size(nodes,1)*6;

%% free dofs
freedofs = setdiff(1:ndof,bounddofs);
%% generating Fext Matrix
%force contains matrix: #node_id #F(vector representing force perunit length)
%NOTE, no direction is given
% F = [0 0 -100]; % Force per unit length of edge
% F_nodes = readnodeset(mshfile,nbc_ids(1)); % Nodes containing Force / Pressure
[~,nbc_eid,~]=readelements(mshfile,nbc_ids(1)); % NBC elements

%create zero vector which is 6 times length of node number
Fext=sparse(ndof,1); % size(nodes,1) returns no. of rows in matrix 'nodes' i.e. no. of nodes, At every node, six types of force inputs can be given

% force = Assembly_Force(F_nodes,nbc,nodes,F);
% %this sets the force amplitude at every 3th(zdirection) of loads vector
% Fext(force(:,1)*6-5) = force(:,2); % 1st column of Force matrix has node nos. and 2nd column has force (in Z direction) value in that node. 
% Fext(force(:,1)*6-4) = force(:,3); %In our representation, 1st 6 rows
% Fext(force(:,1)*6-3) = force(:,4); %correspond to forces on 1st node. So
%                                     %force in z direction on node 1 is in
%                                     %row 3. Similarly force on 2nd node
%                                     %will be in row 9, nth node will be on
%                                     %6n-3 th row.
                                    
% pressure on NBC area
nbc_el = [];
for n= nbc_eid   
    i = find(el_id==n-1);% duplicate element for force boundary with index = original index + 1
    el_load = zeros(18,1); el_load([3 9 15]) = p*model(i).area/3; % pressure load in y direction
    pload_el=model(i).T'*el_load;
    Fext(model(i).index,1)=Fext(model(i).index,1)+pload_el;
    nbc_el = [nbc_el; i];    
end
model(1).nbc_el = nbc_el;
% def_plot(zeros(length(F_nodes)*6,1),1,nodes(F_nodes,:),elements(nbc_el,:));
    trisurf(elements(:,2:4),nodes(:,2),nodes(:,3),nodes(:,4),0*ones(size(nodes,1),1))    
    axis equal
    axis off
    hold on
    trisurf(elements(nbc_el,2:4),nodes(:,2),nodes(:,3),nodes(:,4),10*ones(size(nodes,1),1))

% el = [];
% for i=1:2
% 
% [con,el_id,~]=readelements(mshfile,dmn_ids(i)); %return elements belonging to particular physical ID. ele_id, eletyp are row matrices
% el = [el;[el_id',con,dmn_ids(i)*ones(size(el_id,2),1)]];
% 
% end
% 
% figure
% trisurf(el(:,2:4),nodes(:,2),nodes(:,3),nodes(:,4),0*ones(size(nodes,1),1))    
% axis equal

%     for i = 1:size(nodes,1)
%      text(nodes(i,2),nodes(i,3),nodes(i,4), num2str(i));
%     end
%     for i = 1:nelm
%         x1 = nodes(elements(i,2),2);
%         x2 = nodes(elements(i,3),2);
%         x3 = nodes(elements(i,4),2);
%          y1 = nodes(elements(i,2),3);
%         y2 = nodes(elements(i,3),3);
%         y3 = nodes(elements(i,4),3);
%          z1 = nodes(elements(i,2),4);
%         z2 = nodes(elements(i,3),4);
%         z3 = nodes(elements(i,4),4);
%      text((x1+x2+x3)/3,(y1+y2+y3)/3,(z1+z2+z3)/3, num2str(i));
%     end
    hold off

end
