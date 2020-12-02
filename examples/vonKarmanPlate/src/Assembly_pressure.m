function [pload]=Assembly_pressure(model,p)

nelem = model(1).nelem;
ndof = model(1).ndof;

pload=sparse(ndof,1);

for i=1:nelem    
    
    el_load = zeros(18,1); el_load([3 9 15]) = p*model(i).area/3;
    pload_el=model(i).T'*el_load;
    pload(model(i).index,1)=pload(model(i).index,1)+pload_el;
    
end
