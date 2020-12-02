function AA = get_AA(O)
if ~O.isupdated.AA
    update_AA(O);
end
AA = O.AA;
end