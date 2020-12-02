function CM = get_ConvMtx(O)
if ~O.isupdated.CML
    update_ConvMtx(O);
end
CM = O.ConvMtx;
end