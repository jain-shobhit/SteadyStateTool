function update_Ft(O)

switch O.sys_order
    case 'first'
        O.Ft = O.F(O.t,O.T);
    case 'second'
        O.Ft = O.f(O.t,O.T);
end
O.isupdated.Ft = true;
end