function not_updated(O)


switch O.sys_order
    case 'first'
        switch O.type
            case 'p'
                update.G = false;
                update.CML = false;
                update.Ft = false;
                O.isupdated = update;
            case 'qp'
                update.H = false;
                update.kappa = false;
                update.theta = false;
                update.F_kappa = false;
                update.E    = false;
                O.isupdated = update;
        end
    case 'second'
        
        switch O.type
            case 'p'
                update.L = false;
                update.J = false;
                update.CML = false;
                update.DL = false;
                update.Ft = false;
                update.AA = false;
                update.Lint = false;
                O.isupdated = update;
            case 'qp'
                update.Q = false;
                update.kappa = false;
                update.theta = false;
                update.f_kappa = false;
                update.E    = false;
                O.isupdated = update;
        end
end