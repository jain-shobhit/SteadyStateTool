function Lint = get_Lint(O)
if ~O.isupdated.Lint
    update_Lint(O);
end
Lint = O.Lint;
end