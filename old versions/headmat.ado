prog def headmat
syntax name, [Rows(int 5)]
confirm matrix `namelist'
tempname headmat
mat `headmat' = `namelist'[1..`rows',.]
di as text "First `rows' rows of `namelist':" _c
mat l `headmat', noheader
end
