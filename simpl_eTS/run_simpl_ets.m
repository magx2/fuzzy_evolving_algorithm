function run_simpl_ets

    % y = 2x_1 + 3x_2 + 4
    x = [
        1 2 3 4;
        1 2 4 8 
    ]

    y = [
        9 14 22 36
    ]

    r = 1
    
    OMEGA = 1;
    
    simpl_ets( x, y, r, OMEGA )
end