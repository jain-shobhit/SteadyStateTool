function startup()
    % initial notes
    disp('SteadyStateTool: Fast computation of Steady-State response of multi-degree of freedom nonlinear mechanical systems')
    disp('Maintained by Shobhit Jain: shjain@ethz.ch')
    disp('The examples from the following reference are available in the examples folder: ')
    disp(['S. Jain, T. Breunung, G. Haller, Fast Computation of Steady-State Response for ' ...
        'Nonlinear Vibrations of High-Degree-of-Freedom Systems, Nonlinear Dyn (2019) 97: ' ...
        '313. https://doi.org/10.1007/s11071-019-04971-1'])
    disp('Please cite this reference if you use this package in your work.')
    disp('You are now ready to play with the existing examples and create your own.')
    
    % plotting preferences
    set(groot,'defaultAxesTickLabelInterpreter', 'latex')
    set(groot,'defaultColorbarTickLabelInterpreter', 'latex')
    set(groot,'defaultLegendInterpreter', 'latex')
    set(groot,'defaultTextInterpreter', 'latex')
    set(groot,'defaultAxesFontSize', 14)
    set(groot,'defaultLegendFontSize', 14)
    
    % adding current folder and subfolders to matlab path
    addpath(genpath(pwd))
end