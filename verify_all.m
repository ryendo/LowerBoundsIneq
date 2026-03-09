
runner = VerificationRunner();  
% results = runner.verifyOmegaUp('J1');
% results = runner.verifyOmegaUp('J2');
results = runner.verifyOmegaMid('J1', [], 'inputs/cell_def_failed_J1.csv');
% results = runner.verifyOmegaMid('J2', [], 'inputs/cell_def.csv');
