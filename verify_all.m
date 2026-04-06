
runner = VerificationRunner();  
results_J1Up = runner.verifyOmegaUp('J1');
results_J2Up = runner.verifyOmegaUp('J2');
results_J1Mid = runner.verifyOmegaMid('J1', [], 'inputs/cell_def.csv');
results_J2Mid = runner.verifyOmegaMid('J2', [], 'inputs/cell_def.csv');
