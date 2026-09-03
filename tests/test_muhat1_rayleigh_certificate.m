function test_muhat1_rayleigh_certificate()
%TEST_MUHAT1_RAYLEIGH_CERTIFICATE Interval replay and tamper smoke test.

project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(project_root,'scripts_run'));
addpath(fullfile(project_root,'src','degenerate'));

if ~strcmp(getenv('RUN_INTLAB_SMOKE'),'1')
    error('test_muhat1_rayleigh_certificate:IntlabRequired', ...
        'Set RUN_INTLAB_SMOKE=1 and INTLAB_ROOT for this proof test.');
end

temporary_root = tempname;
mkdir(temporary_root);
cleanup = onCleanup(@() local_remove_tree(temporary_root));
output_dir = fullfile(temporary_root,'thin-smoke');

[manifest,validation] = run_muhat1_rayleigh_certificate( ...
    'mode','interval', ...
    't','0.38', ...
    'degree',10, ...
    's_cells',20, ...
    's_lo','0', ...
    's_hi','1', ...
    'down_y_max','0.06', ...
    'down_y_cells',10, ...
    'progress_every',0, ...
    'production',false, ...
    'output_dir',output_dir);

assert(~manifest.complete_certificate);
assert(validation.complete);
assert(~validation.manifest_complete_certificate);
assert(validation.number_of_s_cells == 20);
assert(validation.number_of_down_cells == 200);
assert(validation.all_masses_strictly_positive);
assert(validation.all_down_upper_bounds_strictly_negative);
assert(numel(validation.replayed_max_down_scaled_upper_hex) == 16);
assert(manifest.runtime_rechecked_after_run);
assert(ischar(manifest.runtime_sha256) ...
    && ~isempty(regexp(manifest.runtime_sha256, ...
        '^[0-9a-f]{64}$','once')));

manifest_path = fullfile(output_dir,'manifest.json');
did_reject_incomplete = false;
try
    validate_muhat1_rayleigh_certificate(manifest_path);
catch ME
    did_reject_incomplete = strcmp(ME.identifier, ...
        'validate_muhat1_rayleigh_certificate:InvalidCertificate');
end
assert(did_reject_incomplete, ...
    'The validator promoted a nonproduction smoke artifact to complete.');

coefficients_path = fullfile(output_dir,'frozen_coefficients.csv');
fid = fopen(coefficients_path,'at');
assert(fid >= 0);
file_cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'\n');
clear file_cleanup

did_reject = false;
try
    validate_muhat1_rayleigh_certificate( ...
        manifest_path,'require_complete_certificate',false);
catch ME
    did_reject = strcmp(ME.identifier, ...
        'validate_muhat1_rayleigh_certificate:InvalidCertificate');
end
assert(did_reject, ...
    'The validator accepted a frozen-trial file with the wrong hash.');

fprintf('test_muhat1_rayleigh_certificate: PASS\n');
clear cleanup
local_remove_tree(temporary_root);
end


function local_remove_tree(path_value)
if isfolder(path_value)
    rmdir(path_value,'s');
end
end
