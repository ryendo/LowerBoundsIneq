function test_omega_mid_geometric_coverage()
%TEST_OMEGA_MID_GEOMETRIC_COVERAGE Exact-decimal, fail-closed cover test.

if ~strcmp(getenv('RUN_INTLAB_SMOKE'),'1')
    fprintf(['Omega_mid coverage test skipped; set RUN_INTLAB_SMOKE=1 ', ...
        'and INTLAB_ROOT to enable it.\n']);
    return;
end

project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(project_root,'scripts_run'));
addpath(fullfile(project_root,'tools'));
global INTERVAL_MODE
if isempty(INTERVAL_MODE) || ~logical(INTERVAL_MODE)
    omega_up_all_prepare_worker(project_root,'interval');
end

temporary_dir = tempname;
mkdir(temporary_dir);
cleanup = onCleanup(@() local_remove_tree(temporary_dir));

% Numerically equal endpoint spellings must join exactly.  This two-cell
% rectangle is intentionally much larger than the target and therefore is
% a cheap positive fixture.
touching_csv = fullfile(temporary_dir,'touching.csv');
local_write_fixture(touching_csv,{ ...
    '0.500000000000000000','0.700000000000000000','0.03','1.0'; ...
    '0.7','1.000000000000000000','0.0300','1.000'});
touching_json = fullfile(temporary_dir,'touching_coverage.json');
touching = validate_omega_mid_coverage( ...
    touching_csv,'manifest_file',touching_json);
assert(touching.complete);
assert(touching.rigorous);
assert(strcmp(touching.status,'certified'));
assert(touching.proof.decimal_topology_exact);
assert(~touching.proof.binary_double_endpoint_roundtrip);
assert(~touching.proof.arbitrary_tolerance_used);
assert(numel(touching.hashes.proof_payload) == 64);
assert(exist(touching_json,'file') == 2);
decoded = jsondecode(fileread(touching_json));
assert(decoded.complete && decoded.rigorous);
assert(strcmp(decoded.hashes.proof_payload, ...
    touching.hashes.proof_payload));

% This gap is 10^-18 wide.  Its two endpoints normally collapse to the
% same binary double, but the exact decimal arrangement must retain the
% open x slab and reject the cover without a tolerance.
gap_csv = fullfile(temporary_dir,'one_attometre_gap.csv');
local_write_fixture(gap_csv,{ ...
    '0.5','0.700000000000000000','0.03','1'; ...
    '0.700000000000000001','1','0.03','1'});
gap = validate_omega_mid_coverage(gap_csv);
assert(~gap.complete);
assert(strcmp(gap.status,'failed'));
assert(~isempty(gap.witness));
assert(strcmp(gap.witness.kind,'uncovered_elementary_x_slab'));
assert(strcmp(gap.witness.x_inf_decimal,'0.7'));
assert(strcmp(gap.witness.x_sup_decimal,'0.700000000000000001'));

% Finally check the production geometry itself.  This performs only CSV
% parsing, a roughly 1.2-million-entry logical arrangement, and scalar
% boundary interval evaluations; it dispatches no FEM calculation.
canonical_csv = fullfile(project_root,'inputs','cell_def.csv');
started = tic;
canonical = validate_omega_mid_coverage(canonical_csv);
elapsed = toc(started);
assert(canonical.complete);
assert(canonical.input.cell_count == 188623);
assert(strcmp(canonical.input.sha256, ...
    ver10_file_sha256(canonical_csv)));
assert(strcmp(canonical.hashes.input_csv,canonical.input.sha256));
assert(canonical.statistics.unique_x_coordinates == 1637);
assert(canonical.statistics.unique_theta_coordinates == 723);
assert(canonical.statistics.target_x_slabs_checked > 0);
assert(canonical.statistics.minimum_certified_lower_margin >= 0);
assert(canonical.statistics.minimum_certified_upper_margin >= 0);

fprintf(['test_omega_mid_geometric_coverage: PASS ', ...
    '(canonical cover %.3f s; 10^-18 gap rejected)\n'],elapsed);
fprintf('canonical proof payload SHA-256: %s\n', ...
    canonical.hashes.proof_payload);
fprintf('certified theta margins: lower=%.17g, upper=%.17g\n', ...
    canonical.statistics.minimum_certified_lower_margin, ...
    canonical.statistics.minimum_certified_upper_margin);
clear cleanup
local_remove_tree(temporary_dir);
end


function local_write_fixture(filename,rows)
fid = fopen(filename,'w');
if fid < 0
    error('test_omega_mid_geometric_coverage:FixtureOpenFailed', ...
        'Could not create %s.',filename);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'x_inf,x_sup,theta_inf,theta_sup\n');
for k = 1:size(rows,1)
    fprintf(fid,'%s,%s,%s,%s\n', ...
        rows{k,1},rows{k,2},rows{k,3},rows{k,4});
end
clear cleanup
end


function local_remove_tree(pathname)
if exist(pathname,'dir')
    rmdir(pathname,'s');
end
end
