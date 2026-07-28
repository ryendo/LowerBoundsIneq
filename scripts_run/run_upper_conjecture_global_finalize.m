%RUN_UPPER_CONJECTURE_GLOBAL_FINALIZE  Join compact and residual manifests.

compact_manifest = getenv('VER10_COMPACT_UPPER_MANIFEST');
residual_manifest = getenv('VER10_OMEGA_UP_MANIFEST');
certificate_root = getenv('VER10_CERTIFICATE_ROOT');
if isempty(compact_manifest) || isempty(residual_manifest) ...
        || isempty(certificate_root)
    error('run_upper_conjecture_global_finalize:MissingEnvironment', ...
        ['Set VER10_COMPACT_UPPER_MANIFEST, VER10_OMEGA_UP_MANIFEST, ', ...
         'and VER10_CERTIFICATE_ROOT.']);
end
output_path = fullfile( ...
    certificate_root,'siudeja_upper_global_manifest.json');
manifest = verify_upper_conjecture_global_manifests( ...
    compact_manifest,residual_manifest,output_path);
if ~manifest.complete_certificate
    error('run_upper_conjecture_global_finalize:Incomplete', ...
        'The combined global upper certificate is incomplete.');
end
fprintf('global Siudeja upper certificate: %s\n',output_path);
