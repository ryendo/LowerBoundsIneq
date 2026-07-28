function digest = ver10_file_sha256(filename)
%VER10_FILE_SHA256 Compute a streaming SHA-256 digest without shell paths.

fid = fopen(filename,'rb');
if fid < 0
    error('ver10_file_sha256:CannotOpen', ...
        'Cannot open %s for hashing.',filename);
end
cleanup = onCleanup(@() fclose(fid));

try
    engine = javaMethod( ...
        'getInstance','java.security.MessageDigest','SHA-256');
catch ME
    error('ver10_file_sha256:NoSHA256', ...
        'Cannot initialize the SHA-256 engine: %s',ME.message);
end

while true
    bytes = fread(fid,1024*1024,'*uint8');
    if isempty(bytes)
        break;
    end
    engine.update(typecast(bytes(:),'int8'));
end
raw = typecast(engine.digest(),'uint8');
digest = lower(reshape(dec2hex(raw,2).',1,[]));

clear cleanup
end
