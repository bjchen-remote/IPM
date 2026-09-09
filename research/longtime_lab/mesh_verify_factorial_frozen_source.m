function mesh_verify_factorial_frozen_source(sourceDirectory)
%MESH_VERIFY_FACTORIAL_FROZEN_SOURCE Verify registered bytes and path closure.
manifest=jsondecode(fileread(fullfile(sourceDirectory,'sha256_manifest.json')));
for k=1:numel(manifest)
    file=fullfile(sourceDirectory,manifest(k).path);fid=fopen(file,'rb');assert(fid>=0);
    bytes=fread(fid,Inf,'*uint8');fclose(fid);
    digest=java.security.MessageDigest.getInstance('SHA-256');digest.update(bytes);
    actual=reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]);
    assert(strcmpi(actual,manifest(k).sha256),'ipm:BoxSpaceFrozenSource','Frozen source SHA mismatch: %s',file);
end
files=matlab.codetools.requiredFilesAndProducts(which('ipm.solve'));
assert(all(startsWith(files,[sourceDirectory,filesep])), ...
    'ipm:BoxSpaceFrozenSource','Solver dependency escaped the registered source.');
end
