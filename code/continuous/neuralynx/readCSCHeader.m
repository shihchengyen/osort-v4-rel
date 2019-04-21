function header = readCSCHeader(fpath)

fid = fopen(fpath, 'r');
if fid == -1
    error ('IMPORT:fileNotFound','failed to open the file: file not found or permission denied');
end

% read file header (16 kilobyte ASCII)
header = transpose(fread(fid, 16*1024, '*char'));

fclose(fid);

