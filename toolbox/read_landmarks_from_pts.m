function [px,py] = read_landmarks_from_pts(f)


fid = fopen(f);
tline = fgetl(fid);
idx = 0;
px = zeros(68,1);
py = zeros(68,1);
relidx = 0;
while ischar(tline)
    idx = idx+1;
    tline = fgetl(fid);
    if idx <=2 || idx >= 71
        continue
    end
    relidx = relidx+1;
    
    tmp = strsplit(tline, ' ');
    x = str2double(tmp{1});
    y = str2double(tmp{2});
    px(relidx) = x; 
    py(relidx) = y;
    
end

fclose(fid);

