fid = fopen('1GTO.rc');
%fid = fopen('withheader-1GTO.grid');
isovalue = 0.1;
gridsize= fscanf(fid, '%d',1);
rawdata =  fscanf(fid, '%f');
fclose(fid);

rawindex = 1;
volumedata = zeros(gridsize, gridsize, gridsize);
for x=1:gridsize
    for y=1:gridsize
        for z=1:gridsize
            volumedata(x,y,z) = rawdata(rawindex);
            rawindex = rawindex + 1;
        end
    end
end

isosurface(volumedata, isovalue);