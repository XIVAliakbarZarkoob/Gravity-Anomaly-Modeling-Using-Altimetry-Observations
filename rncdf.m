%% Aliakbar Zarkoob, AKA "XIV"
%  Gmail: XIV.Aliakbar.Zarkoob@gmail.com
%  Telegram: @XIVAliakbar

function data = rncdf(file)

    ncid = netcdf.open(file, 'NC_NOWRITE');
    [~, numvars, ~, ~] = netcdf.inq(ncid);
    
    varname = string();
    for i = 0:numvars-1
        [varname(i+1), ~, ~, ~] = netcdf.inqVar(ncid, i);     
        data.(strrep(varname(i+1),'.','_')) = netcdf.getVar(ncid, i);
    end
    netcdf.close(ncid);
    data.varname = varname;
    
end