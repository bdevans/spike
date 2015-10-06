function genReport(rdir)
% Function to generate an html report for a batch of simulations

extPath = '/Network/Servers/mac0.cns.ox.ac.uk/Volumes/Data/Users/evans/Sites/';

if ~strcmp(rdir(end),filesep); rdir = strcat(rdir,filesep); end
cd(rdir);
%system(['cd ',rdir]);
fid = fopen([rdir,'report.html'],'w'); % Change to last dir of rdir

ind = find(rdir == filesep);
sdir = rdir(ind(end-1)+1:end-1);

% Create header
fprintf(fid,'<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">\n');
fprintf(fid,'<html xmlns="http://www.w3.org/1999/xhtml">\n');
fprintf(fid,'<head>\n');
fprintf(fid,'<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />\n');
fprintf(fid,'<title>%s</title>\n',sdir);
fprintf(fid,'<script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jquery/1.6/jquery.min.js"></script>\n');
fprintf(fid,'<script type="text/javascript" src="%stoggleDisplay.js"></script>\n',extPath);
fprintf(fid,'<link rel="stylesheet" type="text/css" href="%stheme.css" />\n',extPath);
fprintf(fid,'</head>\n');
% Favicons (16X16)
%<link rel="icon" type="image/x-ico" href="http://www.yourwebsite.com/favicon.ico" />
%<link rel="shortcut icon" type="image/x-icon" href="http://www.yourwebsite.com/favicon.ico" />

% Begin body
fprintf(fid,'<body>\n');
fprintf(fid,'<h1 id="-1">Batch result summary: %s</h1>\n',rdir);

% List simulation directories and generate
content = dir(rdir);
    tmpcell = struct2cell(content);
    [~,ind] = sort_nat(tmpcell(1,:));
    clear tmpcell;
    for d = 1:length(content)
        simrun = content(ind(d)).name;
        if (content(ind(d)).isdir)
            ignoreDir = any(strcmpi(simrun, {'private','CVS','.','..'}));
            ignorePrefix = any(strncmp(simrun, {'@','.'}, 1));
            if ~(ignoreDir || ignorePrefix)
                createSimReport(fid,str2double(simrun)); % Call sub routine
            end
        end
    end

fprintf(fid,'</body>\n');
fprintf(fid,'</html>\n');
fclose(fid);

disp(['Report created in ','<a href="file://',pwd,filesep,'report.html">',rdir,'report.html</a>']);
%system(['ln -s ',rdir,'report.html ',extPath,sdir,'.html']); % Open
%http://mac0.cns.ox.ac.uk/~evans/<sdir>.html instead - put in private
%directory and change the links to js, css and pngs. 
web(['file://',pwd,filesep,'report.html'],'-browser');


function createSimReport(fid,simrun)
% Create hidden simulation reports

if exist([int2str(simrun),filesep,'matlab_workspace.mat'],'file')
    load([int2str(simrun),filesep,'matlab_workspace.mat'],'MP'); %Slow
else
    return
end

fprintf(fid,'<div class="parent" id="parent_%d">\n',simrun);
fprintf(fid,'\t<div class="simTitle" id="%d" onclick="showHide(this.id)">Simulation %d: ',simrun,simrun); % Print CLIparams here
if exist([int2str(simrun),filesep,'CLIparams.m'],'file')% Print CLIparams here
pID = fopen([int2str(simrun),filesep,'CLIparams.m'],'r');
textcell = textscan(pID,'%s','delimiter','\n');
for p=1:length(textcell{1}); fprintf(fid,'%s ',char(textcell{1}(p))); end
fprintf(fid,'</div>\n');
fclose(pID);
else
    fprintf(fid,'[No CLI arguments!] </div>\n');
end
fprintf(fid,'\t<div class="simulation" id="sim_%d" style="display:none;">\n',simrun);

% Information analysis
fprintf(fid,'\t\t<div class="analysis analysis1">\n');
fprintf(fid,'\t\t\t<div class="title">Information analysis</div>\n');
if exist([int2str(simrun),filesep,'infoSingle',int2str(simrun),'.png'],'file')
fprintf(fid,'\t\t\t<div class="plot"><a href="%d/infoSingle%d.png"><img class="png" src="%d/infoSingle%d.png" alt="Single" /></a></div>\n',simrun,simrun,simrun,simrun);
end
if exist([int2str(simrun),filesep,'infoMultiple',int2str(simrun),'.png'],'file')
fprintf(fid,'\t\t\t<div class="plot"><a href="%d/infoMultiple%d.png"><img class="png" src="%d/infoMultiple%d.png" alt="Multiple" /></a></div>\n',simrun,simrun,simrun,simrun);
end
fprintf(fid,'\t\t</div>\n');

% Input rasters
fprintf(fid,'\t\t<div class="analysis analysis2">\n');
fprintf(fid,'\t\t\t<div class="title">Input rasters</div>\n');
if exist([int2str(simrun),filesep,'inputs',int2str(simrun),'.png'],'file')
fprintf(fid,'\t\t\t<div class="plot"><a href="%d/inputs%d.png"><img class="png" src="%d/inputs%d.png" alt="Testing Inputs" /></a></div>\n',simrun,simrun,simrun,simrun);
end
if exist([int2str(simrun),filesep,'E0inputs',int2str(simrun),'.png'],'file')
fprintf(fid,'\t\t\t<div class="plot"><a href="%d/E0inputs%d.png"><img class="png" src="%d/E0inputs%d.png" alt="Training Inputs (Epoch 0)" /></a></div>\n',simrun,simrun,simrun,simrun);
end
fprintf(fid,'\t\t</div>\n');

% Output analysis
fprintf(fid,'\t\t<div class="analysis analysis3">\n');
fprintf(fid,'\t\t\t<div class="title">Output analysis</div>\n');
if exist([int2str(simrun),filesep,'ptOutputAnalysis',int2str(simrun),'.png'],'file')
fprintf(fid,'\t\t\t<div class="plot"><a href="%d/ptOutputAnalysis%d.png"><img class="png" src="%d/ptOutputAnalysis%d.png" alt="PreTraining Output analysis" /></a></div>\n',simrun,simrun,simrun,simrun);
end
if exist([int2str(simrun),filesep,'OutputAnalysis',int2str(simrun),'.png'],'file')
fprintf(fid,'\t\t\t<div class="plot"><a href="%d/OutputAnalysis%d.png"><img class="png" src="%d/OutputAnalysis%d.png" alt="PostTraining Output analysis" /></a></div>\n',simrun,simrun,simrun,simrun);
end
fprintf(fid,'\t\t</div>\n');

if MP.SOM
    fprintf(fid,'\t\t<div class="analysis analysis4">\n');
    fprintf(fid,'\t\t\t<div class="title">SOM output layer</div>\n');
    if exist([int2str(simrun),filesep,'ptL',int2str(MP.nLayers-1),'SOMfRates.png'],'file')
    fprintf(fid,'\t\t\t<div class="plot"><a href="%d/ptL%dSOMfRates.png"><img class="png" src="%d/ptL%dSOMfRates.png" alt="PreTraining SOM outputs" /></a></div>\n',simrun,MP.nLayers-1,simrun,MP.nLayers-1);
    end
    if exist([int2str(simrun),filesep,'L',int2str(MP.nLayers-1),'SOMfRates.png'],'file')
    fprintf(fid,'\t\t\t<div class="plot"><a href="%d/L%dSOMfRates.png"><img class="png" src="%d/L%dSOMfRates.png" alt="PostTraining SOM outputs" /></a></div>\n',simrun,MP.nLayers-1,simrun,MP.nLayers-1);
    end
    fprintf(fid,'\t\t</div>\n');
end

% Close  <div> tags
fprintf(fid,'\t</div>\n'); % Simulation div
fprintf(fid,'</div>\n'); % Parent div