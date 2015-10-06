function fname = saveFig(fh, fname, plotLevel, opt) %fh, method) %plotLevel

% Replace strcat() with []

% For an image, you can call imshow with the following option:
% >> imshow(<image>, 'Border', 'tight');

% To expand a figure to fill an entire window (remove gray border & labels):
% >> set(gca, 'Position', [0 0 1 1]);
% To crop the figure tight up to the labels and title (eps2pdf does this): 
% set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
% http://nibot-lab.livejournal.com/73290.html

% Tips: http://blogs.mathworks.com/loren/2007/12/11/making-pretty-graphs/
% \usepackage{epsfig}: http://andrewjpage.com/index.php?/archives/27-Using-an-EPS-in-Latex.html
if nargin < 4
    % Publications require high quality 800DPI eps
    opt.format = 'eps2';
    opt.dpi = 800;
    %opt.width = 6;
    opt.FontMode = 'scaled';
    opt.LineMode = 'scaled';
    %opt.fontmode = 'fixed';
    %opt.fontsize = 8;
    opt.color = 'cmyk';
    opt.renderer = 'painters'; %[];
    args = {};
    
    if nargin < 3
        plotLevel = 1; % PNG
        if nargin < 2
            error('Insufficient arguments');
        end
    end
else
    if isfield(opt,'FontSize') || isfield(opt,'fontsize')
        opt.FontMode = 'fixed';
    end
    if isfield(opt,'LineMode') || isfield(opt,'linemode')
        opt.LineMode = 'fixed';
    end
end

gspath='/usr/local/bin/gs';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plotLevel > 0
    saveas(fh,fname,'png');
    if plotLevel > 1
        saveas(fh, [fname,'-orig'], 'epsc2');
        eps2pdf([fname,'-orig','.eps'],gspath);
        
        if (~isempty(opt.renderer))
            args = {args{:}, 'renderer', opt.renderer};
        end
        exportfig(fh, fname, 'Format',opt.format, 'Resolution',opt.dpi, ...
            'Color',opt.color, 'FontMode',opt.FontMode, 'LineMode',opt.LineMode, args{:});
        eps2pdf([fname,'.eps'],gspath);
    end
end


% N.B. exportfig has been modified so that 'Color', 'cmyk' simply enables
% colour for the figure, (not cmyk instead of the default rgb) - line 282


% savePDF=true; 
% %%%%%%%%%%%%%%%%%%%%%%%%
% 
% if nargin < 2
%     % Default method
%     method = 0;
%     %ah = get(gcf,'CurrentAxes'); % Used for getting axes within a figure
%     fh = get(0,'CurrentFigure'); % eqivalent to fh = gcf; but will not create a figure if one does not exist
% end


% %if plotLevel > 0
% if method%plotLevel > 1 %if (savePDF)
%     saveas(fh,['[M1]',fname],'epsc');
%     eps2pdf(['[M1]',fname,'.eps'],gspath);
% end
% saveas(fh,fname,'png');
% %end
% 
% % Replace epsc wih epsc2 in saveas commands 
% saveas(fh,['[M1b]',fname],'epsc2');
% 
% 
% % Use exportfig() (in ~/Matlab/figures/) N.B. Lots of parameters to play
% % with and used at fMRIB: http://www.mathworks.co.uk/company/newsletters/digest/june00/export/index.html http://www.mathworks.co.uk/company/newsletters/digest/december00/export.html
% %  e.g. exportfig(gcf,'test.eps','color','rgb','resolution',300)
% exportfig(fh,['[M2]',fname],'resolution',dpi);
% % Other companion functions available
% 
% 
% % Compare print() to saveas() http://cvlab.epfl.ch/~ksmith/tips.html
% print(['-f',int2str(fh)],'-depsc2',['-r',int2str(dpi)],['[M3]',fname]);
% % -noui
% % Renderers: -painters, -zbuffer, -opengl
% 
% % export_fig 
% addpath('~/Matlab/figures/export_fig/');
% export_fig(['[M4]',fname], fh, ['-r',int2str(dpi)], '-transparent', '-eps', '-pdf');
% 
% % savefig
% addpath('~/Matlab/figures/savefig/');
% savefig(['[M5]',fname], fh, ['-r',int2str(dpi)], '-lossless', 'eps');
% savefig(['[M5]',fname], fh, ['-r',int2str(dpi)], '-lossless', 'pdf');
% 
% 
% % EPS and TEX file outputs
% % Try laprint for latex: http://www.uni-kassel.de/fb16/rat/matlab/laprint/
% 
% Matlabfrag
% % http://www.mathworks.com/matlabcentral/fileexchange/21286-matlabfrag
% addpath('~/Matlab/figures/matlabfrag/');




% http://win.ua.ac.be/~nschloe/content/matlab2tikz [Matlab -> LaTeX plots]

% http://pgfplots.sourceforge.net/ [LaTeX package]
