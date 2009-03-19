function raster(spikes,TMAX,DT)
%for f=1:(size(spikes,1)/nps)
    %figure();
    %TMAX=T*DT*1000;
    for p=1:size(spikes,1)%nps
        s = 2; % Could set back to 1 and remove the first column of spikes
        %subplot(nps,1,p);
        hold on;
        if mod(p,2)==0
            while (spikes(p,s) ~= 0)
                plot([(spikes(p,s)*DT*1000),(spikes(p,s)*DT*1000)],[p-0.5,p+0.5],'-b');
                plot([0,TMAX],[p,p],':k');
                s = s + 1;
            end
        else
            while (spikes(p,s) ~= 0)
                plot([(spikes(p,s)*DT*1000),(spikes(p,s)*DT*1000)],[p-0.5,p+0.5],'-r');
                plot([0,TMAX],[p,p],':k');
                s = s + 1;
            end
        end
        % axes('position',[left,bottom,width,height]);
        %axis off;
    end
    title('Raster plots of output layer');
    xlabel('Time (msec)');
    ylabel('Cell Number');
    %ylabel(['Neuron=',int2str(((f-1)*nps)+p)]);
    xlim([1,TMAX]);
    ylim([0.5,size(spikes,1)+0.5]);
    hold off;
    %set(gca, 'YGrid','on','YMinorGrid','on');
    %box off;
    %axes('position',[0.1,0.1,0.9,0.9]);
%end
