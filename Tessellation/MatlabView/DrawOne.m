% Created on 13.8.2017
% author: Milan Pultar, email: milan.pultar@gmail.com
% script for a visualisation of the Laguerre tessellation in a special format

cc=hsv(12);
C = importdata(strcat('..\test10uncon.txt'));
% change to whatever directory the 'out.txt' file is in
h = figure;
ct = 1;
zrno=1;
for i=1:size(C,1)
    if ct==2
        % This section displays centers of grains
        pp = scatter3(C(i,1), C(i,2),C(i,3),'go');
        hold on;
        scatter3(C(i,1), C(i,2),C(i,3),'r+');
        hold on;
        currentunits = get(gca,'Units');
        set(gca, 'Units', 'Points');
        axpos = get(gca,'Position');
        set(gca, 'Units', currentunits);
        markerWidth = C(i,4)/diff(xlim)*axpos(3); % Calculate Marker width in points
        set(pp, 'SizeData', markerWidth^2)
    elseif not(isnan(C(i,1)))
        %plot3([C(i,1) C(i,4)],[C(i,2) C(i,5)],[C(i,3) C(i,6)],'Color',cc(mod(zrno,12),:));
        % change to the line above if you want each grain drawn by
        % different color
        plot3([C(i,1) C(i,4)],[C(i,2) C(i,5)],[C(i,3) C(i,6)],'Color','b');
        hold on;
        if ct==1
            ct=0;
        end
    else
        ct=ct+1;
        zrno=zrno+1;
        if ct==2
            zrno=zrno-2;
        end
    end
end
disp('Pocet zrn: ');
disp(zrno);
