function DrawMonosac(centralLoc, MonoSacType)

pos = [centralLoc(1)-0.4 centralLoc(2)-0.4 0.8 0.8];

if strcmp(MonoSacType,'M')
    rectangle('Position',pos,'Curvature',[1 1],'FaceColor','#77AC30','EdgeColor','#77AC30');
elseif strcmp(MonoSacType,'F')
    patch('XData',[centralLoc(1)-0.4 centralLoc(1)+0.4 centralLoc(1)+0.4],'YData',[centralLoc(2) centralLoc(2)+0.4 centralLoc(2)-0.4],'FaceColor','#A2142F','EdgeColor','#A2142F');
elseif strcmp(MonoSacType,'A')
    rectangle('Position',pos,'Curvature',[1 1],'FaceColor','#EDB120','EdgeColor','#EDB120');
elseif strcmp(MonoSacType,'GN')
    rectangle('Position',pos,'FaceColor','#0072BD','EdgeColor','#0072BD');
elseif strcmp(MonoSacType,'NN')
    displacement = 0.7/sqrt(2);
    patch('XData',[centralLoc(1)-displacement centralLoc(1) centralLoc(1)+displacement centralLoc(1)],...
        'YData',[centralLoc(2) centralLoc(2)+displacement centralLoc(2) centralLoc(2)-displacement],'FaceColor','#7E2F8E','EdgeColor','#7E2F8E');
end

end