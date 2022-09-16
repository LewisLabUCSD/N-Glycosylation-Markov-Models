function DrawGlycanStructure(str,filename,mzVals)

if ischar(str)
    str = cellstr(str);
end
strs = str;

if nargin==1
    filename = [];
    mzVals = 1:length(strs);
end

plotWidth = ceil(length(strs)/5);

% Initiate the plot
h = figure;
t = tiledlayout(plotWidth,5,'TileSpacing','none','Padding','tight');


for k = 1:length(strs)

    nexttile
    hold on
    str = strs{k};

    % Process and Decompose linear code of the branches
    [BranchFlag,~,Branches] = IdentifyEndMoiety(str);
    str = strrep(str,';Asn','');
    Fut8Flag = contains(str,'GNb4(Fa6)GN');



    % Plot root glycan
    DrawMonosac([0 0], 'GN');
    DrawMonosac([0 1], 'GN');
    if Fut8Flag
        DrawMonosac([1 0], 'F');
    end
    DrawMonosac([0 2], 'M');
    DrawMonosac([-1 2.7], 'M');
    DrawMonosac([1 2.7], 'M');

    % Define the coordinates based on the number of branches
    BranchXPos = [-1.5 -0.5 0.5 1.5];
    if (BranchFlag(1) && ~BranchFlag(2)) || (BranchFlag(2) && ~BranchFlag(1))
        BranchXPos([1 2]) = -1;
    end
    if (BranchFlag(3) && ~BranchFlag(4)) || (BranchFlag(4) && ~BranchFlag(3))
        BranchXPos([3 4]) = 1;
    end
    YPosStart = 3.7;

    % Plot each branches
    for a = 1:length(BranchFlag)
        if BranchFlag(a)
            YPosStart_temp = YPosStart;
            linecode = Branches{a};
            [startIdx,endIdx] = regexp(linecode,'(A)|(GN)|(NN)|(F)|(M)');
            for b = 0:length(startIdx)-1
                MonoSacType = linecode(startIdx(end-b):endIdx(end-b));
                DrawMonosac([BranchXPos(a) YPosStart_temp], MonoSacType);
                YPosStart_temp = YPosStart_temp+1;
            end
        end
    end

    % set visuals
    globmax = max(abs([xlim ylim]));
    xlim([-globmax globmax]);
    ylim([-globmax globmax]);
    axis square
    axis off
    % set(gca,'visible','off');

    if nargin == 3
        title(sprintf('%0.1f',mzVals(k)));
    end

    hold off
end

h=axes(h,'visible','off'); 
h.XLabel.Visible='on';

% save as png
if ~isempty(filename)
xlabel(h, ['Major Glycoform at annotated m/z for ',strrep(filename,'_','/')],'FontWeight','bold');
saveas(gcf,['Figures/GlycanFigures/',filename,'.tif']);
else
    xlabel(h,'Glycoforms');
end

end