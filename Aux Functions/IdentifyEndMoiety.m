function [BranchFlag,EndVec,Branches] = IdentifyEndMoiety(GlycanStr)

%Branch1: GNb2Ma3
%Branch2: GNb4Ma3
%Branch3: GNb2Ma6
%Branch4: GNb6Ma6

% Initiation
BranchFlag = [0 0 0 0];
EndVec = {'','','',''};
idx = strfind(GlycanStr,')Mb4');
GlycanStr = GlycanStr(1:idx+3);
B1 = '';B2 = '';B3 = '';B4 = '';

% Presence of the branches
BranchFlag(4) = any(regexp(GlycanStr,'(GNb6\)Ma6)|(Ma6(\))?Ma6)'));
BranchFlag(2) = any(regexp(GlycanStr,'(?<=GNb2\(\w*GNb4)\)Ma3\(.*Ma6\)Mb4'));
BranchFlag(3) = any(regexp(GlycanStr,'(?<=Ma3\()\w*((Ma3)|(GNb2))(\(\w*\))?Ma6\)Mb4'));
BranchFlag(1) = any(regexp(GlycanStr,'(?<=(GNb2)|(Ma2))(\(\w*GNb4\))?Ma3\(.*Ma6\)Mb4'));


if BranchFlag(1)     
    idx = regexp(GlycanStr,'(?<=(GNb2)|(Ma2))(\(\w*GNb4\))?Ma3\(.*Ma6\)Mb4');
    B1 = GlycanStr(2:idx-1);
    EndVec{1} = B1(1:regexp(B1,'\d','once'));
end

if BranchFlag(2)     
    Endidx = regexp(GlycanStr,'(?<=GNb2\(\w*GNb4)\)Ma3\(.*Ma6\)Mb4');
    Startidx = regexp(GlycanStr,'(?<=GNb2\()\w*GNb4\)Ma3\(.*Ma6\)Mb4');
    B2 = GlycanStr(Startidx:Endidx-1);
    EndVec{2} = B2(1:regexp(B2,'\d','once'));
end

if BranchFlag(3)     
    Startidx = regexp(GlycanStr,'(?<=Ma3\()\w*((Ma3)|(GNb2))(\(\w*\))?Ma6\)Mb4');
    Endidx = regexp(GlycanStr,'(?<=Ma3\(\w*((Ma3)|(GNb2)))(\(\w*\))?Ma6\)Mb4');
    B3 = GlycanStr(Startidx:Endidx-1);
    EndVec{3} = B3(1:regexp(B3,'\d','once'));
end

if BranchFlag(4)     
    Startidx = regexp(GlycanStr,'\(\w*6(\))?Ma6\)Mb4');
    Endidx = regexp(GlycanStr,'(?<=\(\w*6)(\))?Ma6\)Mb4');
    B4 = GlycanStr(Startidx+1:Endidx-1);
    EndVec{4} = B4(1:regexp(B4,'\d','once'));
end

Branches = {B1;B2;B3;B4};

end