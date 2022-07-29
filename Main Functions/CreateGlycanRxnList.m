function [AllrxnList,Glys] = CreateGlycanRxnList(K,RxnSel,targetGlycans)

% set target
if nargin == 2
    targetGlycans = {};
end
targetflag = true;

% Enzyme localization (standard)
Enzymes = {'N/A' 'ManI' 'ManII' 'GnTI' 'GnTII' 'GnTIV' 'GnTV' 'a6FucT' 'b4GalT' 'a3SiaT' 'iGnT' 'a6SiaT'};
EnzLocs = {'N/A' '[cg]' '[mg]' '[cg]' '[mg]' '[mg]' '[tg]' '[mg]' '[tg]' '[tg]' '[tg]' '[tg]'};
[~,idx] = ismember(RxnSel,Enzymes);
Enzymes = Enzymes(idx);EnzLocs = EnzLocs(idx);

% initiate variables
counter = 1;
RxnCounter = 0;
CurrentGlycans = {'(Ma2Ma2Ma3(Ma2Ma3(Ma2Ma6)Ma6)Mb4GNb4GN);Asn'};
ProcessedGlycans = {};

% Initiate output variable
Rxns = cell(1,6); % column 1: reactant ; column 2: product; column 3: gene 4: compartments 5: complexity level 6. Rxn Names

% Wait Bar
fprintf('Construct a generic N-glycosylation network...\n');
f1 = waitbar(0,['Add N-Glycotransferases & Glycosidases Reactions (Complexity Level ',num2str(counter),') ...']);

while counter <= K && targetflag
    
    % Get Current length of Network.mets
    
    % Update Wait Bar
    waitbar((counter-1)/K,f1,['Step1: Add N-Glycotransferases & Glycosidases Reactions (Complexity Level ',num2str(counter),') ...'])
    
    %% ManI
    %----------------------------------------------------
    % Identify glycotransferase features
    Enzyme = 'ManI';
    if any(strcmp(Enzyme,Enzymes))
        SubSystem = EnzLocs{strcmp(Enzymes,Enzyme)};
        SubstrateString = '\(Ma2Ma';SubstrateString_alt = '(Ma2Ma';
        ProductString = 'Ma';
        
        % get corresponding substrates and add reactions
        Reactants_temp = CurrentGlycans(cellfun(@(x) contains(x,SubstrateString_alt),CurrentGlycans));
        
        if ~isempty(Reactants_temp)
            for k = 1:length(Reactants_temp)
                [startIndex,endIndex] = regexp(Reactants_temp{k},SubstrateString);
                Reactants = repmat(Reactants_temp(k),length(startIndex),1);
                Products = Reactants;
                for k1 = 1:length(startIndex)
                    Products{k1} = [Products{k1}(1:startIndex(k1)),ProductString,Products{k1}((endIndex(k1)+1):end)];
                end
                
                % document new glycan conversion
                for a = 1:length(Reactants)
                    % add new glycotransferase reaction
                    RxnCounter = RxnCounter + 1;
                    Rxns{RxnCounter,1} = Reactants{a};Rxns{RxnCounter,2} = Products{a};Rxns{RxnCounter,3} = Enzyme;Rxns{RxnCounter,4} = SubSystem;Rxns{RxnCounter,5} = counter;Rxns{RxnCounter,6} = ['GlyRxn_',num2str(RxnCounter)];
                end
            end
        end
    end
    
    %% GnTI
    %----------------------------------------------------
    % Identify glycotransferase features
    Enzyme = 'GnTI';
    if any(strcmp(Enzyme,Enzymes))
        SubSystem = EnzLocs{strcmp(Enzymes,Enzyme)};
        SubstrateString = '\(Ma3\(Ma3\(Ma6\)Ma6\)Mb4';SubstrateString_alt = '(Ma3(Ma3(Ma6)Ma6)Mb4';
        ProductString = 'GNb2Ma3(Ma3(Ma6)Ma6)Mb4';
        
        % get corresponding substrates and add reactions
        Reactants_temp = CurrentGlycans(cellfun(@(x) contains(x,SubstrateString_alt),CurrentGlycans));
        
        if ~isempty(Reactants_temp)
            for k = 1:length(Reactants_temp)
                [startIndex,endIndex] = regexp(Reactants_temp{k},SubstrateString);
                Reactants = repmat(Reactants_temp(k),length(startIndex),1);
                Products = Reactants;
                for k1 = 1:length(startIndex)
                    Products{k1} = [Products{k1}(1:startIndex(k1)),ProductString,Products{k1}((endIndex(k1)+1):end)];
                end
                
                
                % document new glycan conversion
                for a = 1:length(Reactants)
                    % add new glycotransferase reaction
                    RxnCounter = RxnCounter + 1;
                    Rxns{RxnCounter,1} = Reactants{a};Rxns{RxnCounter,2} = Products{a};Rxns{RxnCounter,3} = Enzyme;Rxns{RxnCounter,4} = SubSystem;Rxns{RxnCounter,5} = counter;Rxns{RxnCounter,6} = ['GlyRxn_',num2str(RxnCounter)];
                end
            end
        end
    end
    
    
    %% ManII_1
    %----------------------------------------------------
    % Identify glycotransferase features
    Enzyme = 'ManII';
    if any(strcmp(Enzyme,Enzymes))
        SubSystem = EnzLocs{strcmp(Enzymes,Enzyme)};
        SubstrateString = '\(Ma3\(Ma6\)Ma6'; SubstrateString_alt = '(Ma3(Ma6)Ma6';
        ProductString = 'Ma6Ma6';
        
        % Identify constraints
        contraints = regexp(CurrentGlycans,'(\(GNb2Ma3)|(\(GNb2\(\w+\)Ma3)', 'once');
        
        % get corresponding substrates and add reactions
        Reactants_temp = CurrentGlycans(cellfun(@(x) contains(x,SubstrateString_alt), CurrentGlycans) & ~cellfun(@isempty,contraints));
        
        if ~isempty(Reactants_temp)
            for k = 1:length(Reactants_temp)
                [startIndex,endIndex] = regexp(Reactants_temp{k},SubstrateString);
                Reactants = repmat(Reactants_temp(k),length(startIndex),1);
                Products = Reactants;
                for k1 = 1:length(startIndex)
                    Products{k1} = [Products{k1}(1:startIndex(k1)),ProductString,Products{k1}((endIndex(k1)+1):end)];
                end
                
                
                % document new glycan conversion
                for a = 1:length(Reactants)
                    % add new glycotransferase reaction
                    RxnCounter = RxnCounter + 1;
                    Rxns{RxnCounter,1} = Reactants{a};Rxns{RxnCounter,2} = Products{a};Rxns{RxnCounter,3} = Enzyme;Rxns{RxnCounter,4} = SubSystem;Rxns{RxnCounter,5} = counter;Rxns{RxnCounter,6} = ['GlyRxn_',num2str(RxnCounter)];
                end
            end
        end
    end
    
    
    %% ManII_2
    %----------------------------------------------------
    % Identify glycotransferase features
    Enzyme = 'ManII';
    if any(strcmp(Enzyme,Enzymes))
        SubSystem = EnzLocs{strcmp(Enzymes,Enzyme)};
        SubstrateString = '\(Ma6Ma6'; SubstrateString_alt = '(Ma6Ma6';
        ProductString = 'Ma6';
        
        % Identify constraints
        contraints = regexp(CurrentGlycans,'(\(GNb2Ma3)|(\(GNb2\(\w+\)Ma3)', 'once');
        
        % get corresponding substrates and add reactions
        Reactants_temp = CurrentGlycans(cellfun(@(x) contains(x,SubstrateString_alt), CurrentGlycans) & ~cellfun(@isempty,contraints));
        
        if ~isempty(Reactants_temp)
            for k = 1:length(Reactants_temp)
                [startIndex,endIndex] = regexp(Reactants_temp{k},SubstrateString);
                Reactants = repmat(Reactants_temp(k),length(startIndex),1);
                Products = Reactants;
                for k1 = 1:length(startIndex)
                    Products{k1} = [Products{k1}(1:startIndex(k1)),ProductString,Products{k1}((endIndex(k1)+1):end)];
                end
                
                
                % document new glycan conversion
                for a = 1:length(Reactants)
                    % add new glycotransferase reaction
                    RxnCounter = RxnCounter + 1;
                    Rxns{RxnCounter,1} = Reactants{a};Rxns{RxnCounter,2} = Products{a};Rxns{RxnCounter,3} = Enzyme;Rxns{RxnCounter,4} = SubSystem;Rxns{RxnCounter,5} = counter;Rxns{RxnCounter,6} = ['GlyRxn_',num2str(RxnCounter)];
                end
            end
        end
    end
    
    %% GnTII
    
    % Identify glycotransferase features
    Enzyme = 'GnTII';
    if any(strcmp(Enzyme,Enzymes))
        SubSystem = EnzLocs{strcmp(Enzymes,Enzyme)};
        SubstrateString = '((\(GNb2)|(\(GNb2\(\w+\)))Ma3\(Ma6\)Mb4'; SubstrateString_alt = '((\(GNb2)|(\(GNb2\(\w+\)))Ma3\(Ma6\)Mb4';
        ProductString = 'GNb2Ma6)Mb4';
        
        % get corresponding substrates and add reactions
        Reactants_temp = CurrentGlycans(cellfun(@(x) ~isempty(regexp(x,SubstrateString_alt,'once')),CurrentGlycans));
        
        if ~isempty(Reactants_temp)
            for k = 1:length(Reactants_temp)
                [~,endIndex] = regexp(Reactants_temp{k},SubstrateString);
                Reactants = repmat(Reactants_temp(k),length(startIndex),1);
                Products = Reactants;
                for k1 = 1:length(startIndex)
                    Products{k1} = [Products{k1}(1:(endIndex(k1)-7)),ProductString,Products{k1}((endIndex(k1)+1):end)];
                end
                
                
                % document new glycan conversion
                for a = 1:length(Reactants)
                    % add new glycotransferase reaction
                    RxnCounter = RxnCounter + 1;
                    Rxns{RxnCounter,1} = Reactants{a};Rxns{RxnCounter,2} = Products{a};Rxns{RxnCounter,3} = Enzyme;Rxns{RxnCounter,4} = SubSystem;Rxns{RxnCounter,5} = counter;Rxns{RxnCounter,6} = ['GlyRxn_',num2str(RxnCounter)];
                end
            end
        end
    end
    
    %% a6FucT
    %----------------------------------------------------
    % Identify glycotransferase features
    Enzyme = 'a6FucT';
    if any(strcmp(Enzyme,Enzymes))
        SubSystem = EnzLocs{strcmp(Enzymes,Enzyme)};
        SubstrateString = 'GNb4GN)'; SubstrateString_alt = 'GNb4GN)';
        ProductString = 'GNb4(Fa6)GN)';
        
        % get corresponding substrates and add reactions
        Reactants_temp = CurrentGlycans(cellfun(@(x) contains(x,SubstrateString_alt), CurrentGlycans));
        
        if ~isempty(Reactants_temp)
            for k = 1:length(Reactants_temp)
                [startIndex,endIndex] = regexp(Reactants_temp{k},SubstrateString);
                Reactants = repmat(Reactants_temp(k),length(startIndex),1);
                Products = Reactants;
                for k1 = 1:length(startIndex)
                    Products{k1} = [Products{k1}(1:(startIndex(k1)-1)),ProductString,Products{k1}((endIndex(k1)+1):end)];
                end
                
                
                % document new glycan conversion
                for a = 1:length(Reactants)
                    % add new glycotransferase reaction
                    RxnCounter = RxnCounter + 1;
                    Rxns{RxnCounter,1} = Reactants{a};Rxns{RxnCounter,2} = Products{a};Rxns{RxnCounter,3} = Enzyme;Rxns{RxnCounter,4} = SubSystem;Rxns{RxnCounter,5} = counter;Rxns{RxnCounter,6} = ['GlyRxn_',num2str(RxnCounter)];
                end
            end
        end
    end
    
    
    %% GnTIV
    %----------------------------------------------------
    % Identify glycotransferase features
    Enzyme = 'GnTIV';
    if any(strcmp(Enzyme,Enzymes))
        SubSystem = EnzLocs{strcmp(Enzymes,Enzyme)};
        SubstrateString = '\(GNb2Ma3'; SubstrateString_alt = '(GNb2Ma3';
        ProductString = 'GNb2(GNb4)Ma3';
        
        % get corresponding substrates and add reactions
        Reactants_temp = CurrentGlycans(cellfun(@(x) contains(x,SubstrateString_alt),CurrentGlycans));
        
        if ~isempty(Reactants_temp)
            for k = 1:length(Reactants_temp)
                [startIndex,endIndex] = regexp(Reactants_temp{k},SubstrateString);
                Reactants = repmat(Reactants_temp(k),length(startIndex),1);
                Products = Reactants;
                for k1 = 1:length(startIndex)
                    Products{k1} = [Products{k1}(1:startIndex(k1)),ProductString,Products{k1}((endIndex(k1)+1):end)];
                end
                
                
                % document new glycan conversion
                for a = 1:length(Reactants)
                    % add new glycotransferase reaction
                    RxnCounter = RxnCounter + 1;
                    Rxns{RxnCounter,1} = Reactants{a};Rxns{RxnCounter,2} = Products{a};Rxns{RxnCounter,3} = Enzyme;Rxns{RxnCounter,4} = SubSystem;Rxns{RxnCounter,5} = counter;Rxns{RxnCounter,6} = ['GlyRxn_',num2str(RxnCounter)];
                end
            end
        end
    end
    
    
    %% GnTV
    %----------------------------------------------------
    % Identify glycotransferase features
    Enzyme = 'GnTV';
    if any(strcmp(Enzyme,Enzymes))
        SubSystem = EnzLocs{strcmp(Enzymes,Enzyme)};
        SubstrateString = '\(GNb2Ma6'; SubstrateString_alt = '(GNb2Ma6';
        ProductString = 'GNb2(GNb6)Ma6';
        
        % get corresponding substrates and add reactions
        Reactants_temp = CurrentGlycans(cellfun(@(x) contains(x,SubstrateString_alt),CurrentGlycans));
        
        if ~isempty(Reactants_temp)
            for k = 1:length(Reactants_temp)
                [startIndex,endIndex] = regexp(Reactants_temp{k},SubstrateString);
                Reactants = repmat(Reactants_temp(k),length(startIndex),1);
                Products = Reactants;
                for k1 = 1:length(startIndex)
                    Products{k1} = [Products{k1}(1:startIndex(k1)),ProductString,Products{k1}((endIndex(k1)+1):end)];
                end
                
                % document new glycan conversion
                for a = 1:length(Reactants)
                    % add new glycotransferase reaction
                    RxnCounter = RxnCounter + 1;
                    Rxns{RxnCounter,1} = Reactants{a};Rxns{RxnCounter,2} = Products{a};Rxns{RxnCounter,3} = Enzyme;Rxns{RxnCounter,4} = SubSystem;Rxns{RxnCounter,5} = counter;Rxns{RxnCounter,6} = ['GlyRxn_',num2str(RxnCounter)];
                end
            end
        end
    end
    
    %% b4GalT
    %----------------------------------------------------
    % Identify glycotransferase features
    Enzyme = 'b4GalT';
    if any(strcmp(Enzyme,Enzymes))
        SubSystem = EnzLocs{strcmp(Enzymes,Enzyme)};
        SubstrateString = '\(GN';SubstrateString_alt = '(GN';
        ProductString = 'Ab4GN';
        
        % get corresponding substrates and add reactions
        Reactants_temp= CurrentGlycans(cellfun(@(x) contains(x,SubstrateString_alt),CurrentGlycans));
        
        if ~isempty(Reactants_temp)
            for k = 1:length(Reactants_temp)
                [startIndex,endIndex] = regexp(Reactants_temp{k},SubstrateString);
                Reactants = repmat(Reactants_temp(k),length(startIndex),1);
                Products = Reactants;
                for k1 = 1:length(startIndex)
                    Products{k1} = [Products{k1}(1:startIndex(k1)),ProductString,Products{k1}((endIndex(k1)+1):end)];
                end
                
                % document new glycan conversion
                for a = 1:length(Reactants)
                    % add new glycotransferase reaction
                    RxnCounter = RxnCounter + 1;
                    Rxns{RxnCounter,1} = Reactants{a};Rxns{RxnCounter,2} = Products{a};Rxns{RxnCounter,3} = Enzyme;Rxns{RxnCounter,4} = SubSystem;Rxns{RxnCounter,5} = counter;Rxns{RxnCounter,6} = ['GlyRxn_',num2str(RxnCounter)];
                end
            end
        end
    end
    
    
    %% iGnT
    % Identify glycotransferase features
    Enzyme = 'iGnT';
    if any(strcmp(Enzyme,Enzymes))
        SubSystem = EnzLocs{strcmp(Enzymes,Enzyme)};
        SubstrateString = '\(Ab4GN'; SubstrateString_alt = '(Ab4GN';
        ProductString = 'GNb3Ab4GN';
        
        % get corresponding substrates and add reactions
        Reactants_temp = CurrentGlycans(cellfun(@(x) contains(x,SubstrateString_alt),CurrentGlycans));
        
        if ~isempty(Reactants_temp)
            for k = 1:length(Reactants_temp)
                [startIndex,endIndex] = regexp(Reactants_temp{k},SubstrateString);
                
                Reactants = repmat(Reactants_temp(k),length(startIndex),1);
                Products = Reactants;
                for k1 = 1:length(startIndex)
                    Products{k1} = [Products{k1}(1:startIndex(k1)),ProductString,Products{k1}((endIndex(k1)+1):end)];
                end
                
                % add constraint
                flag = cellfun(@(x) isempty(regexp(x,'(GNb3\w*(\))?Ma3)|(GNb3\w*\(\w*\)Ma3)','once')),Products);
                Reactants = Reactants(flag); Products = Products(flag);
                
                % document new glycan conversion
                for a = 1:length(Reactants)
                    % add new glycotransferase reaction
                    RxnCounter = RxnCounter + 1;
                    Rxns{RxnCounter,1} = Reactants{a};Rxns{RxnCounter,2} = Products{a};Rxns{RxnCounter,3} = Enzyme;Rxns{RxnCounter,4} = SubSystem;Rxns{RxnCounter,5} = counter;Rxns{RxnCounter,6} = ['GlyRxn_',num2str(RxnCounter)];
                end
                
                
            end
        end
    end
    
    
    %% a3SiaT
    % Identify glycotransferase features
    Enzyme = 'a3SiaT';
    if any(strcmp(Enzyme,Enzymes))
        SubSystem = EnzLocs{strcmp(Enzymes,Enzyme)};
        SubstrateString = '\(Ab4GN';SubstrateString_alt = '(Ab4GN';
        ProductString = 'NNa3Ab4GN';
        
        % get corresponding substrates and add reactions
        Reactants_temp = CurrentGlycans(cellfun(@(x) contains(x,SubstrateString_alt),CurrentGlycans));
        
        if ~isempty(Reactants_temp)
            for k = 1:length(Reactants_temp)
                [startIndex,endIndex] = regexp(Reactants_temp{k},SubstrateString);
                Reactants = repmat(Reactants_temp(k),length(startIndex),1);
                Products = Reactants;
                for k1 = 1:length(startIndex)
                    Products{k1} = [Products{k1}(1:startIndex(k1)),ProductString,Products{k1}((endIndex(k1)+1):end)];
                end
                
                % document new glycan conversion
                for a = 1:length(Reactants)
                    % add new glycotransferase reaction
                    RxnCounter = RxnCounter + 1;
                    Rxns{RxnCounter,1} = Reactants{a};Rxns{RxnCounter,2} = Products{a};Rxns{RxnCounter,3} = Enzyme;Rxns{RxnCounter,4} = SubSystem;Rxns{RxnCounter,5} = counter;Rxns{RxnCounter,6} = ['GlyRxn_',num2str(RxnCounter)];
                end
            end
        end
    end
    
    %% a6SiaT
    % Identify glycotransferase features
    Enzyme = 'a6SiaT';
    if any(strcmp(Enzyme,Enzymes))
        SubSystem = EnzLocs{strcmp(Enzymes,Enzyme)};
        SubstrateString = '\(Ab4GN';SubstrateString_alt = '(Ab4GN';
        ProductString = 'NNa6Ab4GN';
        
        % get corresponding substrates and add reactions
        Reactants_temp = CurrentGlycans(cellfun(@(x) contains(x,SubstrateString_alt),CurrentGlycans));
        
        if ~isempty(Reactants_temp)
            for k = 1:length(Reactants_temp)
                [startIndex,endIndex] = regexp(Reactants_temp{k},SubstrateString);
                Reactants = repmat(Reactants_temp(k),length(startIndex),1);
                Products = Reactants;
                for k1 = 1:length(startIndex)
                    Products{k1} = [Products{k1}(1:startIndex(k1)),ProductString,Products{k1}((endIndex(k1)+1):end)];
                end
                
                % document new glycan conversion
                for a = 1:length(Reactants)
                    % add new glycotransferase reaction
                    RxnCounter = RxnCounter + 1;
                    Rxns{RxnCounter,1} = Reactants{a};Rxns{RxnCounter,2} = Products{a};Rxns{RxnCounter,3} = Enzyme;Rxns{RxnCounter,4} = SubSystem;Rxns{RxnCounter,5} = counter;Rxns{RxnCounter,6} = ['GlyRxn_',num2str(RxnCounter)];
                end
            end
        end
    end
    
    
    %% Process Update
    % Show processing time for each loop
    %     duration = toc(timerVal);
    %     fprintf(['Construction time for Complexity Level ',num2str(counter),' is: ',num2str(duration),'\n']);
    
    
    % Update Current Glycans
    ProcessedGlycans = union(ProcessedGlycans,CurrentGlycans);
    CurrentGlycans = setdiff(Rxns(:,2),ProcessedGlycans);
    
    % Update Counter
    counter = counter+1;
    
    %% Check if target Glycans Achieved
    if ~isempty(targetGlycans) && all(ismember(targetGlycans,ProcessedGlycans))
        targetflag = false;
    end
    
end

% close waitbar
delete(f1);

%% Remove backfeed reactions
f2 = waitbar(0,'Step 2: Trim backfeed reactions...');
% Step 2: Trim backfeed reactions .

Rxns = [Rxns,cell(size(Rxns,1),1)];
RxnSize = size(Rxns,1);
for a = 1:RxnSize
    
    if mod(a/RxnSize,0.1) == 0
        waitbar(a/RxnSize,f2);
    end
    
    Rct = Rxns{a,1};
    rxnCmp = Rxns{a,4};
    
    cmp = IdentifyBirthComp(Rct);
    if (strcmp(cmp,'mg') && ismember(rxnCmp,{'[mg]','[tg]'})) || ...
            (strcmp(cmp,'tg') && ismember(rxnCmp,{'[tg]'})) || ...
            (strcmp(cmp,'cg'))
        % not backfed
        Rxns{a,7} = 1;
    else
        % backfed
        Rxns{a,7} = 0;
    end
    
end

Rxns = Rxns(logical(cell2mat(Rxns(:,7))),1:6);

delete(f2);

%% Add Transport Reactions
f3 = waitbar(0,'Step 3: Add Golgi-intercompartmental transport ...');

Glys = unique(union(unique(Rxns(:,1)),unique(Rxns(:,2))));
GlysSize = length(Glys);

% Add transportation to reaction list
transportlist = cell(1,6);
counter = 1;

for a = 1:GlysSize
    
    if mod(a/GlysSize,0.1) == 0
        waitbar(a/GlysSize,f3);
    end
    
    % identify first compartment appearance
    cmp = IdentifyBirthComp(Glys{a});
    transportlist(counter,[1,2,3,4,6]) = {Glys{a},Glys{a},'tg2ab','[tg]/[ab]',['tg2ab_',num2str(a)]};
    counter = counter + 1;
    if ~strcmp(cmp,'tg')
        transportlist(counter,[1,2,3,4,6]) = {Glys{a},Glys{a},'mg2tg','[mg]/[tg]',['mg2tg_',num2str(a)]};
        counter = counter + 1;
        if strcmp(cmp,'cg')
            transportlist(counter,[1,2,3,4,6]) = {Glys{a},Glys{a},'cg2mg','[cg]/[mg]',['cg2mg_',num2str(a)]};
            counter = counter + 1;
        end
    end
end

AllrxnList = [Rxns;transportlist];

% remove syntheses of unrealistic glycans
% Rule 1: with excessively long antenna (GNb3 >4)
for a = 1:size(AllrxnList,1)
    [~,~,Branches] = IdentifyEndMoiety(AllrxnList{a,2});
    freq = cellfun(@(x) length(strfind(x,'GNb3')),Branches);
    flag(a) = all(freq<5);
end
AllrxnList = AllrxnList(flag,:);

% Obtain the list of all generated glycans
counter = 1;
for a = 1:size(AllrxnList,1)
    comp = AllrxnList{a,4};
    if contains(comp,'/')
        Glys{counter} = [AllrxnList{a,1},comp(1:4)];
        Glys{counter+1} = [AllrxnList{a,1},comp(6:end)];
    else
        Glys{counter} = [AllrxnList{a,1},comp];
        Glys{counter+1} = [AllrxnList{a,2},comp];
    end
    counter = counter + 2;
end

Glys = unique(Glys);

% reorganize glycan list
Glys_abs = Glys(contains(Glys,'[ab]'));
Glys_real = Glys(~contains(Glys,'[ab]'));
Glys = [Glys_real;Glys_abs];

delete(f3);
end