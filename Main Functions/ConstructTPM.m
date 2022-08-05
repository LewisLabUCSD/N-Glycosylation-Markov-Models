function  [TM, AllrxnList_TMidx, AllrxnList_Br, AllrxnList_RxnTypes, RxnTypes,AllrxnList_steric,stericRxns,AllrxnList_LacNAcLen,AllrxnList_LacNAcLen_idx] = ConstructTPM(AllrxnList,Glys)

fprintf('Construct Transition Probability Matrix from the generic N-glycosylation network...\n');

N = size(AllrxnList,1);
WaitMessage = parfor_wait(N, 'Waitbar', true);
AllrxnList_TMrow = zeros(size(AllrxnList,1),1);
AllrxnList_TMcol = AllrxnList_TMrow;
AllrxnList_Br = AllrxnList_TMrow;
AllrxnList_steric = AllrxnList_TMrow;
AllrxnList_LacNAcLen = AllrxnList_TMrow;

% List reactions to consider steric interactions
stericRxns =  {'a3SiaT','b4GalT','iGnT','GnTIV','GnTV'};
% List reactions to consider branching
branchingRxns = {'a3SiaT','b4GalT','iGnT'};

% Identify transition probability matrix indice for each reaction
parfor k = 1:size(AllrxnList,1)

    WaitMessage.Send;

    TempSel = AllrxnList(k,:);
    Rct = TempSel{1};
    Pro = TempSel{2};
    RxnType = TempSel{3};
    comp1 = TempSel{4};
    if ~contains(comp1,'/')
        comp2 = comp1;
    else
        comp2 = comp1(6:9);
        comp1 = comp1(1:4);
    end

    AllrxnList_TMrow(k) = find(strcmp(Glys,[Rct,comp1]),1);
    AllrxnList_TMcol(k) = find(strcmp(Glys,[Pro,comp2]),1);

    if strcmp(RxnType,'iGnT')
        [~,~,Rctbranches] = IdentifyEndMoiety(Rct);
        [~,~,Probranches] = IdentifyEndMoiety(Pro);
        AllrxnList_LacNAcLen(k) = length(strfind(Rctbranches{~strcmp(Rctbranches,Probranches)},'Ab4GNb3'));
    end

    % Distinguish branches and identify steric hindrance for
    % 'a3SiaT','b4GalT', 'iGnT','GnTIV', & 'GnTV')
    if ismember(RxnType,stericRxns)
        [~,~,Rctbranches] = IdentifyEndMoiety(Rct);
        [~,~,Probranches] = IdentifyEndMoiety(Pro);
        loc = find(~strcmp(Rctbranches,Probranches));

        AllrxnList_Br(k) = loc;

        % a reactant site is called "shielded" if the addition site is on a
        % branch where the neighboring branches are longer than the branch where the addition site is on.
        neighborIdx = intersect(1:4,[loc-1 loc+1]);
        shieldflag = false(1,length(neighborIdx));
        addbranchLen = length(strfind(Rctbranches{loc},'a')) + length(strfind(Rctbranches{loc},'b'));
        for a = 1:length(shieldflag)
            neighborbranch = Rctbranches{neighborIdx(a)};
            neighborbranchLen = length(strfind(neighborbranch,'a')) + length(strfind(neighborbranch,'b'));
            if neighborbranchLen>addbranchLen
                shieldflag(a) = true;
            end
        end

        AllrxnList_steric(k) = sum(shieldflag);

    end
end

WaitMessage.Destroy;

% Construct transition probability matrix
TM = sparse(length(Glys));
for a = 1:length(AllrxnList_TMrow)
    TM(AllrxnList_TMrow(a),AllrxnList_TMcol(a)) = 1;
end
% add absoprtion
AbsFlag = find(contains(Glys,'[ab]'));
TM(AbsFlag,AbsFlag) = eye(length(AbsFlag));

% Normalize TM
TM = sparse(TM./(sum(TM,2)));

% linear indice for each reaction
AllrxnList_TMidx = sub2ind(size(TM),AllrxnList_TMrow,AllrxnList_TMcol);
AllrxnList_LacNAcLen_idx = AllrxnList_TMidx(AllrxnList_LacNAcLen~=0);
AllrxnList_LacNAcLen = AllrxnList_LacNAcLen(AllrxnList_LacNAcLen~=0);

% AllrxnList_RxnTypes & RxnTypes
AllrxnList_RxnTypes = AllrxnList(:,3);
for a = 1:length(AllrxnList_RxnTypes)
    if any(strcmp(branchingRxns,AllrxnList{a,3}))
        if AllrxnList_Br(a)~=0
            AllrxnList_RxnTypes{a} = [AllrxnList{a,3},'_B',num2str(AllrxnList_Br(a))];
        else
            AllrxnList_RxnTypes{a} = AllrxnList{a,3};
        end
    end
end

RxnTypes = unique(AllrxnList_RxnTypes);

end