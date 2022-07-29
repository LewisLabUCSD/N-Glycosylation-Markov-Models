function  [TM, AllrxnList_TMidx, AllrxnList_Br, AllrxnList_RxnTypes, RxnTypes,AllrxnList_steric,stericRxns] = ConstructTPM(AllrxnList,Glys)

fprintf('Construct Transition Probability Matrix from the generic N-glycosylation network...\n');

N = size(AllrxnList,1);
WaitMessage = parfor_wait(N, 'Waitbar', true);
AllrxnList_TMrow = zeros(size(AllrxnList,1),1);
AllrxnList_TMcol = AllrxnList_TMrow;
AllrxnList_Br = AllrxnList_TMrow;
AllrxnList_steric = AllrxnList_TMrow;

% List reactions to consider steric interactions
stericRxns = {'a3SiaT','b4GalT','iGnT','GnTIV','GnTV'};

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
    
    % Distinguish branches and identify steric hindrance for
    % 'a3SiaT','b4GalT', 'iGnT','GnTIV', & 'GnTV')
    if ismember(RxnType,stericRxns)
        [RctFlag,~,Rctbranches] = IdentifyEndMoiety(Rct);
        [~,~,Probranches] = IdentifyEndMoiety(Pro);
        loc = find(~strcmp(Rctbranches,Probranches));
        
        AllrxnList_Br(k) = loc;  
        AllrxnList_steric(k) = sum(RctFlag(1:4~=loc));
    end   
end

WaitMessage.Destroy;

% Construct transition probability matrix
TM = zeros(length(Glys));
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

% AllrxnList_RxnTypes & RxnTypes
AllrxnList_RxnTypes = cell(size(AllrxnList_TMidx));
for a = 1:length(AllrxnList_RxnTypes)
    if AllrxnList_Br(a)~=0
        AllrxnList_RxnTypes{a} = [AllrxnList{a,3},'_B',num2str(AllrxnList_Br(a))];
    else
        AllrxnList_RxnTypes{a} = AllrxnList{a,3};
    end
end

RxnTypes = unique(AllrxnList_RxnTypes);

end