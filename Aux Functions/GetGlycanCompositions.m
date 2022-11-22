function compositions = GetGlycanCompositions(fileName)

if ischar(fileName)
    [~,linearcodes] = xlsread(fileName,'Annotation','B:B');
    alltext = xlsread(fileName,'Annotation');
    linearcodes = linearcodes(2:end);
    annotmat = alltext(:,3:end);
else
    linearcodes = fileName;
end

compositions = cell(size(linearcodes));

for k1 = 1:length(linearcodes)
    str = linearcodes{k1};
    composition_monoCounts_pre(k1,1) = length(strfind(str,'NN'));
    composition_monoCounts_pre(k1,2) = length(strfind(str,'F'));
    composition_monoCounts_pre(k1,3) = length(strfind(str,'A')) + length(strfind(str,'M'))...
        -length(strfind(str,'Asn'))-length(strfind(str,'AN'));% discount for Asn and AN
    composition_monoCounts_pre(k1,4) = length(strfind(str,'GN'));
end

generaltypes = {'NeuAc','dHex','Hex','HexNAc'};

for a = 1:length(compositions)
    str = [];
    for b = 1:length(generaltypes)
        if composition_monoCounts_pre(a,b)~=0
            if composition_monoCounts_pre(a,b) == 1
                str = [str,generaltypes{b}];
            else
                str = [str,generaltypes{b},num2str(composition_monoCounts_pre(a,b))];
            end
        end
    end
    compositions{a} = str;
end

if ischar(fileName)
    compositions_unique = unique(compositions);
    [~,idx] = ismember(compositions,compositions_unique);

    % Write
    unique_Len = length(compositions_unique);
    Len = length(idx);
    z = 'A':'Z';
    signalmat = zeros(unique_Len,size(annotmat,2));
    for a = 1:size(signalmat,1)
        for b = 1:size(signalmat,2)
            signalmat(a,b) = sum(annotmat(a == idx,b));
        end
    end

    out = arrayfun(@(x)z(rem(floor(x*26.^(1-floor(log(x)/log(26)+1):0)),26)),size(signalmat,2)+2,'un',0);
    annotmat = double(annotmat~=0);

    xlswrite(fileName,idx,'Annotation',['A2:A',num2str(Len+1)]);
    xlswrite(fileName,(1:unique_Len)','MS Raw',['A2:A',num2str(unique_Len+1)]);
    xlswrite(fileName,compositions_unique,'MS Raw',['B2:B',num2str(unique_Len+1)]);
    xlswrite(fileName,signalmat,'MS Raw',['C2:',out{1},num2str(unique_Len+1)]);
    xlswrite(fileName,annotmat,'Annotation',['C2:',out{1},num2str(unique_Len+1)]);
end

end