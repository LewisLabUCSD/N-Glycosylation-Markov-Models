function compositions = GetGlycanCompositions(linearcodes)

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

end