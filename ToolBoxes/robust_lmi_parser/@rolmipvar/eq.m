function out = eq(A,B)

if ((isa(A,'rolmipvar')) && (isa(B,'rolmipvar')))
    if (length(A.data) ~= length(B.data))
        out = 0;
        return;
    else
        out = 1;
        for cont=1:length(A.data)
            if (A.data(cont).value ~= B.data(cont).value)
                out = 0;
                return
            end
        end
    end
elseif (isa(A,'rolmipvar'))
    out = 1;
    for cont=1:length(A.data)
        if (A.data(cont).value ~= B)
            out = 0;
            return
        end
    end
elseif (isa(B,'rolmipvar'))
    out = 1;
    for cont=1:length(B.data)
        if (B.data(cont).value ~= A)
            out = 0;
            return
        end
    end
end
return