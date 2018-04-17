function L = combustorDesign(UInf, altitude, engine, fluid)
%This outputs the zonal lengths of a nominal RQL combustor.
    % TODO : compute the length of the evaporation and mixing zone

    % TODO : compute the length of the dilution zone

    % TODO : compute the length of the secondary reaction zone

    % total length
    L.LTotal = L.LEvap + L.LDil + L.LSrz;
end

