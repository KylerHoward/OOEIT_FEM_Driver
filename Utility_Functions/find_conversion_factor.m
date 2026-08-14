function conversion_factor = find_conversion_factor(unit, operation)
    %{
    Determine what conversion factor is needed from metadata units to 
    standard units based on multiplication or division
    8/13/26 - Kyler Howard

    param: unit      - Character array of the input unit
    param: operation - String/Char of multiplication ("*") or division ("/")

    return: conversion_factor - How much to multiply/divide the value by
    %}

    if operation ~= "*" && operation ~= "/"
        error("Did not give multipication ('*') or division ('/') operator")
    end

    unit = char(unit);
    if isscalar(unit)
        conversion_factor = 1;
    elseif strcmp(unit, 'in')
        conversion_factor = 0.0254;
    elseif unit(1) == 'c'
        conversion_factor = 1e-2;
    elseif unit(1) == 'm'
        conversion_factor = 1e-3;
    elseif unit(1) == 'µ' || unit(1) == 'u'
        conversion_factor = 1e-6;
    elseif unit(1) == 'n'
        conversion_factor = 1e-9;
    end

    if operation == "/"
        conversion_factor = 1 / conversion_factor;
    end
end