function [ matched, constsum ] = matchesConstraints( xvalues, fileID )
    
    evalc(['[~, consts] = constraints ' fileID '(xvalues)']);

    constsum = sum(consts.*consts);
    matchcutoff = 1e-4;

    matched = (constsum < matchcutoff);

end
