function [ output ] = trackLattice(homedir)

    % run elegant file
    filename = 'params.ele';
    eledir = [homedir '/elegant/'];
    bindir = [eledir 'bin/'];
    syscmd = [bindir 'elegant ' eledir filename ' -rpnDefns=' bindir 'defns.rpn'];
    system(syscmd);

    % collect and return
    addpath('bin');
    initialBeam = bun2matlab(homedir, 'initial.bun');
    finalBeam = bun2matlab(homedir, 'final.bun');
    [enx0, eny0] = beamStatistics(initialBeam);
    [enx, eny, betx, bety, alfx, alfy] = beamStatistics(finalBeam);
    output = [enx, eny, enx0, eny0, betx, bety, alfx, alfy];

end
