function [ output ] = trackParticle(homedir)

    % run elegant file
    filename = 'params.ele';
    eledir = [homedir '/elegant/'];
    bindir = [eledir 'bin/'];
    syscmd = [bindir 'elegant ' eledir filename ' -rpnDefns=' bindir 'defns.rpn'];
    system(syscmd);

    % collect and
    addpath('bin');
    
    initialBeam = bun2matlab(homedir, 'initial.bun');
    finalBeam = bun2matlab(homedir, 'final.bun');
    output = [finalBeam(1,1:4), initialBeam(1,1:4)];
    

end
