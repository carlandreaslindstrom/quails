function [ beam ] = bun2matlab(homedir, filename)
    
    % read sdds beam file
    eledir = [homedir '/elegant/'];
    bindir = [eledir 'bin/'];
    tmpfile = [tempdir 'stream.tmp'];
    syscmd = [bindir 'sdds2stream ' eledir filename ' -columns=x,xp,y,yp,t,p' ' > ' tmpfile ];
    system(syscmd);
    beam = load(tmpfile);
    delete(tmpfile);
    
    % constants
    SI_c = 299792458;
    SI_me = 9.1093897e-31;
    SI_e = 1.60217733e-19;
    
    % format phase space to t->z[m] and p->E[GeV]
    beam(:,5) = (beam(:,5) - mean(beam(:,5))) * SI_c; % [m]
    beam(:,6) = beam(:,6) / (1e9*SI_e/(SI_me*SI_c^2)); % [GeV]

end
