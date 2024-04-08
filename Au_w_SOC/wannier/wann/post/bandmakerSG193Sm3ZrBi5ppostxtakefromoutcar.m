%KPOINTS were taken from pbe bands calculation from the outcar of the
%calculation. This is specfic to the system.
clear all
format long
lines = 9;
kpoints = 50;
kpvals=load('kpointsweneed');
kpvals=kpvals.';
kgrid = kpvals(1:3,:)
ids = zeros(1,lines*(kpoints));
for i = 1:lines*(kpoints)
    ids(1,i)=i;
end    
kgrid = vertcat(ids,kgrid);
kgridtran=kgrid.';
fileID = fopen('sg193kwalkpostwtakefromoutcar.wannier90_geninterp.kpt','w');
formatSpec = '%d %8.8f %8.8f %8.8f\n';
%for i = 1:lines*(kpoints)+1
fprintf(fileID,formatSpec,kgrid)
%end
fclose(fileID);

    
