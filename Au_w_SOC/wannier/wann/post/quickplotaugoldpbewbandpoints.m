clear all
%READ IN FOR WANNIER
numk=491;
numwann=16;
lenkwann=numk*numwann;
fileID = fopen('wannier90_geninterp.dat','r');
formatSpec='         %d   %f   %f   %f   %f\n';
fgets(fileID);
fgets(fileID);
fgets(fileID);
temp = fscanf(fileID,formatSpec,[5,lenkwann]);
fclose(fileID);
wannread=temp.';

%----------------------------------------------%
[energyd,total_dosd,efermid,pdosd] = import_doscar('DOSCARdos');
[eigenvalues,kpoints,nelect] = import_eigenval('EIGENVAL');
[geometry] = import_poscar('POSCAR');
ktotal = length(kpoints);
dx2=load('./banddata/dx2_up.dat');
dx2=dx2.'; % .' to transpose
dx2=max(dx2,1.e-15);
dxy=load('./banddata/dxy_up.dat');
dxy=dxy.'; 
dxy=max(dxy,1.e-15);
dxz=load('./banddata/dxz_up.dat');
dxz=dxz.';
dxz=max(dxz,1.e-15);
dyz=load('./banddata/dyz_up.dat');
dyz=dyz.';
dyz=max(dyz,1.e-15);
dz2=load('./banddata/dz2_up.dat');
dz2=dz2.';
dz2=max(dz2,1.e-15);
%Set to number of kpoints in the kwalk
%Set to 0 off and 1 on.
soc = 1;
%Set to 0 off and with soc and 1 on.
spin = 0;
%Set to 0 off and 1 on.
felect = 0;
%Read in POSCAR to populate number of species and atoms and get names
%Set numbers of atoms 
natoms = geometry.atomcount.';
%Number of Kpoints per line
%25 for small and 50 for large
%3/19/2024 This can be read in from kpoints line. Think about
%incorporating. Should be able to use vasplab to read in lines.
NKline=50;
%Lines duh.
%3/19/2024 This can be read in from kpoints line. Think about
%incorporating. Should be able to use vasplab to read in lines.
Lines=10;
kinterest = Lines*NKline;
%3/19/2024 For KPOINTS in line mode there is always an extra point that
%needs to be removed which is done in the next line.
z=0;
for i = kinterest-NKline:-1*NKline:NKline
    eigenvalues((ktotal-kinterest+1)+i,:,:)=[];
    dx2((ktotal-kinterest+1)+i,:)=[];
    wannread((ktotal-kinterest+1)+i,:,:)=[]
    dxy((ktotal-kinterest+1)+i,:)=[];
    dxz((ktotal-kinterest+1)+i,:)=[];
    dyz((ktotal-kinterest+1)+i,:)=[];
    dz2((ktotal-kinterest+1)+i,:)=[];
    z=z+1;
end

b=1;
kp=1;
j = 1;
check = 1;
for i = 1:size(wannread,1)
    kp = wannread(i,1);
    if kp == check;
        wannenergy(kp,b)=wannread(i,5);
        b=mod(b,numwann)+1;
    elseif kp ~= check;
        check = kp;
        b = 1;
    end
    % if mod(b,numwann)==0;
    %     kp=mod(kp,numk)+1;
    % end
end



%Max and min of plot
emax = 10;
emin = -10;
% z=0;
% for i = kinterest-NKline:-1*NKline:NKline
%     eigenvalues((ktotal-kinterest+1)+i,:,:)=[];
%     z=z+1
% end
%((1+spin+soc*3)*(felect*5+9))
%Sum out per atoms
sup=zeros(size(pdosd,1),size(natoms,2));
ktofatoms = 0;
for l = 1:size(natoms,2)
    for i = 1:natoms(1,l)
        for k = 1:(size(pdosd,2)-soc*(3/4)*size(pdosd,2))
            for j = 1:size(pdosd,1)
%Weirdly works out for both SOC and ISPIN=1                
                sup(j,l) = sup(j,l)+pdosd(j,k+3*soc*(k-1),i+ktofatoms);
            end
        end
    end
    ktofatoms = ktofatoms+natoms(1,l);
end
%Sum out per orbital and atom. Just gets total not be magnetization
super=zeros(size(pdosd,1),9,size(natoms,2));
ktofatoms = 0;
%1-numdiffspecies
for l = 1:size(natoms,2)
    %1-9 or 1-18
    for k = 0:(size(pdosd,2)-soc*(3/4)*size(pdosd,2))-1
        %1-number of atoms in species
        for i = 1:natoms(1,l)
            %1-energybins
            for j = 1:size(pdosd,1)
%Weirdly works out for both SOC and ISPIN=1          
                super(j,mod(k,9)+1,l) = super(j,mod(k,9)+1,l)+pdosd(j,mod(k,9)+1+3*soc*(mod(k,9)),i+ktofatoms);
            end
        end
    end
    ktofatoms = ktofatoms+natoms(1,l);
end
figure(1);
set(figure(1),'defaultLegendAutoUpdate','off')
subplot(1,2,1);
%Plot totaldos
if spin == 0
    plot(total_dosd,energyd-efermid,'b');
    %plot(total_dosd,energyd,'b');
end
if spin == 1
    plot(total_dosd(:,1)+total_dosd(:,2),energyd-efermid,'b');
end
%Plot element decomposed
hold on
% for i = 1:size(natoms,2)
%      plot(sup(:,i),energyd-efermid);
%      %plot(sup(:,i),energyd);
% end
%Plot element and orbital decomposed
%For s-p. Insert atom
% for i = 1:4
%      plot(super(:,i,1),energyd-efermid);
% end
%For d. Insert atom
for i = 5:9
     plot(super(:,i,1),energyd-efermid);
end
%plot(sup(:,3),energyd-efermid);
hold off
% legnames=[{'Total'},geometry.symbols];
% leg = legend(legnames)
% leg.Title.String = 'Atoms'
% for i = [1,3,5,7,9,11,13,15,17]
% for i = [1]
% for i = [3,5,7]
%     plot(sup(:,i),energyd-efermid);
% end
% leg = legend('total','s','py','pz','px','dxy','dyz','dzz','dxz','dxx-yy')
% leg = legend('total','s')
%leg = legend('total','s','py','pz','px')
leg = legend('total','d_{xy}','d_{yz}','d_{z^2}','d_{xz}','d_{x^2-y^2}');
%leg.Title.String = 'Bi s,p Orbitals';
leg.Title.String = 'Au d Orbitals';
%leg.Title.String = 'Sm s,p Orbitals';
%leg.Title.String = 'Sm d Orbitals';
%leg.Title.String = 'Zr s,p Orbitals';
%leg.Title.String = 'Zr d Orbitals';
leg.FontSize = 20;
leg.Location = 'northwest'
set(gca,'YLim',[emin emax],'XLim',[0,12],'LineWidth',2,'XDir','reverse','FontSize',20,'FontWeight','Bold')
%title('Checker Density of States')
title('DOS')
line([0 12],[0 0],'Marker','.','LineStyle','--','LineWidth',2, 'Color',[0 0 0])
%line([0 1],[10 10],'Marker','.','LineStyle','--','LineWidth',2, 'Color',[1 0 0])
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
%set(gca, 'XTickLabelMode', 'Manual')
set(gca, 'XTick', [2,4,6,8,10,12])
% mkdir dos
% tdos = horzcat(energyd-efermid,total_dosd(:,1)+total_dosd(:,2));
% save dos/tdos.txt tdos -ASCII
% odos = horzcat(energyd-efermid,sup(:,1));
% kdos = horzcat(energyd-efermid,sup(:,2));
% sdos = horzcat(energyd-efermid,sup(:,3));
% hdos = horzcat(energyd-efermid,sup(:,4));
% aldos = horzcat(energyd-efermid,sup(:,5));
% save dos/odos.txt odos -ASCII
% save dos/kdos.txt kdos -ASCII
% save dos/sdos.txt sdos -ASCII
% save dos/hdos.txt hdos -ASCII
% save dos/aldos.txt aldos -ASCII
%%%%%fix the kpoints spacing to be correct.%%%%%%%%%%%
%Hex walk Work on Getting these extracted from POSCAR
%Real space Lattice params
%0.0000000000	5.1617000000	12.4565500000
%2.9484500000	0.0000000000	12.4565500000
%2.9484500000	5.1617000000	0.0000000000
%got k lattice vecs
b1=[-1.5062936360	1.5062936360	1.5062936360]
b2=[1.5062936360	-1.5062936360	1.5062936360]
b3=[1.5062936360	1.5062936360	-1.5062936360]
highsymlat(1,:)=[0.0000 0.0000 0.0000]; %G
highsymlat(2,:)=[0.5000 0.0000 0.5000]; %X
highsymlat(3,:)=[0.5000 0.2500 0.7500]; %W
highsymlat(4,:)=[0.3750 0.3750 0.7500]; %K
highsymlat(5,:)=[0.5000 0.5000 0.5000]; %L
highsymlat(6,:)=[0.6250 0.2500 0.6250]; %U
for i = 1:size(highsymlat,1);
    highsymcart(i,:)= b1*highsymlat(i,1)+b2*highsymlat(i,2)+b3*highsymlat(i,3);
end
pathmag(1)=norm(highsymcart(1,:)-highsymcart(2,:)); %G->X
pathmag(2)=norm(highsymcart(2,:)-highsymcart(3,:)); %X->W
pathmag(3)=norm(highsymcart(3,:)-highsymcart(4,:)); %W->K
pathmag(4)=norm(highsymcart(4,:)-highsymcart(1,:)); %K->G
pathmag(5)=norm(highsymcart(1,:)-highsymcart(5,:)); %G->L
pathmag(6)=norm(highsymcart(5,:)-highsymcart(6,:)); %L->U
pathmag(7)=norm(highsymcart(6,:)-highsymcart(3,:)); %U->W
pathmag(8)=norm(highsymcart(3,:)-highsymcart(5,:)); %W->L
pathmag(9)=norm(highsymcart(5,:)-highsymcart(4,:)); %L->K
pathmag(10)=norm(highsymcart(6,:)-highsymcart(2,:)); %U->X
%This is for band structure made by line mode KPOINTS
%We define this one as the length scale and scale everything else by it
spac=pathmag(1);
points(1:NKline)=linspace(0,NKline-1,NKline);
k=NKline;
for i = 1:Lines-1
%points(x=start 1 after last:x+NKline-1)=linspace(points(last)+pathmag(i)/spac,points(last)+(NKline-1)*pathmag(i)/spac,NKline-1);
    points(k+1:k+NKline-1)=linspace(points(k)+pathmag(i+1)/spac,points(k)+(NKline-1)*pathmag(i+1)/spac,(NKline-1));
    k=k+NKline-1;
end
%%%%%fix the kpoints spacing to be correct.%%%%%%%%%%%
subplot(1,2,2);
if spin == 0
    %plot(points,eigenvalues((ktotal-kinterest+1):ktotal-z,:)-efermid,'b');
    %make black when using scatter ontop
    plot(points,eigenvalues((ktotal-kinterest+1):ktotal-z,:)-efermid,'k');
    %plot(points,eigenvalues((ktotal-kinterest+1):ktotal-z,:),'b');
    %plot(points,eigenvalues((ktotal-kinterest+1):ktotal-z,85:85+128),'b');
    hold on
    scatter(points,wannenergy-efermid,'r','filled');
    %eg
    %scatter(points,eigenvalues((ktotal-kinterest+1):ktotal-z,:)-efermid,dx2((ktotal-kinterest+1):ktotal-z,:)*50,[0.30,0.75,0.93],'filled');
    %scatter(points,eigenvalues((ktotal-kinterest+1):ktotal-z,:)-efermid,dz2((ktotal-kinterest+1):ktotal-z,:)*50,[0.49,0.18,0.56],'filled');
    %t2g
    %scatter(points,eigenvalues((ktotal-kinterest+1):ktotal-z,:)-efermid,dxy((ktotal-kinterest+1):ktotal-z,:)*50,[0.85,0.33,0.10],'filled');
    %scatter(points,eigenvalues((ktotal-kinterest+1):ktotal-z,:)-efermid,dxz((ktotal-kinterest+1):ktotal-z,:)*50,[0.47,0.67,0.19],'filled');
    %scatter(points,eigenvalues((ktotal-kinterest+1):ktotal-z,:)-efermid,dyz((ktotal-kinterest+1):ktotal-z,:)*50,[0.93,0.69,0.13],'filled');
    hold off
end
if spin == 1
    plot(points,eigenvalues((ktotal-kinterest+1):ktotal-z,:,1)-efermid,'b',points,eigenvalues((ktotal-kinterest+1):ktotal-z,:,2)-efermid,'r');
end
set(gca,'YLim',[emin emax],'XLim',[0,points(end)],'LineWidth',1.7,'FontSize',20,'FontWeight','Bold','YAxisLocation','right')
for l = 1:(Lines-1)
    line([points(NKline+((NKline-1))*(l-1)) points(NKline+((NKline-1))*(l-1))],[emin emax],'Marker','.','LineStyle','--','LineWidth',1, 'Color',[0 0 0])
end    
set(gca, 'XTickLabelMode', 'Manual')
set(gca, 'XTick', [])
%sxlabel={'$\Gamma$' 'M' '$\Gamma^{\prime}$' '$\Gamma$'};
%High-sym points on X-axis. May have to adjust position.
H2(1) = 0;
V2(1) = emin/20.0-0.1+emin;
for l = 1:(Lines-1)
    H2(l+1) = points(NKline+(NKline-1)*(l-1));
    V2(l+1) = emin/20.0-0.1+emin;
end
H2(Lines+1) = points(end);
V2(Lines+1) = emin/20.0-0.1+emin;
line([0 points(end)],[0 0],'Marker','.','LineStyle','--','LineWidth',2, 'Color',[0 0 0])
%line([0 points(end)],[10 10],'Marker','.','LineStyle','--','LineWidth',2, 'Color',[1 0 0])
sxlabel={'$\Gamma$', 'X', "W", "K", '$\Gamma$', 'L', 'U', 'W','L','K$|$U','X'}
%sxlabel={'$\Gamma$', 'M', "K", '$\Gamma$', 'A', 'L', 'H','A$|$L','M$|$K','H'}
%sxlabel={'$\Gamma$', 'M', "K", '$\Gamma$'}
%sxlabel={'F$_1$', '$\Gamma$', 'L', 'T$_1|$T', '$\Gamma$', 'F$_2$'};
%sxlabel={'$\Gamma$', 'K', 'M', '$\Gamma$'}
%sxlabel={'$\Gamma$', 'M$_a$', 'Y', 'M$_b$', '$\Gamma$', 'Y$|$M$_a$', 'M$_b$', 'M$_b$`', '$\Gamma$`', 'M$_a$`', 'Y`', 'M$_b$`', 'M$_a$`$|\Gamma$`', 'Y`', 'Y$|\Gamma$', '$\Gamma$`$|$M$_a$', 'M$_a$`'}
%Make title
%title('NdTi3Bi4 frozenmin=2 frozenmax=7')
%title('Au PBE w/SOC Primitive Cell a=2.94955Ang \theta =60deg ENCUT=600eV 15x15x15KPs')
title('Au PBE w/SOC Primitive Cell a=2.94955Ang \theta =60deg ENCUT=600eV 15x15x15KPs')
%title('AFM CHECKERBOARD SCAN + rvv10 Bands')
%title('AFM CHECKERBOARD Different POSCAR SCAN + rvv10 Bands')
%title('DOS GARBO AFM CHECKERBOARD ENCUT 1200eV SCAN + rvv10 Bands')
ylabel('E-E_F (eV)')
figx = text(H2,V2,sxlabel,'HorizontalAlignment','center','Interpreter', 'latex');
set(figx,'FontSize',20,'FontWeight','Bold')
set(findall(gca, 'Type', 'Line'),'LineWidth',1.75);
set(subplot(1,2,1),'position',[.04, .09, .17, .82])
set(subplot(1,2,2),'position',[.22, .09, .7, .82])
set(gca, 'FontWeight','Bold','fontsize',20,'linewidth',2)
% mkdir spinupbands
% for i = 1:size(eigenvalues,2)
%     bands = horzcat(points.',eigenvalues((ktotal-kinterest+1):ktotal-z,i,1)-efermid);
%     name="spinupbands/spinupband"+num2str(i)+".txt";
%     save(name, 'bands', '-ascii')
% end
% mkdir spindownbands
% for i = 1:size(eigenvalues,2)
%     bands = horzcat(points.',eigenvalues((ktotal-kinterest+1):ktotal-z,i,2)-efermid);
%     name="spindownbands/spindownband"+num2str(i)+".txt";
%     save(name, 'bands', '-ascii')
% end
