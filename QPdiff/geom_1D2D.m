%let's visually show how 2D geometries improve with longer diffuion
%lengths

ltes=1;
lfin=[0.25,4]';

for jfin=1:length(lfin)
    
    figure(jfin)
    set(jfin,'position',[50 50 525 475]);
    %---------- 1D QP Propagation ----------
    %let's offset from the origin by 1.1 fin lengths
    dx= -1.1*lfin(jfin)
    dy= 0;
    
    %TES
    xtes = [0,0]';
    ytes = ltes*[-0.5,0.5]';
    
    plot(dx+xtes,dy+ytes,'-m')
    hold on
    %FIN
    yfin = ltes*[-0.5, 0.5, 0.5, -0.5, -0.5]';
    xfin = lfin(jfin)*[-1, -1, 1, 1, -1]';
    plot(dx+xfin,dy+yfin,'-b')
    
    %---------- 2D QP Propagation ----------
    %let's offset from the origin by 1.1 fin lengths
    dx= 1.1*lfin(jfin)
    dy= 0;
    
    %TES
    xtes = [0,0]';
    ytes = ltes*[-0.5,0.5]';   
    plot(dx+xtes,dy+ytes,'-m')
 
    %FIN
    
    % rectangular section
    yfin = ltes*[-0.5, 0.5, 0.5, -0.5, -0.5]';
    xfin = lfin(jfin)*[-1, -1, 1, 1, -1]';
    plot(dx+xfin,dy+yfin,'-b')
    
    % 2D section:
    phid = [0:180]';
    xfin = 0 + lfin(jfin)*cosd(phid);
    yfin = ltes/2+lfin(jfin)*sind(phid);
    plot(dx+xfin,dy+yfin,'-b');
    
    phid = [180:360]';
    xfin = 0 + lfin(jfin)*cosd(phid);
    yfin = -ltes/2+lfin(jfin)*sind(phid);
    plot(dx+xfin,dy+yfin,'-b');
    
    hold off
    
    axis square
    dlim = max(1.1*(ltes/2+lfin(jfin)), 1.1*3.1*lfin(jfin));
    xlim([-dlim,dlim])
    ylim([-dlim,dlim])
    
    xlabel('length [1/tes]')
    ylabel('length [1/tes]')
    title(['1D and 2D  QET Geometries for various fin/tes ratio'])
    
end    
    