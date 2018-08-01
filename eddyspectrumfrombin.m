function [] = eddyspectrumfrombin(fn,TE)
%bmatfrombin This function calculates the b-matrix of the binary external
%waveforms fed to the GE qti sequence (sampling is fixed at 4us, Gmax at
%70mt/m (diffusion for premier), and the refocusing duration of 10ms
close all
%fn = 'NOWbin2/TE60x0y0z20';
fid = fopen([fn '_x1']);
Gx1 = fread(fid,'int16','ieee-be');
fid = fopen([fn '_x2']);
Gx2 = fread(fid,'int16','ieee-be');

fid = fopen([fn '_y1']);
Gy1 = fread(fid,'int16','ieee-be');
fid = fopen([fn '_y2']);
Gy2 = fread(fid,'int16','ieee-be');

fid = fopen([fn '_z1']);
Gz1 = fread(fid,'int16','ieee-be');
fid = fopen([fn '_z2']);
Gz2 = fread(fid,'int16','ieee-be');

dt = 4e-6;
Gmax = 70;

Gx = [Gx1; zeros(round(10e-3/dt),1);Gx2];

Gy = [Gy1; zeros(round(10e-3/dt),1);Gy2];

Gz = [Gz1; zeros(round(10e-3/dt),1);Gz2];

scale = max([abs(Gx(:));abs(Gy(:));abs(Gz(:))]);
Gx = Gx*Gmax./scale;
Gy = Gy*Gmax./scale;
Gz = Gz*Gmax./scale;

G = [Gx,Gy,Gz];
q = 2*pi*42.58e3*cumsum(G)*dt;
%TE = 60;
N = length(G);

figure (1)
for i = 1:3
    for lambda = 1:250
        ec = exp((-[0:N*2-1]*dt)/lambda*1000);
        slew = diff(G)/dt;
        a = conv([slew(:,i); zeros(N,1)],ec)/(2*N);
        s(lambda) = a(round(TE/dt*1e-3));
    end
    if(max(s)<-min(s))
        s = -s;
    end
    hold all
    semilogy(s, 'linewidth', 3)
end

semilogy([1:250],0*[1:250],'--k')
legend('x','y','z')
xlabel('\lambda (ms)','fontsize',24)
ylabel('B_{EC}(TE) (a.u.)','fontsize',24)
set(gca,'fontsize',24)
set(gcf,'color','w')

figure (2)
t = [0:N-1]*dt*1e3;
plot(t, G, 'linewidth', 3);
xlabel('time (ms)','fontsize',24)
ylabel('G (mT/m)','fontsize',24)
%legend('Gx','Gy','Gz')
set(gca,'fontsize',24)
set(gcf,'color','w')
xlim([0,max(t)])
