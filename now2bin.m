function [] = now2bin(result,filename)
%now2bin This function converts the result structure from the NOW gradient
%design framework into binaries.  These binaries can be converted into epic
%readable waveforms using xlatebin and read directly using the
%createextwave macros


dt = result.optimizerProblem.dt*1e-3;
G = result.g;
N = length(G);
Gmax = result.optimizerProblem.gMax;
scale= 32766/Gmax;
zi = result.optimizerProblem.zeroGradientAtIndex;
% cut it into before and after refocusing
% %interpolate to scanner sampling rate
r = round(dt/4e-6);
g1 = G(1:zi(1)+1,:);
g1 = interp1([0:zi(1)],g1,[0:zi(1)*r]/r);


g2 = -G(zi(end)+1:end,:);
g2 = interp1([zi(end):N-1],g2,[zi(end):1/r:N-1]);

% waveforms should be even number of points
if(mod(length(g2),2)~=0)
    g2 = cat(1,g2,[0,0,0]);
end

if(mod(length(g1),2)~=0)
    g1 = cat(1,g1,[0,0,0]);
end

% create binary files...

wave2bin(g1(:,1),[filename '_x1'],scale);
wave2bin(g1(:,2),[filename '_y1'],scale);
wave2bin(g1(:,3),[filename '_z1'],scale);

wave2bin(g2(:,1),[filename '_x2'],scale);
wave2bin(g2(:,2),[filename '_y2'],scale);
wave2bin(g2(:,3),[filename '_z2'],scale);
end



function wave2bin(waveform,file,scale)

if (nargin < 2)
    error('Not enough arguments.');
end

if (nargin <3)
    scale= 32766/max(abs(waveform));
end

% -- Scale waveform
swave = 2*round(waveform * scale/2);

% -- Clip to max/min
f = find(swave > 32766);
swave(f) = 32766;
f = find(swave < -32768);
swave(f) = -32768;
%tt = sprintf('Max value = %g, scale = %g',max(abs(waveform)),scale);
%disp(tt);

swave(end)=swave(end)+1;	% Set EOS bit.
swave;
% write to binary
fid = fopen(file,'w');
for i = 1:length(swave)
    fwrite(fid, swave(i), 'int16','ieee-be');
end
fclose(fid);
end