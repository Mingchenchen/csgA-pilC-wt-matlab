function [ tracks ] = createMTracks2( t, K )

%[ tracks ] = createMTracks( tracks )
%
% Adds:
%   tracks.state which contains the state of the cell as calculated by
%         CellStateDetector()
%
%   tracks.density which contains the cell density from K at the position
%         of the cell
%
%   tracks.count which contins the number of state changes since the
%         beginning of the trajectory
%
%   NOTE: struct values in t are altered in output tracks!
%            ie t.x ~= tracks.x
%
% params
%   t: struct tracks
%      {
%         x
%         y
%         o
%         frame
%         id
%       }
%
% returns
%   struct tracks
%   {
%         x
%         y
%         v
%         o
%         id
%         frame
%         state
%         density
%         count
%         units = 'um'
%     }

BINS = size(K,1);

tracks.x = [];
tracks.y = [];
tracks.o = [];
tracks.v = [];
tracks.id = [];
tracks.frame = [];
tracks.state = [];
tracks.density = [];
tracks.count = [];
tracks.inten=[];
p = Progress(max(t.id));
for i = unique(t.id)'
    p.d(i)

    %Extract info on this cell
    fi = t.frame(t.id == i);
    xi = t.x(t.id == i);
    yi = t.y(t.id == i);
    oi = t.o(t.id == i);
    idi = t.id(t.id == i);
    ini=t.inten(t.id == i);
    %%%

    [~,S] = CellStateDetector(xi * 0.5136,yi * 0.5136 ,oi * -1);
    
    %Get the density`````   
   
     mX = round( xi(2:end) / 1920 * BINS);
     mY = round( yi(2:end) / 1440 * BINS);
     %sub2ind(size(K),mY,mX,fi(2:end))
     
     rho = K(sub2ind(size(K),mY,mX,fi(2:end)));
%      frame_start = fi(2:end);
%      rho = [];
%      for j = 1:length(mX)
%         rho = [rho; image.k(mY(j),mX(j),frame_start(j))];
%         %pX = round(pXY{j}(:,1) / 1920 * BINS);
%         %pY = round(pXY{j}(:,2) / 1440 * BINS);
%         %frame_mean = sub2ind(size(image.k),pY,pX,pF{j});
%         %image.rhoMean = [image.rhoMean; mean(image.k(frame_mean))];
%         %image.rhoDiff = [image.rhoDiff; (image.k(frame_mean(1)) - image.k(frame_mean(end)))];
%         %image.rhoAveDiff = [image.rhoAveDiff; mean(diff(image.k(frame_mean)))];
%      end
%     
    count = zeros(size(rho));
    cur_state = S(1,5);
    c = 1;
    for j = 1:size(count,1)
        if(S(j,5) ~= cur_state)
            c = c + 1;
            cur_state = S(j,5);
        end
        count(j) = c;
    end

    tracks.x = [tracks.x; S(:,1)];
    tracks.y = [tracks.y; S(:,2)];
    tracks.v = [tracks.v; S(:,3)];
    tracks.o = [tracks.o; S(:,4)];
    tracks.state = [tracks.state; S(:,5)];
    tracks.count = [tracks.count; count];
    tracks.density = [tracks.density; rho];
    tracks.frame = [tracks.frame; fi(2:end)];
    tracks.id = [tracks.id; idi(2:end)];
    tracks.inten = [tracks.inten; ini(2:end)];
    tracks.units = 'um';
end
p.done();

end

