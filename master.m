clc; clear;
close all;
%% Define parameters
freq = 1900e6;


% DEBUG = 1;
DEBUG = 0;
SymbRate = 3.84e6;
BlocksPerSymb = 768;
lightspeed = 3e8;
lambda = lightspeed/freq;
TSLOTS = 5;

NORTH = 90;
NORTHEAST = 75; %trial and error values
NORTHWEST = 180 - NORTHEAST;
CENTERBS = 500;
OFF = 250; %trial and error values
HEIGHT = 200;

Count = [5 3 5];
Space = [ ...
    .300 ...
    .4 ...
    .300 ...
    ];
if numel(find(Space <= .2) ~= 0) || ...
    numel(find(Count*.2 + (Count - 1).*Space >= 2.4) ~= 0)
     ME = MException('MYFUN:incorrectSize','Double-check spacing constraints');
     throw(ME);
end

alphatable = [  10, 65, 35
                -10, -65, -35;
              -75 -40 75;
              40, 65, -40;
              40, -65, -40];
cow.antcount = struct('b1',Count(1),'b2',Count(2),'b3',Count(3));
cow.antspace = struct('b1',Space(1),'b2',Space(2),'b3',Space(3)); 
cow.antpower = struct('b1',0,'b2',0,'b3',0); % default 0 dbm.
cow.oriented = struct('b1',deg2rad(NORTHEAST),'b2',deg2rad(NORTH),'b3',deg2rad(NORTHWEST)); % default 90 degrees.
cow.location = struct('b1',CENTERBS-OFF+HEIGHT*i,'b2',CENTERBS+180i,'b3',CENTERBS+OFF+HEIGHT*i);

display(['Timeslot#     ' 'MS Served      ' 'Pxi chosen     ' 'alphai chosen      ' 'SINR at each MS']);
for t = 1:TSLOTS
%     for angle1 = -90:2:90
%         for angle2 = -90:2:90
%             for angle3 = -90:2:90
%                 alphatable(t,:) = [angle1, angle2, angle3];
                angles = alphatable(t,:);
%             a_range = -90:2:90;
%             gaintable = zeros(15,numel(a_range));
%             for aa = 1:numel(a_range);

%                 angles = [-a_range(aa) -80, a_range(aa)];
                cow.antalpha = struct('b1',deg2rad(angles(1)),'b2',deg2rad(angles(2)),'b3',deg2rad(angles(3))); % default 0 degrees.
                cow.antgains = struct('b1',[],'b2',[],'b3',[]);
                cow.chmatrix = struct('b1',[],'b2',[],'b3',[]);


                u = [250:100:750 250:100:750];
                v = [500 500 500 500 500 500 600 600 600 600 600 600];
                u2 = 450:50:550;
                v2 = [325 325 325];

                mscords = zeros(1,length(u) + length(u2));
                idx = 1;
                for m = 1:length(u)
                    mscords(m) = u(m) + 1j*v(m);
                end
                l = length(u);
                for m = l+1:length(u2)+l
                    mscords(m) = u2(m-l) + 1j*v2(m-l);
                end
                basestations = fieldnames(cow.location);
            %     scatter(real(mscords),imag(mscords))
            %     hold on
            %     scatter(real(cow.location.(basestations{1})),imag(cow.location.(basestations{1})),[],'r')
            %     scatter(real(cow.location.(basestations{2})),imag(cow.location.(basestations{2})),[],'r')
            %     scatter(real(cow.location.(basestations{3})),imag(cow.location.(basestations{3})),[],'r')
            %     xlim([0 1000])
            %     ylim([0 650])
                %hold off
                if (DEBUG == 0) 
                    overpoints = mscords;
                else
                    overpoints = zeros(1,length(1:10:1000) *length(1:6.5:650));
                    idx = 1;
                    for m = 1:5:1000
                        for n = 1:3.25:650
                            overpoints(idx) = m + 1j*n;
                            idx = idx + 1;
                        end
                    end
                    overpoints = [overpoints mscords];
                end


                l = overpoints - cow.location.(basestations{1});
                [th1, a] = cart2pol(real(l),imag(l));
                % distances from BS2
                l = overpoints - cow.location.(basestations{2});
                [th2, b] = cart2pol(real(l),imag(l));
                % distances from BS3
                l = overpoints -cow.location.(basestations{3});
                [th3, c] = cart2pol(real(l),imag(l));

                dtheta = struct('b1', {a th1}, 'b2', {b th2}, 'b3', {c th3});



%               gaintable(:,aa) = calcGain(dtheta,cow,2,overpoints,mscords);
% 
% 
% 
%             end

                cow.antgains.(basestations{1}) = calcGain(dtheta,cow,1,overpoints,mscords);
                cow.antgains.(basestations{2}) = calcGain(dtheta,cow,2,overpoints,mscords);
                cow.antgains.(basestations{3}) = calcGain(dtheta,cow,3,overpoints,mscords);

                cow.chmatrix.(basestations{1}) = ChannelMat(dtheta,cow,1,overpoints,mscords);
                cow.chmatrix.(basestations{2}) = ChannelMat(dtheta,cow,2,overpoints,mscords);
                cow.chmatrix.(basestations{3}) = ChannelMat(dtheta,cow,3,overpoints,mscords);

                [cow, msserved, SINR] = gPower(dtheta,cow,t,overpoints,mscords);

                display([num(t) '             '...
                    num(msserved) '       '...
                    num([cow.antpower.(basestations{1}) cow.antpower.(basestations{2})...
                        cow.antpower.(basestations{3})]) '      '...
                    num(angles) '         '...
                    num(SINR)]);
%                 if SINR(1) > 24 && SINR(2) > 24 && SINR(3) > 24
%                     display('done');
%                     break;
%                 end
    
%             end
%         end
%     end
%     
   end