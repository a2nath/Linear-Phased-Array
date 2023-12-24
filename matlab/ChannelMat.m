function hmatrix = ChannelMat(dtheta,cow,bs_num,locations, mscords)
   
    lightspeed = 3e8;
    freq = 1900e6;
    lambda = lightspeed/freq;
    Grx = -2; % dBi
    
    basestations = fieldnames(cow.location);

    Gtx = cow.antgains.(basestations{bs_num});
    dist = dtheta(1).(basestations{bs_num});
    pathloss = 20*log10((4*pi)/lambda) + 20*log10(dist);
    hmatrix = Gtx + Grx - pathloss;
    
    % --- debug ----
%     figure;scatter(real(locations),imag(locations),[],hmatrix); hold on;
%     xlim([1 1000]);ylim([1 650]);
%     scatter(real(cow.location.(basestations{1})),imag(cow.location.(basestations{1})),[],'r');
%     scatter(real(cow.location.(basestations{2})),imag(cow.location.(basestations{2})),[],'r');
%     scatter(real(cow.location.(basestations{3})),imag(cow.location.(basestations{3})),[],'r');
%     scatter(real(mscords),imag(mscords),[],'r');  
% 
%     ii = find(hmatrix > -20);
%     angles = dtheta(2).(basestations{bs_num}) - cow.oriented.(basestations{bs_num}); % theta - thetaC
%     angles = angles(ii);
%     gg = hmatrix(ii);
%     [~,jj] = sort(angles,'ascend');
%     figure;mmpolar(sort(angles(jj)),gg(jj))
%     
%     gainvec = zeros(1,15);
%     for i = 1:15
%         gainvec(i) = hmatrix(find(locations == mscords(i)));
%     end
%     gainvec = gainvec';
%     display(gainvec);
%     

% UNCOMMENT !!!!!!!!!!!!
     hmatrix = 10.^(hmatrix./10)';
end