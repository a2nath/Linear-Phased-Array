function garray = calcGain(dtheta, cow, bs_num, locations,mscords)

    % Redefine the constants because I don't want to pass as parameter
    lightspeed = 3e8;
    freq = 1900e6;
    dims = [20e-2 40e-2];
    lambda = lightspeed/freq;
    
%     /cow.location.(basestation{bs_num})
    basestations = fieldnames(cow.location);
    space = cow.antspace.(basestations{bs_num});
    alpha = cow.antalpha.(basestations{bs_num});
    angles = dtheta(2).(basestations{bs_num}) - cow.oriented.(basestations{bs_num}); % theta - thetaC
    N = cow.antcount.(basestations{bs_num});
    
    m = (pi*dims(1)*sin(angles))/lambda;
    singleant = (10*dims(1)*dims(2)/lambda.^2)*((sin(m)./m).^2).*((1 + cos(angles))/2).^2;
    i = find(isnan(singleant));
    singleant(i) = (10*dims(1)*dims(2)/lambda.^2).*((1 + cos(angles(i)))/2).^2;
    
    phee = (2*pi*space/lambda).*sin(angles) + alpha;
    garray = singleant.*(( sin( (N*phee)/2 )./ ( N*sin( phee/2 ) ) ).^2)*N;    
    k = find(isnan(garray));
    garray(k) = singleant(k).*N;    
    garray = 10*log10(garray);
    
%      --- debug ----
%     figure;scatter(real(locations),imag(locations),[],garray); hold on;
%     xlim([0 1000]);ylim([0 650]);
%     scatter(real(cow.location.(basestations{1})),imag(cow.location.(basestations{1})),[],'r');
%     scatter(real(cow.location.(basestations{2})),imag(cow.location.(basestations{2})),[],'r');
%     scatter(real(cow.location.(basestations{3})),imag(cow.location.(basestations{3})),[],'r');
%     scatter(real(mscords),imag(mscords),[],'r');  
% % % 
%     ii = find(garray > -25);
%     aa = angles(ii);
%     gg = garray(ii);
%     [~,jj] = sort(aa,'ascend');
%     figure;mmpolar(sort(aa(jj)),gg(jj))
% 
%     gainvec = zeros(1,15);
%     for i = 1:15
%         gainvec(i) = garray(find(locations == mscords(i)));
%     end
%     gainvec = gainvec';
%     display([(1:15)' gainvec]);
%     
end