function [cow, ms, SINR] = gPower(dtheta,cow,iteration_num,locations, mscords)
    lightspeed = 3e8;
    freq = 1900e6;
    lambda = lightspeed/freq;
    Grx = -2; % dBi
    
    bw = 20e6;
    NF = 5;
    THERMAL_DB = -174 + round(10*log10(bw)) + NF;
    THER = 10^(THERMAL_DB/10);
    
    basestations = fieldnames(cow.location);
    cmatrix = [cow.chmatrix.(basestations{1}) cow.chmatrix.(basestations{2}) ...
            cow.chmatrix.(basestations{3})];
        
%     timeslot.powers = struct('b1',[],'b2',[],'b3',[]);
% stations = [1 14 6];
     stations = [2 13 11;
                8 15 5;
                1 14 6;
                7 4 10;
                9 3 12];

%         find first max bs and ms
    ms = stations(iteration_num,:);   
   
    breakflag = 0;
    powers = zeros(1,3);
    SINR = zeros(1,3);
    table = cmatrix(ms,:);
    idx = 1;
    
    cutoffdB = 24;
    cutofflinear = 10^(cutoffdB/10);
    maxvec = zeros(6,30*30*30);
    
    
    
    for power1 = -30:30
        p1 = 10^((power1-30)/10);
        for power2 = -30:30
            p2 = 10^((power2-30)/10);
            for power3 = -30:30
                p3 = 10^((power3-30)/10); 
                
                ptable = [p1*table(1,1) / addrms([p2*table(1,2) p3*table(1,3) THER])
                    p2*table(1,2) / addrms([p1*table(1,1) p3*table(1,3) THER])
                    p3*table(1,3) / addrms([p2*table(1,2) p1*table(1,1) THER])
                    p1*table(2,1) / addrms([p2*table(2,2) p3*table(2,3) THER])
                    p2*table(2,2) / addrms([p1*table(2,1) p3*table(2,3) THER])
                    p3*table(2,3) / addrms([p2*table(2,2) p1*table(2,1) THER])
                    p1*table(3,1) / addrms([p2*table(3,2) p3*table(3,3) THER])
                    p2*table(3,2) / addrms([p1*table(3,1) p3*table(3,3) THER])
                    p3*table(3,3) / addrms([p2*table(3,2) p1*table(3,1) THER])];
                
                maxvec(1,idx) = power1;
                maxvec(2,idx) = power2;
                maxvec(3,idx) = power3;
                maxvec(4,idx) = max(ptable(1:3));
                maxvec(5,idx) = max(ptable(4:6));
                maxvec(6,idx) = max(ptable(7:9));
                idx = idx + 1;
                
%                 if ( ptable(1) >= cutofflinear || ptable(2) >= cutofflinear || ptable(3) >= cutofflinear )
%                     display('hit')
%                 else
%                     
%                     continue;
%                 end
%                     
% 
%                 if ( ptable(4) >= cutofflinear || ptable(5) >= cutofflinear || ptable(6) >= cutofflinear )
%                     display('hit')
%                 else
%                     
%                     continue;
%                 end
% 
%                 if ptable(7) >= cutofflinear || ptable(8) >= cutofflinear || ptable(9) >= cutofflinear 
% %                     display(['all hit p1:' num2str(power1) ' p2:' num2str(power2) ' p3:' num2str(power3)])
% %                     maxvec(idx) = max(ptable(7:9));                    
%                     breakflag = 1;
%                     powers = [power1 power2 power3];
                  
%                     break;
%                 end
                 
            end    
%             if(breakflag == 1)
%                 break;
%             end
        end
%         if(breakflag == 1)
%             break;
%         end
    end
%     maxvec(4:6,:) = 10.*log10(maxvec(4:6,:));
%     find(
    [m,i] = max(maxvec(4:6,:),[],2);
    SINR = (10.*log10(m))';
%     display(m);
    cow.antpower.(basestations{1}) = maxvec(1,i(1));
    cow.antpower.(basestations{2}) = maxvec(2,i(2));
    cow.antpower.(basestations{3}) = maxvec(3,i(3));
end