function [cow, ms, SINR] = gPower(dtheta,cow,iteration_num,locations, mscords)
    bs_count = cow.size;
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
    powers = zeros(1, bs_count);
    SINR   = zeros(1, bs_count);
    table  = cmatrix(ms,:);
    idx = 1;
    
    cutoffdB = 24;
    cutoff = 10^(cutoffdB/10);
    maxvec = zeros(6,30*30*30);
    
    sinr_max     = cutoff * bs_count;
    sinr_list    = zeros(1, bs_count);
    max_tx_power = zeros(1, bs_count);
    %base_id_num  = zeros(1, bs_count);
    
    for power1 = -30:30
        p1 = 10^((power1-30)/10);
        for power2 = -30:30
            p2 = 10^((power2-30)/10);
            for power3 = -30:30
                p3 = 10^((power3-30)/10); 
                
                power_table = [p1 p2 p3];
                rx_power_possibilities = zeros(6, bs_count);

                rx_power_possibilities(1,:) = [                                                                      ...
                    power_table(1)*table(1,1) / addrms([power_table(2)*table(1,2) power_table(3)*table(1,3) THER]),  ...
                    power_table(2)*table(2,2) / addrms([power_table(1)*table(2,1) power_table(3)*table(2,3) THER]),  ...
                    power_table(3)*table(3,3) / addrms([power_table(2)*table(3,2) power_table(1)*table(3,1) THER])   ...
                ];

                rx_power_possibilities(2,:) = [                                                                      ...
                    power_table(1)*table(1,1) / addrms([power_table(2)*table(1,2) power_table(3)*table(1,3) THER])   ...
                    power_table(3)*table(2,3) / addrms([power_table(2)*table(2,2) power_table(1)*table(2,1) THER]),  ...
                    power_table(2)*table(3,2) / addrms([power_table(1)*table(3,1) power_table(3)*table(3,3) THER])   ...
                ];

                rx_power_possibilities(3,:) = [                                                                      ...
                    power_table(2)*table(1,2) / addrms([power_table(1)*table(1,1) power_table(3)*table(1,3) THER]),  ...
                    power_table(1)*table(2,1) / addrms([power_table(2)*table(2,2) power_table(3)*table(2,3) THER]),  ...
                    power_table(3)*table(3,3) / addrms([power_table(2)*table(3,2) power_table(1)*table(3,1) THER])   ...
                ];

                rx_power_possibilities(4,:) = [                                                                      ...
                    power_table(1)*table(3,1) / addrms([power_table(2)*table(3,2) power_table(3)*table(3,3) THER]),  ...
                    power_table(2)*table(2,2) / addrms([power_table(1)*table(2,1) power_table(3)*table(2,3) THER]),  ...
                    power_table(3)*table(1,3) / addrms([power_table(2)*table(1,2) power_table(1)*table(1,1) THER])   ...
                ];

                rx_power_possibilities(5,:) = [                                                                      ...
                    power_table(3)*table(1,3) / addrms([power_table(2)*table(1,2) power_table(1)*table(1,1) THER]),  ...
                    power_table(1)*table(2,1) / addrms([power_table(2)*table(2,2) power_table(3)*table(2,3) THER]),  ...
                    power_table(2)*table(3,2) / addrms([power_table(1)*table(3,1) power_table(3)*table(3,3) THER])   ...
                ];

                rx_power_possibilities(6,:) = [                                                                      ...
                    power_table(2)*table(1,2) / addrms([power_table(1)*table(1,1) power_table(3)*table(1,3) THER]),  ...
                    power_table(3)*table(2,3) / addrms([power_table(2)*table(2,2) power_table(1)*table(2,1) THER]),  ...
                    power_table(1)*table(3,1) / addrms([power_table(2)*table(3,2) power_table(3)*table(3,3) THER])   ...
                ];
               
              
				/* this change reveals that the current setup and timeslot does not find a combo which will hit the cuttof */
                for row = 1:length(rx_power_possibilities)
                    if (cutoff >= rx_power_possibilities(row,1) && cutoff >= rx_power_possibilities(row,2) && cutoff >= rx_power_possibilities(row,3) )
                        thesum = sum(rx_power_possibilities(row,:));
                        if (thesum > sinr_max)
                            sinr_max     = thesum;
                            sinr_list    = rx_power_possibilities(row,:); % stage 2: don't get the maximum, get the priority queue version of this
                            max_tx_power = power_table;
                            %base_id_num  = row;
                        end
                    end
                end
                
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
    %if (length(unique(base_id_num)) ~= length(base_id_num))
    %    error("assignment is not unique")
    %end

    SINR = 10.*log10(sinr_list);
%     display(m);
    cow.antpower.(basestations{1}) = 10*log10(max_tx_power(1)) + 30;
    cow.antpower.(basestations{2}) = 10*log10(max_tx_power(2)) + 30;
    cow.antpower.(basestations{3}) = 10*log10(max_tx_power(3)) + 30;
end
