% Copyright 2010
% Massachusetts Institute of Technology
% All rights reserved
% Developed by the Aerospace Controls Lab, MIT

%---------------------------------------------------------------------%
% Runs consensus between neighbors
% Checks for conflicts and resolves among agents
% 
% This is a message passing scheme described in Table 1 of:
% "Consensus-Based Decentralized Auctions for Robust Task Allocation", 
% H.-L. Choi, L. Brunet, and J. P. How, 
% IEEE Transactions on Robotics, Vol. 25, (4): 912 � 926, August 2009
%---------------------------------------------------------------------%

function [CBBA_Data t] = CBBA_Communicate(CBBA_Params, CBBA_Data, Graph, old_t, T)

% Copy data
for n = 1:CBBA_Params.N,
    old_z(n,:) = CBBA_Data(n).winners;
    old_y(n,:) = CBBA_Data(n).winnerBids;    
end

z = old_z;
y = old_y;
t = old_t;

epsilon = 10e-6;


% Start communication between agents

% sender   = k
% receiver = i
% task     = j

for k=1:CBBA_Params.N
    for i=1:CBBA_Params.N
        if( Graph(k,i) == 1 )
            for j=1:CBBA_Params.M

                % Implement table for each task
                
                if( old_z(k,j) == k ) % Entries 1 to 4: Sender thinks he has the task
                
                    % Entry 1: Update or Leave
                    if( z(i,j) == i ) 
                        if( old_y(k,j) - y(i,j) > epsilon )  % Update
                            z(i,j) = old_z(k,j);
                            y(i,j) = old_y(k,j);
                        elseif( abs(old_y(k,j) - y(i,j)) <= epsilon )  % Equal scores
                            if( z(i,j) > old_z(k,j) )  % Tie-break based on smaller index
                                z(i,j) = old_z(k,j);
                                y(i,j) = old_y(k,j);
                            end
                        end
                        
                    % Entry 2: Update
                    elseif( z(i,j) == k ) 
                        z(i,j) = old_z(k,j);
                        y(i,j) = old_y(k,j);
                    
                    % Entry 3: Update or Leave
                    elseif( z(i,j) > 0 ) 
                        if( old_t(k,z(i,j)) > t(i,z(i,j)) )  % Update
                            z(i,j) = old_z(k,j);
                            y(i,j) = old_y(k,j);
                        elseif( (old_y(k,j) - y(i,j)) > epsilon )  % Update
                            z(i,j) = old_z(k,j);
                            y(i,j) = old_y(k,j);
                        elseif( abs(old_y(k,j) - y(i,j)) <= epsilon )  % Equal scores
                            if( z(i,j) > old_z(k,j) )  % Tie-break based on smaller index
                                z(i,j) = old_z(k,j);
                                y(i,j) = old_y(k,j);
                            end
                        end
                
                    % Entry 4: Update
                    elseif( z(i,j) == 0 )
                        z(i,j) = old_z(k,j);
                        y(i,j) = old_y(k,j);
                        
                    else
                        disp('Unknown winner value: Should not be here, please revise')
                    end
                    
                elseif( old_z(k,j) == i ) % Entries 5 to 8: Sender thinks receiver has the task

                    % Entry 5: Leave
                    if( z(i,j) == i ) 
                        % Do nothing
                        
                     % Entry 6: Reset
                    elseif( z(i,j) == k ) 
                        z(i,j) = 0;
                        y(i,j) = 0;
                   
                     % Entry 7: Reset or Leave
                    elseif( z(i,j) > 0 ) 
                        if( old_t(k,z(i,j)) > t(i,z(i,j)) )  % Reset
                            z(i,j) = 0;
                            y(i,j) = 0;
                        end
                        
                    % Entry 8: Leave
                    elseif( z(i,j) == 0 )
                        % Do nothing
                        
                    else
                        disp('Unknown winner value: Should not be here, please revise')
                    end
                   
                elseif( old_z(k,j) > 0 ) % Entries 9 to 13: Sender thinks someone else has the task
                    
                    % Entry 9: Update or Leave
                    if( z(i,j) == i ) 
                        if( old_t(k,old_z(k,j)) > t(i,old_z(k,j)) )
                            if ( (old_y(k,j) - y(i,j)) > epsilon )
                                z(i,j) = old_z(k,j);  % Update
                                y(i,j) = old_y(k,j);
                            elseif( abs(old_y(k,j) - y(i,j)) <= epsilon )  % Equal scores
                                if( z(i,j) > old_z(k,j) )  % Tie-break based on smaller index
                                    z(i,j) = old_z(k,j);
                                    y(i,j) = old_y(k,j);
                                end
                            end
                        end
                        
                     % Entry 10: Update or Reset
                    elseif( z(i,j) == k ) 
                        if( old_t(k,old_z(k,j)) > t(i,old_z(k,j)) )  % Update
                            z(i,j) = old_z(k,j);
                            y(i,j) = old_y(k,j);
                        else  % Reset
                            z(i,j) = 0;
                            y(i,j) = 0;
                        end
                        
                    % Entry 11: Update or Leave
                    elseif( z(i,j) == old_z(k,j) ) 
                        if( old_t(k,old_z(k,j)) > t(i,old_z(k,j)) )  % Update
                            z(i,j) = old_z(k,j);
                            y(i,j) = old_y(k,j);
                        end
                     
                    % Entry 12: Update, Reset or Leave
                    elseif( z(i,j) > 0 ) 
                        if( old_t(k,z(i,j)) > t(i,z(i,j)) )
                            if( old_t(k,old_z(k,j)) >= t(i,old_z(k,j)) )  % Update
                                z(i,j) = old_z(k,j);
                                y(i,j) = old_y(k,j);
                            elseif( old_t(k,old_z(k,j)) < t(i,old_z(k,j)) ) % Reset
                                z(i,j) = 0;
                                y(i,j) = 0;
                            else
                                disp('Should not be here, please revise')
                            end
                        else
                            if( old_t(k,old_z(k,j)) > t(i,old_z(k,j)) )
                                if( (old_y(k,j) - y(i,j)) > epsilon )  % Update
                                    z(i,j) = old_z(k,j);
                                    y(i,j) = old_y(k,j);
                                elseif( abs(old_y(k,j) - y(i,j)) <= epsilon )  % Equal scores
                                    if( z(i,j) > old_z(k,j) )   % Tie-break based on smaller index
                                        z(i,j) = old_z(k,j);
                                        y(i,j) = old_y(k,j);
                                    end
                                end
                            end
                        end

                    % Entry 13: Update or Leave
                    elseif( z(i,j) == 0 )
                        if( old_t(k,old_z(k,j)) > t(i,old_z(k,j)) )  % Update
                            z(i,j) = old_z(k,j);
                            y(i,j) = old_y(k,j);
                        end
                        
                    else
                        disp('Unknown winner value: Should not be here, please revise')
                    end
                    
                elseif( old_z(k,j) == 0 ) % Entries 14 to 17: Sender thinks no one has the task

                    % Entry 14: Leave
                    if( z(i,j) == i ) 
                        % Do nothing
                        
                     % Entry 15: Update
                    elseif( z(i,j) == k ) 
                        z(i,j) = old_z(k,j);
                        y(i,j) = old_y(k,j);
                   
                     % Entry 16: Update or Leave
                    elseif( z(i,j) > 0 ) 
                        if( old_t(k,z(i,j)) > t(i,z(i,j)) )  % Update
                            z(i,j) = old_z(k,j);
                            y(i,j) = old_y(k,j);
                        end
                        
                    % Entry 17: Leave
                    elseif( z(i,j) == 0 )
                        % Do nothing
                        
                    else
                        disp('Unknown winner value: Should not be here, please revise')
                    end
                    
                    % End of table
                    
                else
                    disp('Unknown winner value: Should not be here, please revise')
                end
            end
            
            % Update timestamps for all agents based on latest comm
            for n=1:CBBA_Params.N
                if( n ~= i && t(i,n) < old_t(k,n) )
                    t(i,n) = old_t(k,n);
                end
            end
            t(i,k) = T;
            
        end
    end
end

% Copy data
for n = 1:CBBA_Params.N,
    CBBA_Data(n).winners    = z(n,:);
    CBBA_Data(n).winnerBids = y(n,:);
end

end

