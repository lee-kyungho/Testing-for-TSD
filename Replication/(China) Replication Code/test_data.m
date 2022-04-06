
function [result_data,unique_id] = test_data(Data,column)

% This code reformulates data to be suitable for testing
% strategy:
% First, make zeros matrix with desirable dimension
% Second, fill values
% Last, convert zero to nan.

unique_id = unique(Data.stn);

for i = 1:length(unique_id)
    
    pid = unique_id(i);
    
    Data_pid = Data(Data.stn == pid,[column, "time"]);
    
    for j = 1:length(Data_pid.time)
        
        t = Data_pid.time(j);
        
        y_data(i,t) = Data_pid{j,1};
    end
end

result_data = y_data;
