
function [result_data,unique_id] = test_data(Data,column)

% This code reformulates data suitable for testing

% strategy:
% First, make zeros matrix with desirable dimension
% Second, fill values
% Last, convert zero to nan.

unique_years = unique(Data.year);
unique_id = unique(Data.id);

% Setting time index

Data.time = zeros(length(Data.year),1);

for t = 1:length(unique_years)
    yr = unique_years(t);
    mask_yr = Data.year == yr;
    Data(Data.year == yr,'time') = table(repmat(t,sum(mask_yr),1));
end

% pid : id

for i = 1:length(unique_id)
    
    pid = unique_id(i);
    
    Data_pid = Data(Data.id == pid,[column, "time"]);
    
    for j = 1:length(Data_pid.time)
        
        t = Data_pid.time(j);
        
        y_data(i,t) = Data_pid{j,1};
    end
end

result_data = y_data;
