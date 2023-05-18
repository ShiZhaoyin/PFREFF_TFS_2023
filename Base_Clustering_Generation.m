%%
% This Script is used to generate base clusterings by FCM. 
%  You can use it to generate base clusters manually. 
% Please note that the base clusters used in the article are placed in 
% the "BaseClusterings" folder
%%
DataPath = '.\Data\'; % load raw data(.mat)
Files = dir(fullfile([DataPath,'*.mat'])); 
for q = 1:length(Files)
    DataName = Files(q).name;          
    load([DataPath,'\',DataName]);   
    k = length(unique(gnd));
    m = 1.08;   % Fuzzy exponent, free to set
    % The number of base clusterings is 30, which can be free to set
    for i = 1:30  
        [U(:,:,i),label] = MyFCM(X,k,m);
         result(i,:) = ClusteringMeasure(gnd,label);
    end
  
    Re = mean(result);
%     U = permute(U,[2 1 3]);
    save(['.\Base clusterings(New)\',DataName],'U','k','gnd');
    U = [];
end
