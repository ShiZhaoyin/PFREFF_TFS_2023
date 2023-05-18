%%
%The code of the paper 
%"Zhaoyin Shi, et al. Parameter-Free Robust Ensemble Framework of Fuzzy Clustering.
% IEEE TFS, DOI: 10.1109/TFUZZ.2023.3277692, May. 2023". 
%To use it, please cite it. 
% The Inputs are BaseClusterings U in the paper. In the Demo, they are the same with the
% paper and the readers can change the Input.
% The Outputs: F is the ensemble membership matrix; Result is the
% clustering metrics, from top to bottom: PFREFF, FCM(Average), FCM(Std),
% from left to right: ACC, NMI, Purity. gnd is the ground truth. J is the
% value of the objective function. Tsum is the running time.
%%
DataPath = '.\BaseClusterings\';  
Files = dir(fullfile([DataPath,'*.mat']));
for q = 1:length(Files)
    DataName = Files(q).name;          
    load([DataPath,'\',DataName]);    
    [N,~,M] = size(U);
    for i = 1:M
        [~,label] = max(U(:,:,i),[],2);
        re(i,:) = ClusteringMeasure(gnd,label);
    end
    Result = zeros(3,3);
    Result(2,:) = mean(re);
    Result(3,:) = std(re);
    Tstart = tic;
    [F,label,~,Alpha,J] = PFREFF(U);
    Tsum = toc(Tstart);
    Result(1,:) = ClusteringMeasure(gnd,label);
    save(['.\results\',DataName],'F','Result','gnd','J','Alpha','Tsum');
    U = []; Result = []; F = []; J = []; Alpha = []; T = [];
end