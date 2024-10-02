clc
clear

load('E:\HCP\WM-Getm-over\Results\Gene_expression\expression_schaefer400\expression_data.mat');
%% bootstrap to get gene list
genes=gene_names; % this needs to be imported first
geneindex=1:size(genes,2);

X = expression; % 预测变量
index_nan = isnan(X(:,1));
X(index_nan,:) =[];
% 加载MRI的T-map数据或者相关数据
T_static_path = dir('E:\HCP\WM-Getm-over\Results\Gene_expression\Regress\mat\*r_zscore_new_0911.mat');
%number of bootstrap iterations:
bootnum=1000;

for num = 1:size(T_static_path,1)
    load(fullfile(T_static_path(num).folder,T_static_path(num).name));
    %print out results
    output_path = fullfile(T_static_path(num).folder,T_static_path(num).name(1:end-4));
    mkdir(output_path);


    if size(mean_r_feature_zscore,1) == 1
        mean_r_feature_zscore = mean_r_feature_zscore';  % 响应变量
    end
    Y = mean_r_feature_zscore;
    Y(index_nan) =[];
    % Do PLS in 2 dimensions (with 2 components):
    dim=8;
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);

    %store regions' IDs and weights in descending order of weight for both components:
    [R1,p1]=corr(XS,Y);
    [~,index] = max(R1);
    %align PLS components with desired direction for interpretability
    if R1(1,1)<0  %this is specific to the data shape we were using - will need ammending
        stats.W(:,1)=-1*stats.W(:,1);
        XS(:,1)=-1*XS(:,1);
    end
    if R1(index,1)<0 %this is specific to the data shape we were using - will need ammending
        stats.W(:,index)=-1*stats.W(:,index);
        XS(:,index)=-1*XS(:,index);
    end

    [PLS1w,x1] = sort(stats.W(:,1),'descend');
    PLS1ids=genes(x1);
    geneindex1=geneindex(x1);
    [PLS3w,x2] = sort(stats.W(:,index),'descend');
    PLS3ids=genes(x2);
    geneindex2=geneindex(x2);

    %save data ROI scores
    csvwrite([output_path,'\PLS1_ROIscores.csv'],XS(:,1));
    csvwrite([output_path,'\PLS3_ROIscores.csv'],XS(:,index));

    %define variables for storing the (ordered) weights from all bootstrap runs
    PLS1weights=[];
    PLS3weights=[];

    %start bootstrap
    for i=1:bootnum
        i
        myresample = randsample(size(X,1),size(X,1),1);
        res(i,:)=myresample; %store resampling out of interest
        Xr=X(myresample,:); % define X for resampled subjects
        Yr=Y(myresample,:); % define X for resampled subjects
        [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(Xr,Yr,dim); %perform PLS for resampled data

        temp=stats.W(:,1);%extract PLS1 weights
        newW=temp(x1); %order the newly obtained weights the same way as initial PLS
        if corr(PLS1w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
            newW=-1*newW;
        end
        PLS1weights=[PLS1weights,newW];%store (ordered) weights from this bootstrap run

        temp=stats.W(:,index);%extract PLS3 weights
        newW=temp(x2); %order the newly obtained weights the same way as initial PLS
        if corr(PLS3w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
            newW=-1*newW;
        end
        PLS3weights=[PLS3weights,newW]; %store (ordered) weights from this bootstrap run
    end

    %get standard deviation of weights from bootstrap runs
    PLS1sw=std(PLS1weights');
    PLS3sw=std(PLS3weights');

    %get bootstrap weights
    temp1=PLS1w./PLS1sw';
    temp2=PLS3w./PLS3sw';
    %get the p value

    temp_p1 = sum((PLS1weights./ PLS1w > 1),2) ./ size(PLS1weights,2);
    temp_p2 = sum((PLS3weights./ PLS3w > 1),2) ./ size(PLS3weights,2);
    p1_fdr = mafdr(temp_p1);
    p2_fdr = mafdr(temp_p2);
    %order bootstrap weights (Z) and names of regions
    [Z1 ind1] = sort(temp1,'descend');
    PLS1=PLS1ids(ind1);
    geneindex1 = geneindex1(ind1);
    [Z2 ind2] = sort(temp2,'descend');
    PLS3=PLS3ids(ind2);
    geneindex2 = geneindex2(ind2);

    
    %print out results
    % later use first column of these csv files for pasting into GOrilla (for
    % bootstrapped ordered list of genes)
    fid1 = fopen([output_path,'\PLS1_geneWeights.csv'],'w')
    fprintf(fid1,'%s,%s,%s\n', "Gene_names","P_adjusted","Zscore");
    for i=1:length(genes)
        fprintf(fid1,'%s, %d, %f\n', PLS1{i}, p1_fdr(i), Z1(i));
    end
    fclose(fid1)

    fid2 = fopen([output_path,'\PLS3_geneWeights.csv'],'w')
    fprintf(fid2,'%s,%s,%s\n', "Gene_names","P_adjusted","Zscore");
    for i=1:length(genes)
        fprintf(fid2,'%s, %d, %f\n', PLS3{i},p2_fdr(i), Z2(i));
    end
    fclose(fid2)
end