function T2 = EachAnalysisTable_wResidual(T,Eventstable)

%K-mediods clustering (hard)
X = table2array(T(:,2:end));
X = zscore(X)';



%Analysis of variance
X2 = [Eventstable.Month Eventstable.year];
Y = X';
N = size(Y,2);
SSQ = zeros(N,4);
Coefs = zeros(N,25);%constant, months, years
Resid = zeros(267,N);
for i = 1:N

    [p,tbl,stats,terms] = anovan(Y(:,i),X2(:,1:2),'display','off');

    SSQ(i,:)   = cell2mat(tbl(2:5,2));
    Coefs(i,:) = stats.coeffs;
    Resid(:,i) = stats.resid; 
end


T2 = array2table(Resid');
T2.Properties.VariableNames = Eventstable.SampleID;
T2.taxon = T.Properties.VariableNames(:,2:end)';
T2 = movevars(T2,'taxon','before','Newport_SRmetag_110112.NP1');
