%% test modularity matrix of main 3 datasets

clearvars

addpath ../Network_Analysis_Functions/
addpath ../Helper_Functions/
addpath ../Networks/


%% Pearson's R matrix 
load('C:\Users\lpzmdh\Dropbox\Analyses\Networks\datasets\Retinal Ganglion Cell gene expression\RGC_Expression_Matrices.mat','Similarity')

S = Similarity.Pearson;
S(eye(size(S))==1) = 0;
B_RGC = S - expectedA(S);
[Vrgc,Ergc] = eig(B_RGC,'vector');

load Allen_Gene_Leaf
S = A;
S(eye(size(S))==1) = 0;
B_Allen = S - expectedA(S);
[Vallen,Eallen] = eig(B_Allen,'vector');

load cosyneFinalData.mat
[Data.A,Data.ixRetain,Data.Comps,Data.CompSizes] = prep_A(adjMatrix);
B_Cosyne = Data.A - expectedA(Data.A);
[Vcosyne,Ecosyne] = eig(B_Cosyne,'vector');


figure
subplot(311),
plot(Ergc,'o'); hold on
plot(find(Ergc > 0),Ergc(Ergc>0),'ro')
title(['Pos = ' num2str(sum(Ergc > 0))])
ylabel('RGC eigs')

subplot(312),
plot(Eallen,'o'); hold on
plot(find(Eallen > 0), Eallen(Eallen>0),'ro')
title(['Pos = ' num2str(sum(Eallen > 0))])
ylabel('Allen eigs')


subplot(313),
plot(Ecosyne,'o'); hold on
plot(find(Ecosyne > 0), Ecosyne(Ecosyne>0),'ro')
title(['Pos = ' num2str(sum(Ecosyne > 0))])
ylabel('COSYNE eigs')

