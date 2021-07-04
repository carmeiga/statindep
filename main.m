K = load('gromov.txt');
Y=load('enerxia.txt');

n=length(Y);

escollemos=randsample(n,10);

%K=K(escollemos,escollemos);
%Y=Y(escollemos);

params.bootForce=1;
params.shuff=5000; 
params.sigx=-1; 
params.sigy=-1;

[thresh,testStat,HSICarr] = hsicTestBoot(K,Y,0.001,params);

histogram(HSICarr)
hold on
plot([testStat testStat],[0 400],'LineWidth',2);
hold off
xlabel('HSICb(Z)') 
ylabel('Número de ocorrencias') 

;