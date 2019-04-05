%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  N is the number of subjects
%  T is the number of tasks
%  x: 268*268*N*T
%  y: N*1
%  d: number of elements in lower triangular part of 268*268
%  group: 1*N of type subject
%  subject: a class of 9 task-based connectome and single label
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [results] = main()
%     data = csvread('../data/gender.csv');
    dataset = "math"; % LDA on UCLA + ages on 3 bins for HCP
    if dataset=="ucla.236"
        x = load('../data.236/all_mats_6tasks.mat');
        wais = load('../data.236/wais_raw');
        wms = load('../data.236/wms_raw');
        diagnosis = load('../data.236/diagnosis');
        g = buildGroup(x.all_mats,dataset,zeros(236,1),zeros(236,1),false); % mask=false, Bins
        options = [];
        options.thresh=0.05;
        options.seed = randi([1 10000]);
        options.k = 2;
        phenotypes = [phenotype('wais',wais) phenotype('wms',wms)];
        options.phenotypes = phenotypes;
        options.diagnosis = diagnosis;
    elseif dataset =="math"
        % x = load('/home/kailong/Scheinost-Lab/math/data/all_mats.mat');
        % y = load('/home/kailong/Scheinost-Lab/math/data/all_behav.mat'); 
        x = load('/home/kailong/Scheinost-Lab/code_from_javid/rCPM/input','all_mats');
        y = load('/home/kailong/Scheinost-Lab/code_from_javid/rCPM/input','all_behav');
        y=y.all_behav;
        g = buildGroup(x.all_mats,dataset,zeros(132,1),zeros(132,1),false); % mask=false, Bins
        options = [];
        options.thresh1=0.3;
        options.thresh2=0.1;
        options.seed = randi([1 10000]);
        options.k = 10;
        phenotypes = [phenotype('behav',y)];
        options.phenotypes = phenotypes;
        options.diagnosis = zeros(length(y));
    elseif dataset=="ucla.175.antiDepression" % 42 subjects
        x = load('../data.175/all_mats_antiD_MS.mat');
        x=x.all_mats_AD_MS;
        y=load('../data.175/remit_HAMD.mat');
        y =y.REMIT_HAMD;        
        g = buildGroup(x,dataset,zeros(42,1),zeros(42,1),false); % mask=false, Bins
        options = [];
        options.v_alpha=0; % alpha=0 for ridge , alpha=1 for lasso
        options.thresh=0.05;
        options.seed = randi([1 10000]);
        options.k = 2;
        
        phenotypes = [phenotype('HAMD',y)];
        options.phenotypes = phenotypes;
        options.diagnosis = y;
    elseif dataset=="ucla.175.antiPsychotics" % 46 subjects
        x = load('../data.175/all_mats_antiPsychotics.mat');
        x=x.all_mats_AP;
        y = load('../data.175/remit_SCZ.mat');
        y = y.REMIT_SCZ;
        g = buildGroup(x,dataset,zeros(46,1),zeros(46,1),false); % mask=false, Bins
        options = [];
        options.thresh=0.05;
        options.seed = randi([1 10000]);
        options.k = 2;
        phenotypes = [phenotype('scz',y)];
        options.phenotypes = phenotypes;
        options.diagnosis = y;

    elseif dataset=="ucla.175"
        x = load('../data.175/all_mats.mat');
        long_mem = load('../data.175/memory_long.mat');
        short_mem = load('../data.175/memory_long.mat');
        working_mem = load('../data.175/memory_long.mat');
        iq = load('../data.175/IQ.mat');
        exec_fxn = load('../data.175/executivefxn.mat');

        diagnosis = load('../data.175/groups.mat');

        g = buildGroup(x.all_mats,dataset,zeros(175,1),zeros(175,1),false); % mask=false, Bins
        options = [];
        options.thresh=0.05;
        options.seed = randi([1 10000]);
        options.k = 2;
        phenotypes = [phenotype('long_mem',long_mem.behav) phenotype('short_mem',short_mem.behav) ...
            phenotype('working_mem',working_mem.behav) phenotype('iq',iq.behav) phenotype('exec_fxn',exec_fxn.behav)];
        options.phenotypes = phenotypes;
        options.diagnosis = diagnosis.groups;
    elseif dataset =="hcp.515"
        gender = csvread('../data.515/gender.csv');
        age = csvread('../data.515/age.csv');
        x = load('../data.515/all_mats.mat');
        y = load('../data.515/all_behav.mat');
        desc.age = age;
        g = buildGroup(x.all_mats,dataset,gender,age,false); % mask=false, Bins
        options = [];
        options.thresh=0.05;
        options.seed = randi([1 10000]);
        options.k = 2;
        phenotypes = [phenotype('behav',y.all_behav)];
        options.phenotypes = phenotypes;
        options.diagnosis = gender;
    end
    m = factlysis(g,phenotypes,options);
    m.run();
%     m = manova(g,options);
    %       m = cca(g,options);
%     results = m.HottelingT();
%     results = m.run();
    %      results = g.pheno_cca(options);
end

function g = buildGroup(x,dataset,gender,age,mask)
N =size(x,3);
subjects(1,N) = subject(N);
for i=1:N
    subjects(i) = subject(x(:,:,i,:),i,dataset,gender(i),age(i),mask);
end
g = group(subjects);
end
