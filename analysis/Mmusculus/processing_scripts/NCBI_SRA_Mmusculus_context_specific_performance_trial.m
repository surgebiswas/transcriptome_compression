function NCBI_SRA_Mmusculus_context_specific_performance_trial( lY, qt, coef, class, saveFile )

% Data prep parameters
NUMCOMPONENTS = 3;
PCTILECUTOFF = [1 99];
TESTPROPORTION = 0.2;
GAMMA = 5;

% OMP/Tratrain parameters
NMARKERS = 100;
PROPUNEXP = 0;

% Compute PC scores.
s = lY*coef(:,1:NUMCOMPONENTS);

% Extract out samples with the given class. 
qclass = steq(qt.class, class);

% Since the data is only partially annotated, lets see the distribution of
% PC scores for these samples and extract out unannotated samples that fall
% within the same region. Note this only makes sense if you know a priori
% that the specified class is very distinct.
sclass = s(qclass,:);
spct = prctile(sclass, PCTILECUTOFF);
qclass2 = true(size(qclass));
for i = 1 : NUMCOMPONENTS
    qclass2 = qclass2 & (s(:,i) > spct(1,i) & s(:,i) < spct(2,i));
end

% Partion samples of specified into into training and test folds. Set
% aside the test fold.
[~, ~, trainind] = partition_data(lY(qclass2,:), qt(qclass2,:), TESTPROPORTION);
o = get(qt, 'ObsNames');
oclass2 = get(qt(qclass2,:), 'ObsNames');
otest = oclass2(~trainind,:);

ktrain = ~steq(o, otest);
lytrain = lY(ktrain,:);
lytest = lY(~ktrain,:);
qttrain = qt(ktrain,:);
qttest = qt(~ktrain,:);


% Now weight samples using a gaussian kernel according to how close they
% are to samples of the specified class. We are assuming here that class
% separates in the first few principal components (e.g. tissue differences
% are visible in the first few components).
strain = s(qclass2&ktrain,:);
w = mvnpdf(s, mean(strain), GAMMA*cov(strain));


% Now begin the heavy computations.
if exist(saveFile, 'file');
    load(saveFile);
end

% Run weighted Marker OMP
if ~exist('somp', 'var')
    somp = weighted_marker_OMP( lytrain, PROPUNEXP, 'obsweights', w, 'maxfeatures', NMARKERS, 'savememory', true );
    save(saveFile, 'ktrain', 'w', 'somp');
end


% Tratrain.
if ~exist('model', 'var');
    model = tratrain(lytrain, lytrain(:, somp.S));
    save(saveFile, 'ktrain', 'w', 'somp', 'model');
end



end

