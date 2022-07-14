% This file runs method of moments for 
% multiple enhancer files. 

clear all
warning('off')
close all

% parpool
% pctRunOnAll warning('off')

tic

rng('shuffle')

% if you have already generated the data make this true to avoid
% generating the data again
alreadyGenerated = false;

theModelList = {'','','','Dup'}

for topModelIndex = 1:length(theModelList)
% for topModelIndex = 1

    close all

    if topModelIndex == 1 || topModelIndex == 2 || topModelIndex == 3
        modelsExtra = ''
        if topModelIndex == 2
            subAddBool = true;
            superAddBool = false;
        elseif topModelIndex == 3
            subAddBool = false;
            superAddBool = true;
        else
            subAddBool = false;
            superAddBool = false;
        end
    elseif topModelIndex == 7
        subAddBool = true;
        superAddBool = false;
        modelsExtra = char(theModelList(topModelIndex))
    else
        subAddBool = false;
        superAddBool = false;
        modelsExtra = char(theModelList(topModelIndex))
    end
    % does the tfFidelity but with varying maximum values of beta1,
    % generally keep false
    megaBetaMaxBool = false;
    % for supplementary figures mean vs noise, efficiency plots of mean 
    % mRNA vs beta1, etc
    figsBool = false;
    % keep always true
    tfFidelityBool = true;

    mainFolderName = 'storage/'; 

    betaLower = 0.01;
    betaInter = 0.05;
    betaUpper1 = 1;
    betaUpper2 = 1;

    % betaMaxLower = 1.5;
    % betaMaxInter = 1;
    % betaMaxUpper = 10;

    % betaMaxRange = betaMaxLower:betaMaxInter:betaMaxUpper;
    % betaRange = betaLower:betaInter:betaMaxRange(betaMaxIndex);  

 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    %  make it always true
    MIBool = true; 
    % to only do the scatter plots with the text
    scatterOnlyBool = true;
    % use for any case where you need enhancer patterns 
    % plotted as well
    genEnhancerBool = false;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    titleExtra = modelsExtra;
    if contains(modelsExtra,'Dup')
        dupBool =  true;
    else
        dupBool = false;
    end
    % if you wish to run the SSA sims
    SSABool = false;
    % for making an optimal time-delay plot of correlation
    % between Ti and R
    timeDelayBool = false;
    % for introducing a delay to correlations and plotting the patterns
    givenDelayBool = false;
    % for parameter sensitivity analysis
    paramSenBool = false;
    % true for completely random params within 1 
    % order of magnitude
    trulyRandBool = false;

    % for fixing T1 and T2 and plotting changes over 
    % enhancer numbers - normally not needed 
    enhancerBool = false;
    % this makes the subadditivity at the level 
    % of r instead of kon and koff
    addRSpecialBool = false;

    % for doing parameter analysis with sub super 
    % at different levels
    iterSubSuperBool = false;  
    % for calculating fidelity across a variety of beta1
    % and beta2 ranges
    efficiencyBool = false;

    % this just stores the mean mRNA
    newFidelityBool = false;
    % to move text in text scatter to the side of 
    % a plotted dot
    moverBool = false; 

    if megaBetaMaxBool
        efficiencyBool = true;
        tfFidelityBool = true;
    end 

    if newFidelityBool || tfFidelityBool
        efficiencyBool = true;
    end



    if alreadyGenerated
        load('savedNewTfT1.mat')
        load('savedNewFidT1.mat')
        load('savedFinalTfNewT1MI.mat')
        load('savedNewFidMatT1.mat')
        load('savedNoiseTfMatT1.mat')
        load('savedNewFidMatStdT1.mat')
        load('savedEfficiencyMatT1.mat')
        load('savedMeanEfficiencyT1.mat')
        load('savedEfficiencyMatT2Fixed.mat')
        load('savedNewTfT2.mat')
        load('savedNewFidT2.mat')
        load('savedFinalTfNewT2MI.mat')
        load('savedNewFidMatT2.mat')
        load('savedNoiseTfMatT2.mat')
        load('savedNewFidMatStdT2.mat')
        load('savedEfficiencyMatT2.mat')
        load('savedMeanEfficiencyT2.mat')
        load('savedEfficiencyMatT1Fixed.mat')

        load('savedT1Mat.mat')
        load('savedT2Mat.mat')

        betaRange1 = betaLower:betaInter:betaUpper1;  
        betaRange2 = betaLower:betaInter:betaUpper2;   

        if megaBetaMaxBool
            load('savedBetaMaxFidT1.mat')
            load('savedBetaMaxFidT2.mat')
            load('savedBetaMaxFidT1MI.mat')
            load('savedBetaMaxFidT2MI.mat')
            load('savedBetaMaxRange')
        end
    end


    subParamsSingle = [0.75, 0.04, 0.02];
    superParamsSingle = [0.01,0.4,0.3];
    % matrix of param iters for super/sub
    % each row is a different param:
    % 1st row: plus, 2nd: minus1, 3rd minus2
    subParams = [0.05, 0.75, 1.5; ...
                 0.005, 0.04, 0.08; ...
                 0.003, 0.02, 0.04]; 
    superParams = [0.01, 0.15, 0.32; ...
                   0.05, 0.2, 0.4; ...
                   0.03, 0.15, 0.3];

    % values for reference
    % r1 = 120;
    % r2 = 140;
    % params for sub/superadditivity cases for r instead 
    % of kon/koff
    subR1 = 20;
    subR2 = 25;

    superR1 = 20;
    superR2 = 25;
    % for saving the scatter files
    wordsPlus = {'lowPlus','medPlus','highPlus'};
    wordsMinus = {'lowMinus','medMinus','highMinus'};



    % for SSA
    numRuns = 1000;
    % timespan of moment closure method or Gillespie
    numTimePoints = 100;
    t = linspace(0,30,numTimePoints);

    % there are 100 entries in t so makes sense to make the time-delay
    % array up to around ~25.
    timeDelayArray = 0:numTimePoints/4;

    % this is separate from the time-delay calculation. 
    % this actually does all the patters with a given time
    % delay
    givenDelay = 2;


    if timeDelayBool && ~SSABool
        disp('time-delay only works with SSA')
        return
    end

    if givenDelayBool && ~SSABool
        disp('given delay only works with SSA')
        return
    end

    % determines how much the koff and kon goes up/down
    % in the super/sub additive cases when doing MM2


    figFormat = '.svg';
    tempFigFormat = 'svg';
    % font size for the legend
    lgdFontSize = 38;
    % font size 
    ticksFontSize = 60;
    % defines the fontsize for the plots
    fs = 80;
    width = 10;

    % defining colors
    T1BSColor = [238 123 23]./255;
    T2BSColor = [1 64 172]./255;

    T1BSColorSpecial = [200 80 10]./255;
    T2BSColorSpecial = [1 100 250]./255;

    TotalBSColor = [143 15 228]./255;
    enhancerNumColor = [0 0 0]./255;

    % default colors
    colors(1,:) = [0, 0.4470, 0.7410];          
    colors(2,:) = [0.8500, 0.3250, 0.0980];                 
    colors(3,:) = [0.9290, 0.6940, 0.1250];                 
    colors(4,:) = [0.4940, 0.1840, 0.5560];               
    colors(5,:) = [0.4660, 0.6740, 0.1880];               
    colors(6,:) = [0.3010, 0.7450, 0.9330];                 
    colors(7,:) = [0.6350, 0.0780, 0.1840];

    % colors for the lowest number of binding sites
    % to the highest number 

    orangeInit = [240 234 227]./255;
    blueInit = [206 208 210]./255;
    % for total BS
    % purpleInit = [199 208 193]./255;
    purpleInit = [230 230 230]./255;
    % for noise
    pinkInit = [247 235 247]./255;
    % for enhancers
    blackInit = [229 229 229]./255;

    orangeFinal = [238 123 23]./255;
    blueFinal = [1 64 172]./255;
    pinkFinal = [152 23 152]./255;
    purpleFinal = [143 15 228]./255;
    blackFinal = [0 0 0]./255;


    % 0 to 4 binding sites but T1 = 4 and T2 = 4 are never used
    % for the patterns so keeping it at 4 gives us more flexibility
    % also for pattern 3 we dont have T1 + T2 = 0
    myLength = 4;
    colors_T1 = [linspace(orangeInit(1),orangeFinal(1),myLength)', ...
                linspace(orangeInit(2),orangeFinal(2),myLength)', ...
                linspace(orangeInit(3),orangeFinal(3),myLength)'];

    colors_T2 = [linspace(blueInit(1),blueFinal(1),myLength)', ...
                linspace(blueInit(2),blueFinal(2),myLength)', ...
                linspace(blueInit(3),blueFinal(3),myLength)'];

    colors_Noise = [linspace(pinkInit(1),pinkFinal(1),myLength)', ...
                linspace(pinkInit(2),pinkFinal(2),myLength)', ...
                linspace(pinkInit(3),pinkFinal(3),myLength)'];

    colors_total = [linspace(purpleInit(1),purpleFinal(1),myLength)', ...
                linspace(purpleInit(2),purpleFinal(2),myLength)', ...
                linspace(purpleInit(3),purpleFinal(3),myLength)'];

    % 1 to 4 enhancers
    colors_Enhancers = [linspace(blackInit(1),blackFinal(1),myLength)', ...
                linspace(blackInit(2),blackFinal(2),myLength)', ...
                linspace(blackInit(3),blackFinal(3),myLength)'];

    myLength = 5;
    colors_T1Special = [linspace(orangeInit(1),orangeFinal(1),myLength)', ...
                linspace(orangeInit(2),orangeFinal(2),myLength)', ...
                linspace(orangeInit(3),orangeFinal(3),myLength)'];

    colors_T2Special = [linspace(blueInit(1),blueFinal(1),myLength)', ...
                linspace(blueInit(2),blueFinal(2),myLength)', ...
                linspace(blueInit(3),blueFinal(3),myLength)'];


    % We subtract by 2 to get rid of the pointers '.' and '..'
    modelDir = strcat('models', modelsExtra);
    % will store the fileNames without the .m at the end
    shortFileNames = cell(1);
    % getting all the files in the model directory
    fileinfo = dir(modelDir);
    fnames = {fileinfo.name};
    % getting rid of pointers '.' and '..' 
    fnames(1) = [];
    % the repeated command is not a typo, once we 
    % remove the first element in the array '.', then
    % '..' becomes the first element
    fnames(1) = [];
    tracker = 1;
    numModels = 0;
    for i = 1:length(fnames)
        fileName = fnames{i}(1:end);
        fileName = fileName(1:end-2);
        % deleting the '.m' at the end from all file names
        if contains(fileName,'modelDef_')
            shortFileNames{tracker} = fileName;
            tracker = tracker + 1;
            numModels = numModels + 1;
        end
    end

    betaMaxRange = 1;

    if paramSenBool
        % how many sets of parameters you want to analyze
        numParamSets = 10;
        itersEfficiency = 1;
    elseif efficiencyBool && ~alreadyGenerated  
        if megaBetaMaxBool
            betaMaxRange = betaMaxLower:betaMaxInter:betaMaxUpper;
            save('savedBetaMaxRange','betaMaxRange')
            % the new fidelity for all the models 
            betaMaxStorageFidT1 = zeros(numModels,length(betaMaxRange));
            betaMaxStorageFidT2 = zeros(numModels,length(betaMaxRange));

            betaMaxStorageMIT1 = zeros(numModels,length(betaMaxRange));
            betaMaxStorageMIT2 = zeros(numModels,length(betaMaxRange));        
        end 

        itersEfficiency = 2;

    else
        numParamSets = 1;
        itersEfficiency = 1;
    end


    % constant parameter for volume (ignore
    % for the most part unless volume  of 
    % the container holding the chemicals is 
    % relevant)
    kappa = 1;
    % change to 'time-dependent' for time dependent
    % propensities
    options.mode = 'constant'; 
    % change to 'concentration' for concentrations
    options.scale = 'absolute'; 

    if iterSubSuperBool
        subSuperIter = size(subParams,2)
    else
        subSuperIter = 1;
    end

    if topModelIndex > 3 && topModelIndex < 7
        itersEfficiency = 1;
    end



    % 1
    for plusIndex = 1:subSuperIter
        plusConstantKoffSub = subParams(1,plusIndex);
        plusConstantKonSuper = superParams(1,plusIndex);
        % 2
        for minusIndex = 1:subSuperIter
            minusConstantKonSub1 = subParams(2,minusIndex);
            minusConstantKonSub2 = subParams(3,minusIndex);

            minusConstantKoffSuper1 = superParams(2,minusIndex);
            minusConstantKoffSuper2 = superParams(3,minusIndex);
            % here we iterate over T1 or T2 for the efficiency
            % fidelity analysis
            % 3 

            for mainTfIndex = 1:itersEfficiency

                if efficiencyBool && ~alreadyGenerated && ~megaBetaMaxBool

                    if mainTfIndex == 1
                        betaRange1 = betaLower:betaInter:betaUpper1;  
                        betaRange = betaRange1;  
                        save('savedBetaRange1','betaRange1')
                        numParamSets = length(betaRange);
              

                        efficiencyMatT1 = zeros(numModels,length(betaRange));
                        finalEfficiencyT1 = zeros(1,numModels);
                        efficiencyMatT2fixed = zeros(numModels,length(betaRange));
                        % each row is a model and the columns are the fidelities
                        % with respect to T1 and T2 as beta1 and beta2 change
                        % *for the new fidelity, these can store mean/raw mRNA
                        % for each model. Note, however, that new fidelity is not defined
                        % for either T1 or T2 but it's just a general value so we only need
                        % one mat    
                        % one is for changing beta1 and the other for beta2
                        newFidMatT1 = zeros(numModels,length(betaRange));
                        % for plotting error bars
                        newFidStdMatT1 = zeros(numModels,length(betaRange));
                        % holds the final fidelities for all models in the order that
                        % they were iterated. These are the correlations
                        % between the mean mRNA with the TF beta1 and beta2 params
                        finalNewFidT1 = zeros(1,numModels);

                        % stores T1 and T2 vals for increasing beta1 and gamma1
                        tfFidMatT1 = zeros(numModels,length(betaRange));
                        %  stores the correlation between mean mRNA and mean TF
                        finalTfNewT1 = zeros(1,numModels);
                        %  stores the mutual information between mean mRNA and mean TF
                        finalTfNewT1MI = zeros(1,numModels);
                        % stores the noise for multiple values of beta1 and beta2
                        noiseTfMatT1 = zeros(numModels,length(betaRange));
                    else
                        betaRange2 = betaLower:betaInter:betaUpper2;  
                        betaRange = betaRange2 ;  
                        save('savedBetaRange2','betaRange2')
                        numParamSets = length(betaRange);
                        
                        efficiencyMatT2 = zeros(numModels,length(betaRange));
                        finalEfficiencyT2 = zeros(1,numModels);    
                        efficiencyMatT1fixed = zeros(numModels,length(betaRange));
                        newFidMatT2 = zeros(numModels,length(betaRange));
                        newFidStdMatT2 = zeros(numModels,length(betaRange));
                        finalNewFidT2 = zeros(1,numModels);
                        tfFidMatT2 = zeros(numModels,length(betaRange));
                        finalTfNewT2 = zeros(1,numModels);    
                        finalTfNewT2MI = zeros(1,numModels);
                        noiseTfMatT2 = zeros(numModels,length(betaRange));   
                    end 

                    % this is not for beta max
                    numBinsNOTBetaMax = floor(length(betaRange)/3);
                    numBinsMI = numBinsNOTBetaMax; 

                end     
       

                disp('you are at this TF')
                mainTfIndex
                % pause

                if ~iterSubSuperBool
                    plusConstantKoffSub = subParamsSingle(1,1);
                    minusConstantKonSub1 = subParamsSingle(1,2);
                    minusConstantKonSub2 = subParamsSingle(1,3);

                    plusConstantKonSuper = superParamsSingle(1,1);
                    minusConstantKoffSuper1 = superParamsSingle(1,2);
                    minusConstantKoffSuper2 = superParamsSingle(1,3);
                end

                % for analyzing the effects of betaMax 
                % remember that betaRange is 0.01:0.20:20;  
                % 4 going over betaMaxes
                for betaMaxIndex = 1:length(betaMaxRange)

                    if megaBetaMaxBool
                        betaRange = betaLower:betaInter:betaMaxRange(betaMaxIndex);  
                        numParamSets = length(betaRange); 

                        newFidMatT1 = zeros(numModels,length(betaRange));
                        newFidMatT2 = zeros(numModels,length(betaRange));
                        newFidStdMatT1 = zeros(numModels,length(betaRange));
                        newFidStdMatT2 = zeros(numModels,length(betaRange));
                        finalNewFidT1 = zeros(1,numModels);
                        finalNewFidT2 = zeros(1,numModels);
                        efficiencyMatT1 = zeros(numModels,length(betaRange));
                        efficiencyMatT2 = zeros(numModels,length(betaRange));
                        finalEfficiencyT1 = zeros(1,numModels);
                        finalEfficiencyT2 = zeros(1,numModels);    
                        efficiencyMatT1fixed = zeros(numModels,length(betaRange));
                        efficiencyMatT2fixed = zeros(numModels,length(betaRange));
                        tfFidMatT1 = zeros(numModels,length(betaRange));
                        tfFidMatT2 = zeros(numModels,length(betaRange));
                        finalTfNewT1 = zeros(1,numModels);
                        finalTfNewT2 = zeros(1,numModels);   
                        finalTfNewT1MI = zeros(1,numModels);
                        finalTfNewT2MI = zeros(1,numModels);  
                        noiseTfMatT1 = zeros(numModels,length(betaRange));
                        noiseTfMatT2 = zeros(numModels,length(betaRange));   
                        
                        numBinsBetaMax = floor(length(betaRange)/3); 
                        numBinsMI = numBinsBetaMax;

                    end   

                    % for efficiencyBool, the paramSet index has already been set up 
                    % to get the values betaRange
                    % (all the values of beta)
                    for paramSet = 1:numParamSets
                        if ~alreadyGenerated
                            char(theModelList(topModelIndex))
                            subAddBool
                            superAddBool
                            betaRange(paramSet)
                            betaMaxRange(betaMaxIndex)
                            disp('next param set...')
                        end

                        % storing the means
                        mRNAmatMean = cell(1,numModels); 
                        mRNAmatMeanStd = cell(1,numModels); 
                        mRNAmatMeanCV = cell(1,numModels);
                        TFmatMean = cell(1,numModels); 
                        TFmatMeanCV = cell(1,numModels);
                        mRNATFmatMeanCorr = cell(1,numModels);

                        timeDelayMat = zeros(numModels,length(timeDelayArray)); 

                        % keeps track of total BS T1, total BS T2
                        % and number of enhancers
                        bsT1Tracker = zeros(1,numModels);
                        bsT2Tracker = zeros(1,numModels);
                        bsTotalTracker = zeros(1,numModels);
                        numEnhancersTracker = zeros(1,numModels);
                        fileTagTracker = zeros(1,numModels);

                        % keeping track of how many files there are in total
                        fileCounter = 0;

                        % CRN params
                        beta1 = 0.33;
                        betam1 = 2.7;
                        beta2 = 0.29;
                        betam2 = 3.9;
                        koff1 = 1.8;
                        kon1 = 0.36;
                        koff2 = 1.5;
                        kon2 = 0.19;
                        myalpha = 1.96;
                        r1 = 120;
                        r2 = 140;

                        % for the param sen
                        paramsVector = [beta1; betam1; beta2; betam2; koff1; ...
                            kon1; koff2; kon2; myalpha; r1; r2];

                        if paramSet ~= 1 && paramSenBool

                            if trulyRandBool
                                % 10.^((i--j).*rand + -j): for generating guesses between 10^i and 
                                % 10^{-j}: for parameter sensitivity stuff
                                % beta1 = 10.^((1--2).*rand + -2);
                                % betam1 = 10.^((2--1).*rand + -1);
                                % beta2 = 10.^((1--2).*rand + -2);
                                % betam2 = 10.^((2--1).*rand + -1);
                                % koff1 = 10.^((2--1).*rand + -1);
                                % kon1 = 10.^((1--2).*rand + -2);
                                % koff2 = 10.^((2--1).*rand + -1);
                                % kon2 = 10.^((1--2).*rand + -2);

                                beta1 = 10.^((1--2).*rand + -2);
                                betam1 = 10.^((1--1).*rand + -1);
                                beta2 = 10.^((1--2).*rand + -2);
                                betam2 = 10.^((1--1).*rand + -1);
                                koff1 = 10.^((1--1).*rand + -1);
                                kon1 = 10.^((1--2).*rand + -2);
                                koff2 = 10.^((1--1).*rand + -1);
                                kon2 = 10.^((1--2).*rand + -2);
                                % myalpha = ;
                                % r1 = ;
                                % r2 = ;
                     
                            else
                                % Make the spread of the Gaussians be 20% of the a values
                                sigmas = 0.2 * paramsVector; % Also a column vector.
                                % Create the noise values that we'll add to a.
                                randomNoise = randn(length(paramsVector), 1) .* sigmas;
                                % Add noise to a to make an output column vector.
                                noisyParams = paramsVector + randomNoise;

                                beta1 = noisyParams(1);
                                betam1 = noisyParams(2);
                                beta2 = noisyParams(3);
                                betam2 = noisyParams(4);
                                koff1 = noisyParams(5);
                                kon1 = noisyParams(6);
                                koff2 = noisyParams(7);
                                kon2 = noisyParams(8);
                                % myalpha = noisyParams(9);
                                % r1 = noisyParams(10);
                                % r2 = noisyParams(11);
                            end
                        end

                        % fixing one and iterating over the other
                        if efficiencyBool && ~alreadyGenerated
                            if mainTfIndex == 1
                                beta1 = betaRange(paramSet);
                            else
                                beta2 = betaRange(paramSet);
                            end
                        end

                        % parameters relevant to TFs
                        % each row is a differnet TF
                        paramsTFArray = [beta1, betam1, beta2, betam2];

                        % simulation using second-order moment equations (MM2) 
                        % with ZC closures
                        for modelIndex = 1:numModels    

                            % stores the tf locations in the network
                            tfLocs = [];
                            % getting the number of TFs for this network. Note
                            % that the last 3 characters of every file name
                            % follow totalT1BS-totalT2BS-EnhancerNum (totalTiBS
                            % becomes TiBSperEnhancer for dup case). 
                            t1BS = str2double(shortFileNames{modelIndex}(end-2));
                            t2BS = str2double(shortFileNames{modelIndex}(end-1));
                            numEnhancers = str2double(shortFileNames{modelIndex}(end));
                            
                            bsT1Tracker(modelIndex) = t1BS;
                            bsT2Tracker(modelIndex) = t2BS;
                            bsTotalTracker(modelIndex) = t1BS + t2BS;
                            numEnhancersTracker(modelIndex) = numEnhancers;


                            if addRSpecialBool
                                koff2Temp = koff2;
                                kon2Temp = kon2;

                                koff1Temp = koff1;
                                kon1Temp = kon1; 
                                if subAddBool && numEnhancers >= 2
                                    r2Temp = r2 - subR2 * numEnhancers;
                                    r1Temp = r1 - subR1 * numEnhancers;
                                elseif superAddBool && numEnhancers >= 2
                                    r2Temp = r2 + superR2 * numEnhancers;
                                    r1Temp = r1 + superR1 * numEnhancers;                           
                                else
                                    r1Temp = r1;
                                    r2Temp = r2;               
                                end                        
                            else 

                                r1Temp = r1;
                                r2Temp = r2;  
                                if subAddBool && numEnhancers >= 2
                                    koff2Temp = koff2 + plusConstantKoffSub * numEnhancers;
                                    kon2Temp = kon2 - minusConstantKonSub2  * numEnhancers;           

                                    koff1Temp = koff1 + plusConstantKoffSub * numEnhancers;
                                    kon1Temp = kon1 - minusConstantKonSub1 * numEnhancers;
                                elseif superAddBool && numEnhancers >= 2
                                    koff2Temp = koff2 - minusConstantKoffSuper2 * numEnhancers;
                                    kon2Temp = kon2 + plusConstantKonSuper * numEnhancers;

                                    koff1Temp = koff1 - minusConstantKoffSuper1 * numEnhancers;
                                    kon1Temp = kon1 + plusConstantKonSuper * numEnhancers;              
                                else
                                    koff2Temp = koff2;
                                    kon2Temp = kon2;

                                    koff1Temp = koff1;
                                    kon1Temp = kon1;            
                                end
                            end

                            % parameters relevant to enhancers
                            % each row is a different enhancer
                            paramsEnhancerArray = [kon1Temp, koff1Temp, r1Temp, kon2Temp, koff2Temp, r2Temp];

                            % displaying name of file so we know where we are
                            fileName = shortFileNames{modelIndex};
                            disp(fileName)
                            modelDefName =  strcat(fileName);

                            % each file might have a different number 
                            % of enhancers so we have to just add the parameters
                            % that we will use. 
                            tracker = 1;
                            theta = [];
                            if t1BS == 0
                                theta = [theta; kon2Temp; koff2Temp; r2Temp];
                                currentNumTFs = 1;
                            elseif t2BS == 0
                                theta = [theta; kon1Temp; koff1Temp; r1Temp];
                                currentNumTFs = 1;            
                            else
                                theta = [theta; kon1Temp; koff1Temp; r1Temp; kon2Temp; koff2Temp; r2Temp];
                                currentNumTFs = 2;
                            end 

                            % adding myalpha in accordance to the wanted format
                            theta = [theta; myalpha];

                            % going over TFs
                            tracker = 1;
                            if t1BS == 0
                                theta = [theta; beta2; betam2];
                            elseif t2BS == 0      
                                theta = [theta; beta1; betam1];
                            else
                                theta = [theta; beta1; betam1; beta2; betam2];
                            end


                            if contains(modelDefName,'Dup')
                                tempStr = extractBefore(modelDefName, 'Dup');
                            elseif contains(modelDefName,'Adair')
                                tempStr = extractBefore(modelDefName, 'Adair');
                            else
                                tempStr = extractBefore(modelDefName,'Normal');
                            end

                            fileTag = extractAfter(tempStr,'modelDef_');
                            fileTagTracker(modelIndex) = fileTag;

                            %% Simulation using second-order moment equations (MM2) 
                            % with low dispersion closure
                            modelName = char(strcat('enhancer_',fileTag));

                            % the main difference between the SSA runs and MM2 (besides the method)
                            % is that SSA stores the means of each run while MM2 stores the means (over time)
                            % of the closure method
                            if SSABool
                                System_SSA = simulateSSA_matlab(modelDefName, ...
                                            t, theta, kappa, numRuns, options); %#ok<*UNRCH>
                                % storing the stochastic runs 
                                % rows: time
                                % columns: species
                                % z dim: number of runs

                                % holders of raw runs
                                tempmRNA = zeros(numRuns,length(t));

                                % mRNA is always at the end. storing all runs
                                % we iterate over a loop to get rid of the z dim
                                % and have one run per row
                                for runIndex = 1:numRuns
                                    tempmRNAOrig(runIndex,:) = System_SSA.sol.x(:,end,runIndex); %#ok<*SAGROW>
                                end

                                disp('the mean should be:')
                                mean(tempmRNAOrig, 1)

                                disp('the variance should be:')
                                mean(tempmRNAOrig.^2, 1) - mean(tempmRNAOrig, 1).^2

                                disp('or with the matlab func:')
                                var(tempmRNAOrig, 1)
                                
                                disp('according to CERENA the variance is:')
                                System_SSA.sol.var_y
                                pause

                                % calculating mRNA mean across timepoints (runs) of the mean of runs
                                mRNAmatMean{modelIndex} = mean(mean(tempmRNAOrig, 1));
                                
                                mRNAmatMeanStd{modelIndex} = mean(std(tempmRNAOrig, 0, 1));

                                % calculating mean std from std of across traces (rows)
                                mRNAmatMeanCV{modelIndex} = mRNAmatMeanStd{modelIndex}/mRNAmatMean{modelIndex};
                                % getting rid of nans   
                                mRNAmatMeanCV{modelIndex}(isnan(mRNAmatMeanCV{modelIndex})) = 0;

                                % getting the locations of the TFs in this network. TFs are located
                                % between the enhancer states and mRNA in the species array.
                                % The TFs are arranged in increasing order. tfLocs stores
                                % the locations for T1,T2,...Tn in indices 1,2,....n
                                tempTracker = 1;
                                for loc = currentNumTFs:-1:1
                                    tfLocs(tempTracker) = length(System_SSA.state.variable) - loc; 
                                    tempTracker = tempTracker + 1; 
                                end

                                % storing mean and variances for each TF
                                % i indexes the files, so each element in the cell array
                                tempTf = cell(1,length(tfLocs));
                                corrHolder = zeros(1,length(tfLocs));
                                miHolder = zeros(1,length(tfLocs));

                                % we invert the indices compared to MM2 since we are storing
                                % many runs. Each cell index corresponds to a different TF 
                                for myTF = 1:length(tfLocs)
                                    % storing each run for the current TF
                                    % each row is a run
                                    for runIndex = 1:numRuns
                                        tempTfOrig{myTF}(runIndex,:) = ...
                                                    System_SSA.sol.x(:,tfLocs(myTF),runIndex);
                                    end

                                    % calculating the mean for the current TF across runs to obtain
                                    % the mean TF trace
                                    TFmatMean{modelIndex}(myTF) = mean(mean(tempTfOrig{myTF},1));

                                    % calculating the CV for the current TF
                                    holder = mean(std(tempTfOrig{myTF}, 0, 1))/TFmatMean{modelIndex}(myTF);
                                    % getting rid of nans   
                                    holder(isnan(holder)) = 0;
                                    TFmatMeanCV{modelIndex}(myTF) = holder;

                                    if givenDelayBool
                                        tempCorr = zeros(1,numRuns);
                                        for runIndex = 1:numRuns
                                            tfVecHolder = tempTfOrig{myTF}(runIndex,:)';
                                            mRNAVecHolder = tempmRNAOrig(runIndex,:)';
                                            % seems like no need to interpolate 
                                            % since the runs already have t in them
                                            tempHolder = corrcoef(...
                                                    tfVecHolder(1:(end - givenDelay)),...
                                                    mRNAVecHolder(1 + givenDelay:end));
                                            tempCorr(runIndex) = tempHolder(1,2);
                                        end
                                        corrHolder(myTF) = nanmean(tempCorr); 

                                    else
                                        % storing correlation between the current TF runs and 
                                        % the respective mRNA runs.
                                        % for now doing mRNA in general and not specific to the 
                                        % enhancer that binds to this TF.
                                        % tempTf and tempmRNA are the same size with each row
                                        % storing a run and each column storing a timepoint.
                                        % we use corr on the transpose of these matrices
                                        % since corr takes the pairwise correlation between
                                        % all columns (but our runs are in rows so we transpose)
                                        % and returns a matrix where (i,j) is the correlation
                                        % between ith and jth column of the first and second
                                        % matrix. However, we only care about the correlation 
                                        % between run i of the current TF run and run i of the
                                        % mRNA run so we only keep the diagonal of this matrix.
                                        corrHolder(myTF) = nanmean(diag(corr(tempTfOrig{myTF}',tempmRNAOrig')));
                                    end 

                                    if MIBool
                                        % doing mutual information
                                        tempMI = zeros(1,numRuns);
                                        % storing MI for each run (correlation uses the predefined matlab fcn that
                                        % does it for matrices)
                                        for runIndex = 1:numRuns
                                            tempMI(runIndex) = mi(tempTfOrig{myTF}(runIndex,:)',tempmRNAOrig(runIndex,:)', numBinsMI); 
                                        end
                                        % taking the mean MI
                                        miHolder(myTF) = nanmean(tempMI);
                                    end

                                    % doing it a bit different than with the corr matrix
                                    if timeDelayBool
                                        % we try many time delays to find the optimal
                                        for delayIndex = 1:length(timeDelayArray)
                                            tempCorr = zeros(1,numRuns);
                                            for runIndex = 1:numRuns
                                                tfVecHolder = tempTfOrig{myTF}(runIndex,:)';
                                                mRNAVecHolder = tempmRNAOrig(runIndex,:)';
                                                delayHolder = timeDelayArray(delayIndex);
                                                % seems like no need to interpolate 
                                                % since the runs already have t in them
                                                tempHolder = corrcoef(...
                                                        tfVecHolder(1:(end-delayHolder)),...
                                                        mRNAVecHolder(1 + delayHolder:end));
                                                tempCorr(runIndex) = tempHolder(1,2);
                                            end
                                            % each row is a model and the columns are different
                                            % time delays
                                            timeDelayMat(modelIndex,delayIndex) = nanmean(tempCorr); 
                                        end
                                    end

                                end

                                % doing this so it looks the same as the results from MM:
                                % T1       T2      T3      ......
                                % [corrT1 corrT2  corrT3   ......] 
                                if MIBool
                                    % mean mutual information in the same from as correlation
                                    mRNATFmatMeanCorr{modelIndex} = miHolder;
                                else
                                    mRNATFmatMeanCorr{modelIndex} = corrHolder;
                                end

                            % end of SSA statement
                            else

                                % MEC stands for the “central moment equations”, 
                                % and should be kept unchanged.
                                % XO = 2 specifies the truncation order for the 
                                % moment equations of the state variables.
                                % LD for low dispersion closure technique. 
                                % ZC for zero cumulants closure technique
                                % YO = 2 specifies the order of the output moments
                                % c specifies that the moment equations for the 
                                % concentration of species should be derived.   
                                % method = 'MEC_2_LD_2_c';
                                method = 'MEC_2_ZC_2_c';

                                % checking if the system for this file has already been compiled
                                if exist(strcat(pwd,'/',modelDefName,'.mat'),'file') == 2
                                    currentFileAlreadyGenerated = true;
                                else
                                    currentFileAlreadyGenerated = false;
                                end

                                % for generating the C code for the simulation of ODE 
                                % systems, e.g., FSP, RRE, SSE, and MM
                                % Running these commands will add
                                % model-specific fields to System which contain 
                                % information about the selected modeling approach.
                                if currentFileAlreadyGenerated
                                    load(modelDefName)
                                else
                                    % modelName is the name you give the C code.
                                    % modelDefName is the name of the MATLAB file with the CRN
                                    System_MM2 = genSimFile(modelName,modelDefName,method);
                                    save(modelDefName,'System_MM2')
                                end

                                % numerical simulation of the system:
                                % the model-specific output arguments above will be 
                                % added as different fields in System.sol. In particular,
                                % System.sol.x is the simulation results for the 
                                % state variables while System.sol.y is the simulation results 
                                % for the output variables.
                                System_MM2.sol = feval(strcat('simulate_', modelName),t,theta,kappa);

                                % Indices of means and variances as given in the CERENA
                                % plot MM function
                                ind_mean_x = find(sum(System_MM2.MM.sym.state.order>=1,2) == 1);
                                ind_var_x = find((System_MM2.MM.sym.state.order(:,end-1)== ...
                                    System_MM2.MM.sym.state.order(:,end)).* ...
                                    (sum(System_MM2.MM.sym.state.order~=0,2)==2));

                                % stores the mean of the mean mRNA vectors over time 
                                % from each model.
                                mRNAmatMean{modelIndex} = mean(System_MM2.sol.x(:,ind_mean_x(end)));


                                % mRNA variance
                                tempVar = System_MM2.sol.x(:,ind_var_x(end));
                                mRNAmatMeanStd{modelIndex} = mean(sqrt(tempVar));

                                mRNAmatMeanCV{modelIndex} = mRNAmatMeanStd{modelIndex}/mRNAmatMean{modelIndex};
                                mRNAmatMeanCV{modelIndex}(isnan(mRNAmatMeanCV{modelIndex})) = 0;

                                % getting the locations of the TFs in this network. TFs are located
                                % between the enhancer states and mRNA in the species array.
                                % The TFs are arranged in increasing order. tfLocs stores
                                % the locations for T1,T2,...Tn in indices 1,2,....n
                                tempTracker = 1;
                                for loc = currentNumTFs:-1:1
                                    tfLocs(tempTracker) = length(System_MM2.state.variable) - loc; 
                                    tempTracker = tempTracker + 1; 
                                end

                                % storing mean of means and mean of variances for each TF
                                % i indexes the files, so each element in the cell array 
                                % corresponds to a different file and each row of a matrix in the
                                % cell array corresponds to a TF
                                for myTF = 1:length(tfLocs)
                                    TFmatMean{modelIndex}(myTF) = ...
                                        mean(System_MM2.sol.x(:,ind_mean_x(tfLocs(myTF))));
                                    tempTfSigma = mean(sqrt(System_MM2.sol.x(:,ind_var_x(tfLocs(myTF)))));
                                    % calculating the CV at each time point by dividing the S.D. over time
                                    holder = tempTfSigma./TFmatMean{modelIndex}(myTF);
                                    % since CV can have NaN when dividing 0/0
                                    % we change all NaNs to 0
                                    holder(isnan(holder)) = 0;
                                    TFmatMeanCV{modelIndex}(myTF) = holder;
                                end

                                % getting the correlations
                                % corrMat is a matrix consisting of the maximum correlation coefficients 
                                % across all time points. corrMatAll is a 3-dimensional matrix of all 
                                % correlation coefficients across all time points. The first
                                % two dimensions correspond to variables and the third dimension corresponds to time points.
                                % covMat is a covariance matrix corresponding to maximum correlation 
                                % coefficients across all time points.
                                % The elements in the matrix are ordered as in the CRN file.
                                %      A1  A2 . . .
                                %  A1
                                %  A2
                                %   .
                                %   .
                                %   .
                                [~,tempCorrMatAll,~] = corrmat(System_MM2,System_MM2.sol.x);
                                % we want the mean correlation between all variables across time
                                meanTempCorrMatAll = mean(tempCorrMatAll,3);

                                % constructing the vector with the correlations between
                                % mRNA with every TF. The last row of meanTempCorrMatAll
                                % should have the correlations of mRNA with every other
                                % species because mRNA is always added last. 
                                % T1       T2      T3      ......
                                % [corrT1 corrT2  corrT3   ......] 
                                mRNAlastRowCorrs =  meanTempCorrMatAll(end,:);
                                mRNATFmatMeanCorr{modelIndex} = mRNAlastRowCorrs(tfLocs); 
                            end
                        % end of looping throgh models
                        end
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                        % data analysis beginning
                        if ~alreadyGenerated
                            % storing everything from these runs into the respective column 
                            if efficiencyBool   
                                % making the vector where each entry is a model and the respective 
                                % fidelity to Ti
                                % * for new fidelity we just need one vector
                                % new fidelity is how mean mRNA changes w.r.t to beta1/beta2
                             
                                newFidHolder = zeros(1,numModels);
                                noiseHolder = zeros(1,numModels);
                        
                                % this holds the R TF corr....
                                t1Holder = zeros(1,numModels);
                                t2Holder = zeros(1,numModels);

                                % these hold the actual mean TFs....
                                newT1Holder = zeros(1,numModels);
                                newT2Holder = zeros(1,numModels);

                                newFidStdHolder = zeros(1,numModels);
                     
                                for modelIndex = 1:numModels         
                                    t1BS = str2double(shortFileNames{modelIndex}(end-2));
                                    t2BS = str2double(shortFileNames{modelIndex}(end-1));
                                    if t1BS == 0
                                        t1Holder(modelIndex) = nan;
                                        t2Holder(modelIndex) = mRNATFmatMeanCorr{modelIndex}(1);
                                        
                                        newT1Holder(modelIndex) = nan;
                                        newT2Holder(modelIndex) = TFmatMean{modelIndex}(1);

                                    elseif t2BS == 0
                                        t1Holder(modelIndex) = mRNATFmatMeanCorr{modelIndex}(1);
                                        t2Holder(modelIndex) = nan;

                                        newT1Holder(modelIndex) = TFmatMean{modelIndex}(1);
                                        newT2Holder(modelIndex) = nan;
                                    else
                                        t1Holder(modelIndex) = mRNATFmatMeanCorr{modelIndex}(1);
                                        t2Holder(modelIndex) = mRNATFmatMeanCorr{modelIndex}(2);

                                        newT1Holder(modelIndex) = TFmatMean{modelIndex}(1);
                                        newT2Holder(modelIndex) = TFmatMean{modelIndex}(2);
                                    end

                                    newFidHolder(modelIndex) = mRNAmatMean{modelIndex};
                                    newFidStdHolder(modelIndex) = mRNAmatMeanStd{modelIndex}; 
                                    
                                    noiseHolder(modelIndex) = mRNAmatMeanCV{modelIndex};

                                end 

                                % this is where  you change beta1
                                if mainTfIndex == 1 
                                     newFidMatT1(:,paramSet) = newFidHolder;
                                     newFidStdMatT1(:,paramSet) = newFidStdHolder;

                                     noiseTfMatT1(:,paramSet) = noiseHolder; 

                                     tfFidMatT1(:,paramSet) = newT1Holder;

                                     efficiencyMatT1(:,paramSet) = t1Holder;
                                     efficiencyMatT2fixed(:,paramSet) = t2Holder; 
                                % this is as you change gamma 1/beta2   
                                else
                                    newFidMatT2(:,paramSet) = newFidHolder;
                                    newFidStdMatT2(:,paramSet) = newFidStdHolder;

                                    noiseTfMatT2(:,paramSet) = noiseHolder; 

                                    tfFidMatT2(:,paramSet) = newT2Holder;
                                
                                    efficiencyMatT2(:,paramSet) = t2Holder;
                                    efficiencyMatT1fixed(:,paramSet) = t1Holder;
                                end

                                % when we are done with all parameter sets then we want to make 
                                % the usual plots 
                                if paramSet ~= numParamSets
                                    % going to the next loop since we are not plotting anything here 
                                    % except after all the loops are over
                                    continue 
                                end   
                            end

                            % building final fidelities
                            for modelIndex = 1:numModels  
                                if mainTfIndex == 1
                                    % doing the new fidelity: correlation between beta1 and mean mRNA
                                    finalHolderT1 = corrcoef(betaRange,newFidMatT1(modelIndex,:));
                                    finalNewFidT1(modelIndex) = finalHolderT1(1,2);
                                    
                                    finalHolderNewT1 = corrcoef(tfFidMatT1(modelIndex,:),newFidMatT1(modelIndex,:),'rows','pairwise');
                                    finalTfNewT1(modelIndex) = finalHolderNewT1(1,2);

                                    finalHolderNewT1 = mi(tfFidMatT1(modelIndex,:)',newFidMatT1(modelIndex,:)', numBinsMI);
                                    finalTfNewT1MI(modelIndex) = finalHolderNewT1;
                                   
                                    finalEfficiencyT1(modelIndex) = nanmean(efficiencyMatT1(modelIndex,:));
                                else
                                    finalHolderT2 = corrcoef(betaRange,newFidMatT2(modelIndex,:));
                                    finalNewFidT2(modelIndex) = finalHolderT2(1,2);

                                    % correlation between mean TF and mean mRNA while ignoring nans
                                    finalHolderNewT2 = corrcoef(tfFidMatT2(modelIndex,:),newFidMatT2(modelIndex,:),'rows','pairwise');
                                    finalTfNewT2(modelIndex) = finalHolderNewT2(1,2);

                                    finalHolderNewT2 = mi(tfFidMatT2(modelIndex,:)',newFidMatT2(modelIndex,:)', numBinsMI);
                                    finalTfNewT2MI(modelIndex) = finalHolderNewT2;

                                    finalEfficiencyT2(modelIndex) = nanmean(efficiencyMatT2(modelIndex,:));
                                end
                            end

                            if mainTfIndex == 1
                                finalNewFidT1(isnan(finalNewFidT1)) = 0;
                                finalTfNewT1(isnan(finalTfNewT1)) = 0;
                                finalTfNewT1MI(isnan(finalTfNewT1)) = 0;
                            else
                                % convert all nans to zero
                                % these are the R - beta1 correlations
                                finalNewFidT2(isnan(finalNewFidT2)) = 0;
                                % these are the mean R - mean T1 correlations
                                finalTfNewT2(isnan(finalTfNewT2)) = 0;
                                % these are the mean R - mean T1 MIs
                                finalTfNewT2MI(isnan(finalTfNewT2)) = 0;
                            end

                          % storing the fidelity for this given values of betaMax
                            if megaBetaMaxBool
                                if mainTfIndex == 1
                                    betaMaxStorageFidT1(:,betaMaxIndex) = finalTfNewT1;
                                    betaMaxStorageMIT1(:,betaMaxIndex) = finalTfNewT1MI;
                                else
                                    betaMaxStorageFidT2(:,betaMaxIndex) = finalTfNewT2;
                                    betaMaxStorageMIT2(:,betaMaxIndex) = finalTfNewT2MI;    
                                end
                            end


                            if topModelIndex == 1
                                  modelsExtra = 'ADD';
                            end
                            if topModelIndex == 2
                                    modelsExtra = 'SUBADD';
                            elseif topModelIndex == 3
                                    modelsExtra = 'SUPERADD';
                            end

                            if ~megaBetaMaxBool
                                if mainTfIndex == 1 
                                    save(strcat(mainFolderName, modelsExtra, '/saved/','savedNewTfT1'),'finalTfNewT1')
                                    save(strcat(mainFolderName, modelsExtra, '/saved/','savedNewFidT1'),'finalNewFidT1')
                                    save(strcat(mainFolderName, modelsExtra, '/saved/','savedFinalTfNewT1MI'),'finalTfNewT1MI')
                                    save(strcat(mainFolderName, modelsExtra, '/saved/','savedNewFidMatT1'),'newFidMatT1')
                                    save(strcat(mainFolderName, modelsExtra, '/saved/','savedNoiseTfMatT1'),'noiseTfMatT1')
                                    save(strcat(mainFolderName, modelsExtra, '/saved/','savedNewFidMatStdT1'),'newFidStdMatT1')
                                    save(strcat(mainFolderName, modelsExtra, '/saved/','savedEfficiencyMatT1'),'efficiencyMatT1')
                                    save(strcat(mainFolderName, modelsExtra, '/saved/','savedMeanEfficiencyT1'),'finalEfficiencyT1')
                                    save(strcat(mainFolderName, modelsExtra, '/saved/','savedEfficiencyMatT2Fixed'),'efficiencyMatT2fixed')
                                    save(strcat(mainFolderName, modelsExtra, '/saved/','savedT1Mat'),'tfFidMatT1')

                                else
                                    save(strcat(mainFolderName, modelsExtra, '/saved/', 'savedNewTfT2'),'finalTfNewT2')
                                    save(strcat(mainFolderName, modelsExtra, '/saved/', 'savedNewFidT2'),'finalNewFidT2')
                                    save(strcat(mainFolderName, modelsExtra, '/saved/', 'savedFinalTfNewT2MI'),'finalTfNewT2MI')
                                    save(strcat(mainFolderName, modelsExtra, '/saved/', 'savedNewFidMatT2'),'newFidMatT2')
                                    save(strcat(mainFolderName, modelsExtra, '/saved/', 'savedNoiseTfMatT2'),'noiseTfMatT2')
                                    save(strcat(mainFolderName, modelsExtra, '/saved/', 'savedNewFidMatStdT2'),'newFidStdMatT2')
                                    save(strcat(mainFolderName, modelsExtra, '/saved/', 'savedEfficiencyMatT2'),'efficiencyMatT2')
                                    save(strcat(mainFolderName, modelsExtra, '/saved/', 'savedMeanEfficiencyT2'),'finalEfficiencyT2')
                                    save(strcat(mainFolderName, modelsExtra, '/saved/', 'savedEfficiencyMatT1Fixed'),'efficiencyMatT1fixed')
                                    save(strcat(mainFolderName, modelsExtra, '/saved/', 'savedT2Mat'),'tfFidMatT2')
                            
                                end
                            end
                       end


                       if alreadyGenerated
                            % plotting time delay results
                            if timeDelayBool 
                                for modelIndex = 1:numModels
                                    plot(timeDelayArray,timeDelayMat(modelIndex,:), ...
                                        'LineWidth',4)
                                    hold on 
                                end
                                ax = gca;
                                ax.FontSize = ticksFontSize; 
                                xlabel('\textbf{Time delay}','fontweight', ...
                                'bold','fontsize',fs,'Interpreter','latex') 
                                ylabel('\textbf{Fidelity}','fontweight', ...
                                    'bold','fontsize', fs,'Interpreter','latex')
                                % ending the program 
                                return
                            end

                            saveStuff = {mRNAmatMean, mRNAmatMeanStd, mRNAmatMeanCV, mRNATFmatMeanCorr, TFmatMean, TFmatMeanCV};
                            save('latestSave','saveStuff')

                            % plotting mRNA vs Noise
                            % making vectors of the CV (y axis)
                            % and the means (x axis)

                            meansVec = zeros(1,length(mRNAmatMean));
                            noiseVec = zeros(1,length(mRNAmatMean));
                            stdVec = zeros(1,length(mRNAmatMean));
                            corrT1Vec = zeros(1,length(mRNAmatMean));
                            corrT2Vec = zeros(1,length(mRNAmatMean));

                            meansT1Vec = zeros(1,length(mRNAmatMean));
                            meansT2Vec = zeros(1,length(mRNAmatMean));

                            tracker = 1;

                            for i = 1:length(mRNAmatMean)
                                meansVec(tracker) = mRNAmatMean{i};
                                noiseVec(tracker) = mRNAmatMeanCV{i};
                                stdVec(tracker) = mRNAmatMeanStd{i};
                                if bsT2Tracker(i) == 0
                                    corrT1Vec(tracker) = mRNATFmatMeanCorr{i}(1);
                                    corrT2Vec(tracker) = 0;

                                    meansT1Vec(tracker) = TFmatMean{i}(1);  
                                    meansT2Vec(tracker) = 0;  
                                elseif bsT1Tracker(i) == 0
                                    corrT1Vec(tracker) = 0;
                                    corrT2Vec(tracker) = mRNATFmatMeanCorr{i}(1);

                                    meansT1Vec(tracker) = 0;  
                                    meansT2Vec(tracker) = TFmatMean{i}(1); 
                                else
                                    corrT1Vec(tracker) = mRNATFmatMeanCorr{i}(1);
                                    corrT2Vec(tracker) = mRNATFmatMeanCorr{i}(2);

                                    meansT1Vec(tracker) = TFmatMean{i}(1);
                                    meansT2Vec(tracker) = TFmatMean{i}(2);
                                end        
                                tracker = tracker + 1;
                            end
                            
                            % plotting
                            figTracker = 1;

                            if figsBool && ~efficiencyBool
                                [meansVec, firstsortorder] = sort(meansVec);

                                noiseVec = noiseVec(firstsortorder);
                                enhancerVec = numEnhancersTracker(firstsortorder);
                                totalVec = bsTotalTracker(firstsortorder);
                                stdVec = stdVec(firstsortorder);
                                % fileTagVecM = fileTagTracker(firstsortorder);
                                % t1VecM = bsT1Tracker(firstsortorder);
                                % t2VecM = bsT2Tracker(firstsortorder);

                                figure(figTracker)
                                figTracker = figTracker + 1;

                                for i = 1:length(meansVec)
                                    s = scatter(meansVec(i),noiseVec(i),'LineWidth',width);
                                    s.MarkerEdgeColor = colors_Enhancers(enhancerVec(i),:);
                                    s.MarkerFaceColor = colors_Enhancers(enhancerVec(i),:);
                                    hold on
                                end

                                ax = gca;
                                ax.FontSize = ticksFontSize; 
                                xlabel('\textbf{Mean mRNA}','fontweight', ...
                                'bold','fontsize',fs,'Interpreter','latex') 
                                ylabel('\textbf{Noise}','fontweight', ...
                                    'bold','fontsize', fs,'Interpreter','latex')
                                tempFigName = string(strcat(mainFolderName, modelsExtra, '/meanVsNoiseE', titleExtra, ...
                                                 num2str(numParamSets), figFormat));
                                set(gcf, 'Position', get(0, 'Screensize'));      
                                saveas(gcf,char(tempFigName),tempFigFormat)

                                figure(figTracker)
                                figTracker = figTracker + 1;

                                for i = 1:length(meansVec)
                                    s = scatter(meansVec(i),noiseVec(i),'LineWidth',width);
                                    % add 1 since total BS for T1 can be zero
                                    s.MarkerEdgeColor = colors_total(totalVec(i),:);
                                    s.MarkerFaceColor = colors_total(totalVec(i),:);
                                    hold on
                                end

                                ax = gca;
                                ax.FontSize = ticksFontSize; 
                                xlabel('\textbf{Mean mRNA}','fontweight', ...
                                'bold','fontsize',fs,'Interpreter','latex') 
                                ylabel('\textbf{Noise}','fontweight', ...
                                    'bold','fontsize', fs,'Interpreter','latex')
                                tempFigName = string(strcat(mainFolderName, modelsExtra, '/meanVsNoiseT1', titleExtra, ...
                                                 num2str(numParamSets), figFormat));
                                set(gcf, 'Position', get(0, 'Screensize'));      
                                saveas(gcf,char(tempFigName),tempFigFormat)

                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                figure(figTracker)
                                figTracker = figTracker + 1;

                                for i = 1:length(stdVec)
                                    s = scatter(meansVec(i), stdVec(i),'LineWidth',width);
                                    s.MarkerEdgeColor = colors_Enhancers(enhancerVec(i),:);
                                    s.MarkerFaceColor = colors_Enhancers(enhancerVec(i),:);
                                    hold on
                                end

                                ax = gca;
                                ax.FontSize = ticksFontSize; 
                                xlabel('\textbf{Mean mRNA}','fontweight', ...
                                'bold','fontsize',fs,'Interpreter','latex') 
                                ylabel('\textbf{$\sigma$ mRNA}','fontweight', ...
                                    'bold','fontsize', fs,'Interpreter','latex')
                                tempFigName = string(strcat(mainFolderName, modelsExtra, '/meanvsstdE', titleExtra, figFormat));
                                set(gcf, 'Position', get(0, 'Screensize')); 
                                saveas(gcf,char(tempFigName),tempFigFormat)
                                hold on 


                                figure(figTracker)
                                figTracker = figTracker + 1; %#ok<*NASGU>

                                for i = 1:length(stdVec)
                                    s = scatter(meansVec(i),stdVec(i),'LineWidth',width);
                                    % add 1 since total BS for T1 can be zero
                                    s.MarkerEdgeColor = colors_total(totalVec(i),:);
                                    s.MarkerFaceColor = colors_total(totalVec(i),:);
                                    hold on
                                end
                                disp('jsdf')
                                pause
                                ax = gca;
                                ax.FontSize = ticksFontSize; 
                                xlabel('\textbf{Mean mRNA}','fontweight', ...
                                'bold','fontsize',fs,'Interpreter','latex') 
                                ylabel('\textbf{$\sigma$ mRNA}','fontweight', ...
                                    'bold','fontsize', fs,'Interpreter','latex')
                                tempFigName = string(strcat(mainFolderName, modelsExtra, '/meanvsstdT1', titleExtra, figFormat));
                                set(gcf, 'Position', get(0, 'Screensize')); 
                                saveas(gcf,char(tempFigName),tempFigFormat)
                                hold on
                     
                            end

                 
                            % all models cap at 4 bs T1, 4 bs T2, and 1 to 4 enhancers
                            % the 5 is because you can have 0 1 2 3 4 total BS Ti
                            CVMat = zeros(5,5,4);
                            fidelityMatT1 = zeros(5,5,4);
                            fidelityMatT2 = zeros(5,5,4);

                            % this tracks 
                            CVMatSpecial = zeros(5,5,4);
                            fidelityMatT1Special = zeros(5,5,4);
                            fidelityMatT2Special = zeros(5,5,4);

                            % each row is a batch
                            for modelIndex = 1:numModels
                                t1BS = str2double(shortFileNames{modelIndex}(end-2));
                                t2BS = str2double(shortFileNames{modelIndex}(end-1));
                                numEnhancers = str2double(shortFileNames{modelIndex}(end));

                                if t1BS == 0 || t2BS == 0
                                    currentNumTFs = 1;
                                else
                                    currentNumTFs = 2;
                                end

                                % for the sorted and peppered models we store at the same place but in a 
                                % different matrix
                                if CVMat(t1BS + 1,t2BS + 1,numEnhancers) ~= 0
                                    CVMatSpecial(t1BS + 1,t2BS + 1,numEnhancers) = mRNAmatMeanCV{modelIndex};
                                else 
                                    CVMat(t1BS + 1,t2BS + 1,numEnhancers) = mRNAmatMeanCV{modelIndex};
                                end

                                % building the fidelity mats
                                if newFidelityBool
                                    if currentNumTFs == 1
                                        if t1BS == 0
                                            if fidelityMatT2(t1BS + 1,t2BS + 1,numEnhancers) ~= 0
                                                fidelityMatT2Special(t1BS + 1,t2BS + 1,numEnhancers) = finalNewFidT2(modelIndex);
                                            else
                                                fidelityMatT2(t1BS + 1,t2BS + 1,numEnhancers) = finalNewFidT2(modelIndex);
                                            end
                                        elseif t2BS == 0
                                            if fidelityMatT1(t1BS + 1,t2BS + 1,numEnhancers) ~= 0
                                                fidelityMatT1Special(t1BS + 1,t2BS + 1,numEnhancers) = finalNewFidT1(modelIndex);
                                            else
                                                fidelityMatT1(t1BS + 1,t2BS + 1,numEnhancers) = finalNewFidT1(modelIndex);
                                            end
                                        else
                                            disp('something is wrong with plotting mats')
                                        end
                                    else
                                        if fidelityMatT1(t1BS + 1,t2BS + 1,numEnhancers) ~= 0  
                                            fidelityMatT1Special(t1BS + 1,t2BS + 1,numEnhancers) = finalNewFidT1(modelIndex);
                                            fidelityMatT2Special(t1BS + 1,t2BS + 1,numEnhancers) = finalNewFidT2(modelIndex); 
                                        else
                                            fidelityMatT1(t1BS + 1,t2BS + 1,numEnhancers) = finalNewFidT1(modelIndex);
                                            fidelityMatT2(t1BS + 1,t2BS + 1,numEnhancers) = finalNewFidT2(modelIndex); 
                                        end
                                    end
                                elseif tfFidelityBool
                                    % this is the last place this mat is used, 
                                    % if you use it later again this code needs to update
                                    % if MIBool
                                    %     finalTfNewT1 = finalTfNewT1MI;
                                    %     finalTfNewT2 = finalTfNewT2MI;
                                    % end

                                    if currentNumTFs == 1
                                        if t1BS == 0
                                            if fidelityMatT2(t1BS + 1,t2BS + 1,numEnhancers) ~= 0
                                                fidelityMatT2Special(t1BS + 1,t2BS + 1,numEnhancers) = finalTfNewT2(modelIndex);
                                            else
                                                fidelityMatT2(t1BS + 1,t2BS + 1,numEnhancers) = finalTfNewT2(modelIndex);
                                            end
                                        elseif t2BS == 0
                                            if fidelityMatT1(t1BS + 1,t2BS + 1,numEnhancers) ~= 0
                                                fidelityMatT1Special(t1BS + 1,t2BS + 1,numEnhancers) = finalTfNewT1(modelIndex);
                                            else
                                                fidelityMatT1(t1BS + 1,t2BS + 1,numEnhancers) = finalTfNewT1(modelIndex);
                                            end
                                        else
                                            disp('something is wrong with plotting mats')
                                        end
                                    else
                                        if fidelityMatT1(t1BS + 1,t2BS + 1,numEnhancers) ~= 0  
                                            fidelityMatT1Special(t1BS + 1,t2BS + 1,numEnhancers) = finalTfNewT1(modelIndex);
                                            fidelityMatT2Special(t1BS + 1,t2BS + 1,numEnhancers) = finalTfNewT2(modelIndex); 
                                        else
                                            fidelityMatT1(t1BS + 1,t2BS + 1,numEnhancers) = finalTfNewT1(modelIndex);
                                            fidelityMatT2(t1BS + 1,t2BS + 1,numEnhancers) = finalTfNewT2(modelIndex); 
                                        end
                                    end
                                elseif efficiencyBool
                                    if currentNumTFs == 1
                                        if t1BS == 0
                                            if fidelityMatT2(t1BS + 1,t2BS + 1,numEnhancers) ~= 0
                                                fidelityMatT2Special(t1BS + 1,t2BS + 1,numEnhancers) = finalEfficiencyT2(modelIndex);
                                            else
                                                fidelityMatT2(t1BS + 1,t2BS + 1,numEnhancers) = finalEfficiencyT2(modelIndex);
                                            end
                                        elseif t2BS == 0
                                            if fidelityMatT1(t1BS + 1,t2BS + 1,numEnhancers) ~= 0
                                                fidelityMatT1Special(t1BS + 1,t2BS + 1,numEnhancers) = finalEfficiencyT1(modelIndex);
                                            else
                                                fidelityMatT1(t1BS + 1,t2BS + 1,numEnhancers) = finalEfficiencyT1(modelIndex);
                                            end
                                        else
                                            disp('something is wrong with plotting mats')
                                        end
                                    else
                                        if fidelityMatT1(t1BS + 1,t2BS + 1,numEnhancers) ~= 0  
                                            fidelityMatT1Special(t1BS + 1,t2BS + 1,numEnhancers) = finalEfficiencyT1(modelIndex);
                                            fidelityMatT2Special(t1BS + 1,t2BS + 1,numEnhancers) = finalEfficiencyT2(modelIndex); 
                                        else
                                            fidelityMatT1(t1BS + 1,t2BS + 1,numEnhancers) = finalEfficiencyT1(modelIndex);
                                            fidelityMatT2(t1BS + 1,t2BS + 1,numEnhancers) = finalEfficiencyT2(modelIndex); 
                                        end
                                    end                        
                                else
                                    if currentNumTFs == 1
                                        if t1BS == 0
                                            if fidelityMatT2(t1BS + 1,t2BS + 1,numEnhancers) ~= 0
                                                fidelityMatT2Special(t1BS + 1,t2BS + 1,numEnhancers) = mRNATFmatMeanCorr{modelIndex}(1);
                                            else
                                                fidelityMatT2(t1BS + 1,t2BS + 1,numEnhancers) = mRNATFmatMeanCorr{modelIndex}(1);
                                            end
                                        elseif t2BS == 0
                                            if fidelityMatT1(t1BS + 1,t2BS + 1,numEnhancers) ~= 0
                                                fidelityMatT1Special(t1BS + 1,t2BS + 1,numEnhancers) = mRNATFmatMeanCorr{modelIndex}(1);
                                            else
                                                fidelityMatT1(t1BS + 1,t2BS + 1,numEnhancers) = mRNATFmatMeanCorr{modelIndex}(1);
                                            end
                                        else
                                            disp('something is wrong with plotting mats')
                                        end
                                    else
                                        if fidelityMatT1(t1BS + 1,t2BS + 1,numEnhancers) ~= 0  
                                            fidelityMatT1Special(t1BS + 1,t2BS + 1,numEnhancers) = mRNATFmatMeanCorr{modelIndex}(1);
                                            fidelityMatT2Special(t1BS + 1,t2BS + 1,numEnhancers) = mRNATFmatMeanCorr{modelIndex}(2); 
                                        else
                                            fidelityMatT1(t1BS + 1,t2BS + 1,numEnhancers) = mRNATFmatMeanCorr{modelIndex}(1);
                                            fidelityMatT2(t1BS + 1,t2BS + 1,numEnhancers) = mRNATFmatMeanCorr{modelIndex}(2); 
                                        end
                                    end
                                end
                            end


                            % generate 2D plots
                            if SSABool || subAddBool || superAddBool || genEnhancerBool
                                scatterSize = 50;
                                fsText = 12;   
                            else
                                scatterSize = 300;
                                fsText = 28;        
                            end

                            % Pattern 1:
                            % Increase TF being measured and fix the other one results 
                            % in decreasing noise and increasing fidelity to the TF being measured.

                            % Pattern 2:
                            % Fix TF being measured, and increase the other one. Decreasing 
                            % the noise and fidelity of the TF being measured. Notice that there 
                            % is a small number of exceptions involving one of the TF binding sites being 0. 

                            % Pattern 3:
                            % Increase the one being measured and decrease the other one. 
                            % Fidelity to TF being measured always goes up but noise is mixed.

                            figTracker = 100;
                            % so that the text labels appear next to the dots
                            if moverBool
                                textMoverFid = 0.009;
                                textMoverNoise = 0.007;
                            else
                                textMoverFid = 0;
                                textMoverNoise = 0;    
                            end

                            if ~paramSenBool
                                % total binding sites T1
                                for i = 1:size(fidelityMatT1,1)
                                    % total binding sites T2
                                    for j = 1:size(fidelityMatT1,2)
                                        % enhancer nums, but these shouldn't
                                        % matter when keeping total binding sites constant
                                        for k = 1:size(fidelityMatT1,3)
                                            % - text numbers on fidelity T1 vs noise 
                                            figure(figTracker)
                                            % we just set to 1 since all the models are 
                                            % equivalent when splitting 
                                            if SSABool || subAddBool || superAddBool || dupBool || genEnhancerBool
                                                currentText = strcat('[',num2str(i-1),',',num2str(j-1), ...
                                                                ',', num2str(k),']');
                                            else
                                                % the additive case of MM2 where enhancer numbers don't matter
                                                currentText = strcat('[',num2str(i-1),',',num2str(j-1),']');
                                            end

                                            % case for sorted and peppered models
                                            specialBool = false;
                                            if CVMatSpecial(i,j,k) ~= 0 && ~dupBool
                                                locXSpecial = CVMatSpecial(i,j,k);
                                                % if there is something in i,j,k of the special
                                                % CV matrix then there is something in i,j,k
                                                % of the fidelity matrix
                                                locYSpecial = fidelityMatT1Special(i,j,k);
                                                specialBool = true;
                                            end

                                            locX = CVMat(i,j,k);
                                            locY = fidelityMatT1(i,j,k);
                                            
                                            if locY ~= 0
                                                if moverBool || true
                                                    s = scatter(locX, locY, scatterSize,'filled');
                                                    s.MarkerEdgeColor = T1BSColor;
                                                    s.MarkerFaceColor = T1BSColor;  
                                                end
                                                text(locX + textMoverNoise,locY + textMoverFid, currentText,...
                                                            'Color',T1BSColor,'FontSize',fsText) 
                                                hold on
                                                if specialBool && false
                                                    s = scatter(locXSpecial, locYSpecial, scatterSize,'filled');
                                                    s.MarkerEdgeColor = T1BSColorSpecial;
                                                    s.MarkerFaceColor = T1BSColorSpecial; 
                                                    text(locXSpecial + textMoverNoise,locYSpecial + textMoverFid, currentText,...
                                                            'Color',T1BSColorSpecial,'FontSize',fsText) 
                                                end                       
                                            end
                                            hold on
                                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                            figure(figTracker + 1)
                                            % text numbers on fidelity T2 vs noise 
                                            locY = fidelityMatT2(i,j,k);

                                            if specialBool
                                                locYSpecial = fidelityMatT2Special(i,j,k);
                                            end

                                            if locY ~= 0
                                                if moverBool || true
                                                    s = scatter(locX, locY, scatterSize,'filled');
                                                    s.MarkerEdgeColor = T2BSColor;
                                                    s.MarkerFaceColor = T2BSColor; 
                                                end
                                                text(locX + textMoverNoise,locY + textMoverFid, currentText,...
                                                            'Color',T2BSColor,'FontSize',fsText)

                                                hold on
                                                if specialBool && false
                                                    s = scatter(locXSpecial, locYSpecial, scatterSize,'filled');
                                                    s.MarkerEdgeColor = T2BSColorSpecial;
                                                    s.MarkerFaceColor = T2BSColorSpecial; 
                                                    text(locXSpecial + textMoverNoise,locYSpecial + textMoverFid, currentText,...
                                                            'Color',T2BSColorSpecial,'FontSize',fsText) 
                                                end    
                                            end
                                            hold on
                                        end
                                    end
                                end
                            end


                            if MIBool

                                fidelityMatT1 = zeros(5,5,4);
                                fidelityMatT2 = zeros(5,5,4);

                 
                                fidelityMatT1Special = zeros(5,5,4);
                                fidelityMatT2Special = zeros(5,5,4);

                                % each row is a batch
                                for modelIndex = 1:numModels
                                    t1BS = str2double(shortFileNames{modelIndex}(end-2));
                                    t2BS = str2double(shortFileNames{modelIndex}(end-1));
                                    numEnhancers = str2double(shortFileNames{modelIndex}(end));

                                    if t1BS == 0 || t2BS == 0
                                        currentNumTFs = 1;
                                    else
                                        currentNumTFs = 2;
                                    end

                                    if currentNumTFs == 1
                                        if t1BS == 0
                                            if fidelityMatT2(t1BS + 1,t2BS + 1,numEnhancers) ~= 0
                                                fidelityMatT2Special(t1BS + 1,t2BS + 1,numEnhancers) = finalTfNewT2MI(modelIndex);
                                            else
                                                fidelityMatT2(t1BS + 1,t2BS + 1,numEnhancers) = finalTfNewT2MI(modelIndex);
                                            end
                                        elseif t2BS == 0
                                            if fidelityMatT1(t1BS + 1,t2BS + 1,numEnhancers) ~= 0
                                                fidelityMatT1Special(t1BS + 1,t2BS + 1,numEnhancers) = finalTfNewT1MI(modelIndex);
                                            else
                                                fidelityMatT1(t1BS + 1,t2BS + 1,numEnhancers) = finalTfNewT1MI(modelIndex);
                                            end
                                        else
                                            disp('something is wrong with plotting mats')
                                        end
                                    else
                                        if fidelityMatT1(t1BS + 1,t2BS + 1,numEnhancers) ~= 0  
                                            fidelityMatT1Special(t1BS + 1,t2BS + 1,numEnhancers) = finalTfNewT1MI(modelIndex);
                                            fidelityMatT2Special(t1BS + 1,t2BS + 1,numEnhancers) = finalTfNewT2MI(modelIndex); 
                                        else
                                            fidelityMatT1(t1BS + 1,t2BS + 1,numEnhancers) = finalTfNewT1MI(modelIndex);
                                            fidelityMatT2(t1BS + 1,t2BS + 1,numEnhancers) = finalTfNewT2MI(modelIndex); 
                                        end
                                    end   
                                end 

                                if ~paramSenBool
                                    % total binding sites T1
                                    for i = 1:size(fidelityMatT1,1)
                                        % total binding sites T2
                                        for j = 1:size(fidelityMatT1,2)
                                            % enhancer nums, but these shouldn't
                                            % matter when keeping total binding sites constant
                                            for k = 1:size(fidelityMatT1,3)
                                                % - text numbers on fidelity T1 vs noise 
                                                figure(20000)
                                                % we just set to 1 since all the models are 
                                                % equivalent when splitting 
                                                if SSABool || subAddBool || superAddBool || dupBool || genEnhancerBool
                                                    currentText = strcat('[',num2str(i-1),',',num2str(j-1), ...
                                                                    ',', num2str(k),']');
                                                else
                                                    % the additive case of MM2 where enhancer numbers don't matter
                                                    currentText = strcat('[',num2str(i-1),',',num2str(j-1),']');
                                                end

                                                % case for sorted and peppered models
                                                specialBool = false;
                                                if CVMatSpecial(i,j,k) ~= 0 && ~dupBool
                                                    locXSpecial = CVMatSpecial(i,j,k);
                                                    % if there is something in i,j,k of the special
                                                    % CV matrix then there is something in i,j,k
                                                    % of the fidelity matrix
                                                    locYSpecial = fidelityMatT1Special(i,j,k);
                                                    specialBool = true;
                                                end

                                                locX = CVMat(i,j,k);
                                                locY = fidelityMatT1(i,j,k);
                                                
                                                if locY ~= 0
                                                    if moverBool || true
                                                        s = scatter(locX, locY, scatterSize,'filled');
                                                        s.MarkerEdgeColor = T1BSColor;
                                                        s.MarkerFaceColor = T1BSColor;  
                                                    end
                                                    text(locX + textMoverNoise,locY + textMoverFid, currentText,...
                                                                'Color',T1BSColor,'FontSize',fsText) 
                                                    hold on
                                                    if specialBool && false
                                                        s = scatter(locXSpecial, locYSpecial, scatterSize,'filled');
                                                        s.MarkerEdgeColor = T1BSColorSpecial;
                                                        s.MarkerFaceColor = T1BSColorSpecial; 
                                                        text(locXSpecial + textMoverNoise,locYSpecial + textMoverFid, currentText,...
                                                                'Color',T1BSColorSpecial,'FontSize',fsText) 
                                                    end                       
                                                end
                                                hold on
                                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                figure(20001)
                                                % text numbers on fidelity T2 vs noise 
                                                locY = fidelityMatT2(i,j,k);

                                                if specialBool
                                                    locYSpecial = fidelityMatT2Special(i,j,k);
                                                end

                                                if locY ~= 0
                                                    if moverBool || true
                                                        s = scatter(locX, locY, scatterSize,'filled');
                                                        s.MarkerEdgeColor = T2BSColor;
                                                        s.MarkerFaceColor = T2BSColor; 
                                                    end
                                                    text(locX + textMoverNoise,locY + textMoverFid, currentText,...
                                                                'Color',T2BSColor,'FontSize',fsText)

                                                    hold on
                                                    if specialBool && false
                                                        s = scatter(locXSpecial, locYSpecial, scatterSize,'filled');
                                                        s.MarkerEdgeColor = T2BSColorSpecial;
                                                        s.MarkerFaceColor = T2BSColorSpecial; 
                                                        text(locXSpecial + textMoverNoise,locYSpecial + textMoverFid, currentText,...
                                                                'Color',T2BSColorSpecial,'FontSize',fsText) 
                                                    end    
                                                end
                                                hold on
                                            end
                                        end
                                    end
                                end                            
                            end            


                            if ~scatterOnlyBool

                                figTracker = figTracker + 2;
                                % vectors for pattern 3
                                fidT1Pat3 = zeros(5,5,4);
                                fidT2Pat3 = zeros(5,5,4);
                                noisePat3 = zeros(5,5,4);            

                                fidT1Pat3Special = zeros(5,5,4);
                                fidT2Pat3Special = zeros(5,5,4);
                                noisePat3Special = zeros(5,5,4);

                                % all Mats have the same size for now 
                                % total binding sites T1
                                for enhancerIndex = 1:size(fidelityMatT1,3)
                                    trackerDown = size(fidT1Pat3,2);
                                    trackerUpMain = 1;
                                    for i = 1:size(fidT1Pat3,1) 
                                        trackerUp = trackerUpMain;
                                        % total binding sites T2
                                        for j = 1:trackerDown 
                                            fidT1Pat3(trackerUp, i, enhancerIndex) = fidelityMatT1(i,j,enhancerIndex);
                                            fidT2Pat3(trackerUp, j, enhancerIndex) = fidelityMatT2(i,j,enhancerIndex);
                                            noisePat3(trackerUp, i, enhancerIndex) = CVMat(i,j,enhancerIndex);

                                            fidT1Pat3Special(trackerUp, i, enhancerIndex) = fidelityMatT1Special(i,j,enhancerIndex);
                                            fidT2Pat3Special(trackerUp, j, enhancerIndex) = fidelityMatT2Special(i,j,enhancerIndex);
                                            noisePat3Special(trackerUp, i, enhancerIndex) = CVMatSpecial(i,j,enhancerIndex);

                                            trackerUp = trackerUp + 1;
                                        end

                                        trackerDown = trackerDown - 1;
                                        trackerUpMain = trackerUpMain + 1;
                                    end
                                end

                                % making the pattern 3 for iterating over enhancer numbers
                                if subAddBool || superAddBool || dupBool || genEnhancerBool
                                    % that keep T1 + T2 fixed and iterate over enhancer/
                                    % We need to include T1 (or T2) sites in the third dimension
                                    % since multiple models may have same number of total binding sites.
                                    % Adding T1 or T2 as the third dimension reduces the ambiguity to the 
                                    % to the level of sorted and peppered models.   
                                    fidT1Pat3E = zeros(4,4,5);
                                    fidT2Pat3E = zeros(4,4,5);
                                    noisePat3E = zeros(4,4,5);                

                                    fidT1Pat3ESpecial = zeros(4,4,5);
                                    fidT2Pat3ESpecial = zeros(4,4,5);
                                    noisePat3ESpecial = zeros(4,4,5);

                                    for modelIndex = 1:numModels
                                        t1BS = str2double(shortFileNames{modelIndex}(end-2));
                                        t2BS = str2double(shortFileNames{modelIndex}(end-1));
                                        numEnhancers = str2double(shortFileNames{modelIndex}(end));
                                        totalBS = t1BS + t2BS;

                                        if noisePat3E(totalBS,numEnhancers, t1BS + 1) ~= 0
                                            noisePat3ESpecial(totalBS,numEnhancers,t1BS + 1) = mRNAmatMeanCV{modelIndex};
                                        else
                                            noisePat3E(totalBS,numEnhancers,t1BS + 1) = mRNAmatMeanCV{modelIndex};                       
                                        end

                                        if t1BS == 0
                                            if fidT2Pat3E(totalBS,numEnhancers, t1BS + 1) ~= 0
                                                fidT2Pat3ESpecial(totalBS,numEnhancers,t1BS + 1) = finalTfNewT2MI(modelIndex);
                                            else
                                                fidT2Pat3E(totalBS,numEnhancers,t1BS + 1) = finalTfNewT2MI(modelIndex);
                                            end
                                        elseif t2BS == 0
                                            if fidT1Pat3E(totalBS,numEnhancers, t1BS + 1) ~= 0
                                                fidT1Pat3ESpecial(totalBS,numEnhancers,t1BS + 1) = finalTfNewT1MI(modelIndex);
                                            else
                                                fidT1Pat3E(totalBS,numEnhancers,t1BS + 1) = finalTfNewT1MI(modelIndex);
                                            end
                                        else
                                            if fidT1Pat3E(totalBS,numEnhancers, t1BS + 1) ~= 0 
                                                fidT1Pat3ESpecial(totalBS,numEnhancers,t1BS + 1) = finalTfNewT1MI(modelIndex);
                                                fidT2Pat3ESpecial(totalBS,numEnhancers,t1BS + 1) = finalTfNewT2MI(modelIndex); 
                                            else
                                                fidT1Pat3E(totalBS,numEnhancers,t1BS + 1) = finalTfNewT1MI(modelIndex);
                                                fidT2Pat3E(totalBS,numEnhancers,t1BS + 1) = finalTfNewT2MI(modelIndex); 
                                            end
                                        end   
                                    end 
                                end

                                if subAddBool || superAddBool || dupBool || genEnhancerBool
                                    enhancerMain = size(fidelityMatT1,3);
                                else
                                    enhancerMain = 1;
                                end

                                for enhancerIndex = 1:enhancerMain
                                    % total binding sites T1
                                    % we can do all plots in 1 loop
                                    % since total bs T1 and T2 are the same size
                                    % and enhancer number does not matter
                                    origFig = figTracker;
                                    % keeps track of how many figures are added
                                    % in the next loop because we iterate over the same
                                    % set of figures many times
                                    tracker = 0;
                                    % each row is the total bs for T1 
                                    for i = 1:size(fidelityMatT1,1)
                                        figTracker =  origFig;
                                        % pattern 1 
                                        if i ~= size(fidelityMatT1,2)
                                            yPointsOrig = fidelityMatT1(:,i,enhancerIndex);
                                            yPointsOrigSpecial = fidelityMatT1Special(:,i,enhancerIndex);
                                            PatternPlotter(figTracker,yPointsOrig,yPointsOrigSpecial,colors_T2(i,:), ...
                                                width, 'enhancer', enhancerIndex, enhancerBool);
                                        end
                                        if i == 1
                                            tracker = tracker + 1;
                                        end
                                        figTracker = figTracker + 1;
                                        

                                        % pattern 2
                                        if i ~= size(fidelityMatT1,1)
                                            yPointsOrig = fidelityMatT1(i,:,enhancerIndex);
                                            yPointsOrigSpecial = fidelityMatT1Special(i,:,enhancerIndex);
                                            PatternPlotter(figTracker,yPointsOrig,yPointsOrigSpecial,colors_T1(i,:), ...
                                                width, 'enhancer', enhancerIndex, enhancerBool);
                                        end
                                        if i == 1
                                            tracker = tracker + 1;
                                        end
                                        figTracker = figTracker + 1;
                                        

                                        % pattern 3
                                        % pattern 3 gets a i - 1 in the colors 
                                        % since we don't use the first row 
                                        % which corresponds to T1 + T2 = 0.
                                        if any(fidT1Pat3(i,:,enhancerIndex))
                                            yPointsOrig = fidT1Pat3(i,:,enhancerIndex);
                                            yPointsOrigSpecial = fidT1Pat3Special(i,:);
                                            PatternPlotter(figTracker,yPointsOrig,yPointsOrigSpecial,colors_total(i-1,:),...
                                             width, 'enhancer', enhancerIndex, enhancerBool)              
                                        end
                                        if i == 1
                                            tracker = tracker + 1;
                                        end
                                        figTracker = figTracker + 1;
                                        
                                        
                                        % pattern 1
                                        if i ~= size(fidelityMatT2,1)
                                            yPointsOrig = fidelityMatT2(i,:,enhancerIndex);
                                            yPointsOrigSpecial = fidelityMatT2Special(i,:,enhancerIndex);
                                            PatternPlotter(figTracker,yPointsOrig,yPointsOrigSpecial,colors_T1(i,:), ...
                                                width, 'enhancer', enhancerIndex, enhancerBool)              
                                        end
                                        if i == 1
                                            tracker = tracker + 1;
                                        end
                                        figTracker = figTracker + 1;
                                        

                                        % pattern 2
                                        if i ~= size(fidelityMatT2,2)
                                            yPointsOrig = fidelityMatT2(:,i,enhancerIndex);
                                            yPointsOrigSpecial = fidelityMatT2Special(:,i,enhancerIndex);
                                            PatternPlotter(figTracker,yPointsOrig,yPointsOrigSpecial,colors_T2(i,:), ...
                                                width, 'enhancer', enhancerIndex, enhancerBool)              
                                        end
                                        if i == 1
                                            tracker = tracker + 1;
                                        end
                                        figTracker = figTracker + 1;
                                        

                                        % pattern 3
                                        if any(fidT1Pat3(i,:,enhancerIndex))
                                            yPointsOrig = fidT2Pat3(i,:,enhancerIndex);
                                            yPointsOrigSpecial = fidT2Pat3Special(i,:);
                                            PatternPlotter(figTracker,yPointsOrig,yPointsOrigSpecial,colors_total(i-1,:),...
                                             width, 'enhancer', enhancerIndex, enhancerBool)              
                                        end
                                        if i == 1
                                            tracker = tracker + 1;
                                        end
                                        figTracker = figTracker + 1;
                                        

                                        % pattern 1
                                        % the case when T1 = 0 (first row) and T2 = 0 (first col)
                                        % is not defined and we ignore
                                        if i ~= size(fidelityMatT1,2)       
                                            yPointsOrig = CVMat(:,i,enhancerIndex);
                                            yPointsOrigSpecial = CVMatSpecial(:,i,enhancerIndex);
                                            PatternPlotter(figTracker,yPointsOrig,yPointsOrigSpecial,colors_T2(i,:), ...
                                                width, 'enhancer', enhancerIndex, enhancerBool)              
                                        end
                                        if i == 1
                                            tracker = tracker + 1;
                                        end
                                        figTracker = figTracker + 1;
                                        

                                        % pattern 2
                                        if i ~= size(fidelityMatT1,1)
                                            yPointsOrig = CVMat(i,:,enhancerIndex);
                                            yPointsOrigSpecial = CVMatSpecial(i,:,enhancerIndex);
                                            PatternPlotter(figTracker,yPointsOrig,yPointsOrigSpecial,colors_T1(i,:), ...
                                                width, 'enhancer', enhancerIndex, enhancerBool)              
                                        end

                                        if i == 1
                                            tracker = tracker + 1;
                                        end
                                        figTracker = figTracker + 1;

                                        if any(noisePat3(i,:,enhancerIndex))
                                            % pattern 3
                                            yPointsOrig = noisePat3(i,:,enhancerIndex);
                                            yPointsOrigSpecial = noisePat3Special(i,:);
                                            PatternPlotter(figTracker,yPointsOrig,yPointsOrigSpecial,colors_total(i-1,:),...
                                             width, 'enhancer', enhancerIndex,enhancerBool)              
                                        end

                                        if i == 1
                                            tracker = tracker + 1;
                                        end
                                    % end of i loop
                                    end
                              
                                    figTracker = origFig + tracker;
                                    
                                % end of enhancer index loop
                                end

                                if subAddBool || superAddBool || dupBool || genEnhancerBool
                                    enhancerBool = true;
                                    % iterating over total T2 sites
                                    for T2Index = 1:size(fidelityMatT1,2) 
                                        % this is only because all matrices are 5 by 5
                                        T1Index = T2Index;
                                        origFig = figTracker;
                                        tracker = 0;
                                        % iterating over total T1 sites
                                        for i = 1:size(fidelityMatT1,1)
                                            figTracker = origFig;
                                            
                                            % pattern 1  
                                            yPointsOrig = fidelityMatT1(i,T2Index,:);
                                            yPointsOrigSpecial = fidelityMatT1Special(i,T2Index,:);
                                            PatternPlotter(figTracker,yPointsOrig,yPointsOrigSpecial,colors_T1Special(i,:), ...
                                                width, 'T2 sites ', T2Index - 1, enhancerBool)
                                            if i == 1
                                                tracker = tracker + 1;
                                            end
                                            figTracker = figTracker + 1;
                                            

                                            % pattern 2
                                            yPointsOrig = fidelityMatT1(T1Index,i,:);
                                            yPointsOrigSpecial = fidelityMatT1Special(T1Index,i,:);
                                            PatternPlotter(figTracker,yPointsOrig,yPointsOrigSpecial,colors_T2Special(i,:), ...
                                                width, 'T1 Sites ', T1Index - 1, enhancerBool)              
                                            if i == 1
                                                tracker = tracker + 1;
                                            end
                                            figTracker = figTracker + 1;
                                        
                                            
                                            % pattern 1      
                                            yPointsOrig = fidelityMatT2(i,T2Index,:);
                                            yPointsOrigSpecial = fidelityMatT2Special(i,T2Index,:);
                                            PatternPlotter(figTracker,yPointsOrig,yPointsOrigSpecial,colors_T1Special(i,:), ...
                                                width,  'T2 sites  ', T2Index - 1, enhancerBool)                                                                 
                                            if i == 1
                                                tracker = tracker + 1;
                                            end
                                            figTracker = figTracker + 1;
                                            

                                            % pattern 2
                                            yPointsOrig = fidelityMatT2(T1Index,i,:);
                                            yPointsOrigSpecial = fidelityMatT2Special(T1Index,i,:);
                                            PatternPlotter(figTracker,yPointsOrig,yPointsOrigSpecial,colors_T2Special(i,:), ...
                                                width, 'T1 Sites ', T1Index - 1, enhancerBool)                                               
                                            if i == 1
                                                tracker = tracker + 1;
                                            end
                                            figTracker = figTracker + 1;
                                            
                                            % pattern 1
                                            yPointsOrig = CVMat(i,T2Index,:);
                                            yPointsOrigSpecial = CVMatSpecial(i,T2Index,:);
                                            PatternPlotter(figTracker,yPointsOrig,yPointsOrigSpecial,colors_T1Special(i,:), ...
                                                width, 'T2 sites  ', T2Index - 1, enhancerBool)                                                
                                            if i == 1
                                                tracker = tracker + 1;
                                            end
                                            figTracker = figTracker + 1;
                                            

                                            % pattern 2
                                            yPointsOrig = CVMat(T1Index,i,:);
                                            yPointsOrigSpecial = CVMatSpecial(T1Index,i,:);
                                            PatternPlotter(figTracker,yPointsOrig,yPointsOrigSpecial,colors_T2Special(i,:), ...
                                                width, 'T1 Sites ', T1Index -1, enhancerBool)                                                
                                            if i == 1
                                                tracker = tracker + 1;
                                            end
                                            figTracker = figTracker + 1;
                                        % end of rows for loop
                                        end
                                       figTracker = origFig + tracker;
                                       
                                    % end of column for loop
                                    end
                                % end if subaddbool enhancer loop statement
                                end
                            % end of scatter only bool
                            end

                            figStart = 100;

                            figure(figStart)
                            

                            if iterSubSuperBool
                                subSupWords = strcat(wordsPlus{plusIndex},wordsMinus{minusIndex});
                            else
                                subSupWords = '';
                            end

                            if topModelIndex == 1
                                  modelsExtra = 'ADD';
                            end
                            if topModelIndex == 2
                                    modelsExtra = 'SUBADD';
                            elseif topModelIndex == 3
                                    modelsExtra = 'SUPERADD';
                            end

                            if ~paramSenBool
                                ax = gca;
                                ax.FontSize = ticksFontSize; 
                                xlabel('\textbf{Noise}','fontweight', ...
                                'bold','fontsize',fs,'Interpreter','latex') 
                                ylabel('\textbf{Fidelity $T_1$}','fontweight', ...
                                    'bold','fontsize', fs,'Interpreter','latex')
                                tempFigName = string(strcat(mainFolderName, modelsExtra, '/T1ScatterText',modelsExtra,  subSupWords, ...
                                titleExtra, figFormat));
                                set(gcf, 'Position', get(0, 'Screensize'));
                                saveas(gcf,char(tempFigName),tempFigFormat)
                                % to see if this fixes a bug 
                                pause(0.2)           
                                close(figure(figStart))
                                % ylim([0 0.25])
                                % xlim([1.3 2.01])
                            end
                            figStart = figStart + 1;
                            figure(figStart)
                       

                            if ~paramSenBool
                                ax = gca;
                                ax.FontSize = ticksFontSize; 
                                xlabel('\textbf{Noise}','fontweight', ...
                                'bold','fontsize',fs,'Interpreter','latex') 
                                ylabel('\textbf{Fidelity $T_2$}','fontweight', ...
                                    'bold','fontsize', fs,'Interpreter','latex')
                                tempFigName = string(strcat(mainFolderName, modelsExtra, '/T2ScatterText', modelsExtra, subSupWords, ...
                                        titleExtra, figFormat));
                                set(gcf, 'Position', get(0, 'Screensize'));
                                saveas(gcf,char(tempFigName), tempFigFormat)
                                % to see if this fixes a bug where figures are being plotted too quickly
                                pause(1)            
                                close(figure(figStart))
                                % ylim([0 0.25])
                                % xlim([1.3 2.01])
                            end
                            figStart = figStart + 1;

                            if MIBool

                                figure(20000)
                 
                                if ~paramSenBool
                                    ax = gca;
                                    ax.FontSize = ticksFontSize; 
                                    xlabel('\textbf{Noise}','fontweight', ...
                                    'bold','fontsize',fs,'Interpreter','latex') 
                                    ylabel('\textbf{Fidelity $T_1$ (MI)}','fontweight', ...
                                        'bold','fontsize', fs,'Interpreter','latex')
                                    tempFigName = string(strcat(mainFolderName, modelsExtra, '/T1ScatterTextMI',modelsExtra,  subSupWords, ...
                                    titleExtra, figFormat));
                                    set(gcf, 'Position', get(0, 'Screensize'));
                                    saveas(gcf,char(tempFigName),tempFigFormat)
                                    % to see if this fixes a bug 
                                    pause(0.2)           
                                    close(figure(20000))
                                    % ylim([0 0.25])
                                    % xlim([1.3 2.01])
                                end
                                
                                figure(20001)
                           

                                if ~paramSenBool
                                    ax = gca;
                                    ax.FontSize = ticksFontSize; 
                                    xlabel('\textbf{Noise}','fontweight', ...
                                    'bold','fontsize',fs,'Interpreter','latex') 
                                    ylabel('\textbf{Fidelity $T_2$ (MI)}','fontweight', ...
                                        'bold','fontsize', fs,'Interpreter','latex')
                                    tempFigName = string(strcat(mainFolderName, modelsExtra, '/T2ScatterTextMI', modelsExtra, subSupWords, ...
                                            titleExtra, figFormat));
                                    set(gcf, 'Position', get(0, 'Screensize'));
                                    saveas(gcf,char(tempFigName), tempFigFormat)
                                    % to see if this fixes a bug where figures are being plotted too quickly
                                    pause(1)            
                                    close(figure(20001))
                                    % ylim([0 0.25])
                                    % xlim([1.3 2.01])
                                end
                            end

                            if ~scatterOnlyBool
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                % fidelity T1
                                % for enhancerIndex = 1:enhancerMain
                                for enhancerIndex = 1
                                    if dupBool
                                        xLabel = 'Binding Sites for $T_1$/enhancer';
                                    else
                                        xLabel = 'Total Binding Sites for $T_1$';
                                    end
                                    yLabel = 'Fidelity $T_1$';
                                    tempFigName = string(strcat(mainFolderName, modelsExtra, '/T1FidelityPat1enhancers', num2str(enhancerIndex), ...
                                                    titleExtra, figFormat));
                                    figStart = LabelPutter(figStart, fs, xLabel, yLabel, tempFigName, ticksFontSize);

                                    
                                    if dupBool
                                        xLabel = 'Binding Sites for $T_2$/enhancer';
                                    else
                                        xLabel = 'Total Binding Sites for $T_2$';
                                    end                        
                                    yLabel = 'Fidelity $T_1$ ';
                                    tempFigName = string(strcat(mainFolderName, modelsExtra, '/T1FidelityPat2enhancers',num2str(enhancerIndex), ...
                                                    titleExtra, figFormat));
                                    figStart = LabelPutter(figStart, fs, xLabel, yLabel, tempFigName, ticksFontSize);


                                    if dupBool
                                       xLabel = 'Binding Sites for $T_1$/enhancer';
                                    else
                                       xLabel = 'Total Binding Sites for $T_1$';
                                    end 
                                    yLabel = 'Fidelity $T_1$';
                                    tempFigName = string(strcat(mainFolderName, modelsExtra, '/T1FidelityPat3enhancers',num2str(enhancerIndex), ...
                                                    titleExtra, figFormat));
                                    figStart = LabelPutter(figStart, fs, xLabel, yLabel, tempFigName, ticksFontSize);

                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                    % fidelity T2
                                    if dupBool
                                        xLabel = 'Binding Sites for $T_2$/enhancer';
                                    else
                                        xLabel = 'Total Binding Sites for $T_2$';
                                    end                        
                                    yLabel = 'Fidelity $T_2$';
                                    tempFigName = string(strcat(mainFolderName, modelsExtra, '/T2FidelityPat1enhancers',num2str(enhancerIndex), ...
                                                    titleExtra, figFormat));
                                    figStart = LabelPutter(figStart, fs, xLabel, yLabel, tempFigName, ticksFontSize);


                                    if dupBool
                                        xLabel = 'Binding Sites for $T_1$/enhancer';
                                    else
                                        xLabel = 'Total Binding Sites for $T_1$';
                                    end  
                                    yLabel = 'Fidelity $T_2$';
                                    tempFigName = string(strcat(mainFolderName, modelsExtra, '/T2FidelityPat2enhancers',num2str(enhancerIndex), ...
                                                    titleExtra, figFormat));
                                    figStart = LabelPutter(figStart, fs, xLabel, yLabel, tempFigName, ticksFontSize);

                                    if dupBool
                                        xLabel = 'Binding Sites for $T_2$/enhancer';
                                    else
                                        xLabel = 'Total Binding Sites for $T_2$';
                                    end  
                                    yLabel = 'Fidelity $T_2$';
                                    tempFigName = string(strcat(mainFolderName, modelsExtra, '/T2FidelityPat3enhancers',num2str(enhancerIndex), ...
                                                    titleExtra, figFormat));
                                    figStart = LabelPutter(figStart, fs, xLabel, yLabel, tempFigName, ticksFontSize);

                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                    % noise 
                                    if dupBool
                                        xLabel = 'Binding Sites for $T_1$/enhancer';
                                    else
                                        xLabel = '$T_1$';
                                    end                        
                                    yLabel = 'Noise';
                                    tempFigName = string(strcat(mainFolderName, modelsExtra, '/NoisePat1enhancers',num2str(enhancerIndex), ...
                                                    titleExtra, figFormat));
                                    figStart = LabelPutter(figStart, fs, xLabel, yLabel, tempFigName, ticksFontSize);

                                    if dupBool
                                        xLabel = 'Binding Sites for $T_2$/enhancer';
                                    else
                                        xLabel = '$T_2$';
                                    end  
                                    yLabel = 'Noise';
                                    tempFigName = string(strcat(mainFolderName, modelsExtra, '/NoisePat2enhancers',num2str(enhancerIndex), ...
                                                    titleExtra, figFormat));
                                    figStart = LabelPutter(figStart, fs, xLabel, yLabel, tempFigName, ticksFontSize);
                                   
                                    if dupBool
                                        xLabel = 'Binding Sites for $T_1$/enhancer';
                                    else
                                        xLabel = '$T_1$';
                                    end  
                                    yLabel = 'Noise';
                                    tempFigName = string(strcat(mainFolderName, modelsExtra, '/NoisePat3enhancers',num2str(enhancerIndex), ...
                                                    titleExtra, figFormat));
                                    figStart = LabelPutter(figStart, fs, xLabel, yLabel, tempFigName, ticksFontSize);
                                end

                                if subAddBool || superAddBool || dupBool || genEnhancerBool
                                   % reminder that we are using T2 index here but in reality there should be 
                                   % two loops: one for T1 and one for T2. However we can get away with only one
                                   % because the upper bound for both binding sites is the same in this case. i.e.
                                   % the matrices are 4 x 4 in the first two dimensions so we save time and make 
                                   % only one loop
                                   for T2Index = 1:size(fidelityMatT1,2)    
                                        xLabel = 'Enhancers';
                                        yLabel = 'Fidelity $T_1$';
                                        tempFigName = string(strcat(mainFolderName, modelsExtra, '/enhancerT1FidelityT2sites',num2str(T2Index - 1),titleExtra, figFormat));
                                        figStart = LabelPutter(figStart, fs, xLabel, yLabel, tempFigName, ticksFontSize);

                                        tempFigName = string(strcat(mainFolderName, modelsExtra, '/enhancerT1FidelityT1sites',num2str(T2Index - 1),titleExtra, figFormat));
                                        figStart = LabelPutter(figStart, fs, xLabel, yLabel, tempFigName, ticksFontSize);

                                        yLabel = 'Fidelity $T_2$';
                                        tempFigName = string(strcat(mainFolderName, modelsExtra, '/enhancerT2FidelityT2sites',num2str(T2Index - 1),titleExtra, figFormat));
                                        figStart = LabelPutter(figStart, fs, xLabel, yLabel, tempFigName, ticksFontSize);

                                        tempFigName = string(strcat(mainFolderName, modelsExtra, '/enhancerT2FidelityT1sites',num2str(T2Index - 1),titleExtra, figFormat));
                                        figStart = LabelPutter(figStart, fs, xLabel, yLabel, tempFigName, ticksFontSize);

                                        yLabel = 'Noise';
                                        tempFigName = string(strcat(mainFolderName, modelsExtra, '/enhancerNoiseT2sites',num2str(T2Index - 1),titleExtra, figFormat));
                                        figStart = LabelPutter(figStart, fs, xLabel, yLabel, tempFigName, ticksFontSize);

                                        yLabel = 'Noise';
                                        tempFigName = string(strcat(mainFolderName, modelsExtra, '/enhancerNoiseT1sites',num2str(T2Index - 1),titleExtra, figFormat));
                                        figStart = LabelPutter(figStart, fs, xLabel, yLabel, tempFigName, ticksFontSize);
                                   end

                                   % for T1Index = 1:size(fidT1Pat3E,3)
                                   %      xLabel = 'Enhancers';
                                   %      yLabel = 'Fidelity $T_1$';
                                   %      tempFigName = string(strcat(mainFolderName, modelsExtra, '/enhancerfidT1Pat3T1sites',num2str(T1Index - 1),titleExtra, figFormat));
                                   %      figStart = LabelPutter(figStart, fs, xLabel, yLabel, tempFigName, ticksFontSize);

                                   %      yLabel = 'Fidelity $T_2$';
                                   %      tempFigName = string(strcat(mainFolderName, modelsExtra, '/enhancerfidT2Pat3T1sites',num2str(T1Index - 1),titleExtra, figFormat));
                                   %      figStart = LabelPutter(figStart, fs, xLabel, yLabel, tempFigName, ticksFontSize);
                                       
                                   %      yLabel = 'Noise';
                                   %      tempFigName = string(strcat(mainFolderName, modelsExtra, '/enhancerNoisePat3T1sites',num2str(T1Index - 1),titleExtra, figFormat));
                                   %      figStart = LabelPutter(figStart, fs, xLabel, yLabel, tempFigName, ticksFontSize);
                                   %  end
                                % end of if statement right above
                                end
                            % end of if statement above that
                            end
                        % end of if alreadyGenerated, no need
                        % for these plots if just doing analysis
                        end
                    % 5 end looping over paramSets
                    end
                % 4 end of betaMax loop
                end
            %  3 end for iterating over T1 and 
            %  T2 for efficiency fidelity analysis
            end        
        % 2 end of plus sub super 
        end
    % 1 end of minus sub super 
    end

    if topModelIndex == 1
        modelsExtra = 'ADD';
    end

    if topModelIndex == 2
            modelsExtra = 'SUBADD';
    elseif topModelIndex == 3
            modelsExtra = 'SUPERADD';
    end


    if megaBetaMaxBool
        name1 = strcat(mainFolderName, modelsExtra, ...
            '/saved/', 'savedBetaMaxFidT1');
        name2 = strcat(mainFolderName,'/',modelsExtra, ...
            '/saved/', 'savedBetaMaxMIT1');
        name3 = strcat(mainFolderName,'/',modelsExtra, ...
            '/saved/', 'savedBetaMaxFidT2');
        name4 = strcat(mainFolderName,'/',modelsExtra, ...
            '/saved/', 'savedBetaMaxMIT2');

        save(name1,'betaMaxStorageFidT1')
        save(name2,'betaMaxStorageMIT1')
        save(name3,'betaMaxStorageFidT2')
        save(name4,'betaMaxStorageMIT2')
    end

    close all

    if ~alreadyGenerated 

        if itersEfficiency == 1 
           efficiencyMatT2 = zeros(numModels,length(betaRange));
           finalEfficiencyT2 = zeros(1,numModels);    
           efficiencyMatT1fixed = zeros(numModels,length(betaRange));
           newFidMatT2 = zeros(numModels,length(betaRange));
           newFidStdMatT2 = zeros(numModels,length(betaRange));
           finalNewFidT2 = zeros(1,numModels);
           tfFidMatT2 = zeros(numModels,length(betaRange));
           finalTfNewT2 = zeros(1,numModels);    
           finalTfNewT2MI = zeros(1,numModels);
           noiseTfMatT2 = zeros(numModels,length(betaRange));   
        end

        % betaMaxStorageFidT1 zeros(numModels,length(betaMaxRange));
        % betaMaxStorageFidT2 zeros(numModels,length(betaMaxRange));

        % betaMaxStorageMIT1 zeros(numModels,length(betaMaxRange));
        % betaMaxStorageMIT2 zeros(numModels,length(betaMaxRange));  
        if megaBetaMaxBool

            figure(1)
            % each row is a model
            for modelIndex = 1:numModels
                t1BS = str2double(shortFileNames{modelIndex}(end-2));
                if t1BS ~= 0
                    plot(betaMaxRange,betaMaxStorageFidT1(modelIndex,:),'-o', ...
                            'LineWidth',4, 'Color',colors_T1(t1BS,:))
                    hold on 
                end
            end

            ax = gca;
            ax.FontSize = ticksFontSize; 
            xlabel('Maximum value of \textbf{$\beta_1$}','fontweight', ...
            'bold','fontsize',fs,'Interpreter','latex') 
            ylabel('\textbf{Fidelity to $T_1$}','fontweight', ...
                    'bold','fontsize', fs,'Interpreter','latex')
            set(gcf, 'Position', get(0, 'Screensize'));
            tempFigName = string(strcat(mainFolderName, modelsExtra, '/betaMaxT1', modelsExtra, figFormat));
            saveas(gcf,char(tempFigName),tempFigFormat)


           figure(2)
            for modelIndex = 1:numModels
                t2BS = str2double(shortFileNames{modelIndex}(end-1));
                if t2BS ~= 0 
                    plot(betaMaxRange,betaMaxStorageFidT2(modelIndex,:),'-o','LineWidth',4, 'Color',colors_T2(t2BS,:))     
                    hold on 
                end
            end
            ax = gca;
            ax.FontSize = ticksFontSize; 
            xlabel('Maximum value of \textbf{$\gamma_1$}','fontweight', ...
            'bold','fontsize',fs,'Interpreter','latex') 
            ylabel('\textbf{Fidelity to $T_2$}','fontweight', ...
                    'bold','fontsize', fs,'Interpreter','latex')
            set(gcf, 'Position', get(0, 'Screensize'));
            tempFigName = string(strcat(mainFolderName, modelsExtra, '/betaMaxT2', modelsExtra, figFormat));
            saveas(gcf,char(tempFigName),tempFigFormat)

            figure(3)
            % each row is a model
            for modelIndex = 1:numModels
                t1BS = str2double(shortFileNames{modelIndex}(end-2));
                if t1BS ~= 0
                    plot(betaMaxRange,betaMaxStorageMIT1(modelIndex,:),'-o', ...
                            'LineWidth',4, 'Color',colors_T1(t1BS,:))
                    hold on 
                end
            end

            ax = gca;
            ax.FontSize = ticksFontSize; 
            xlabel('Maximum value of \textbf{$\beta_1$}','fontweight', ...
            'bold','fontsize',fs,'Interpreter','latex') 
            ylabel('\textbf{Fidelity to $T_1$ (MI)}','fontweight', ...
                    'bold','fontsize', fs,'Interpreter','latex')
            set(gcf, 'Position', get(0, 'Screensize'));
            tempFigName = string(strcat(mainFolderName, modelsExtra, '/betaMaxT1MI',modelsExtra, figFormat));
            saveas(gcf,char(tempFigName),tempFigFormat)


           figure(4)
            for modelIndex = 1:numModels
                t2BS = str2double(shortFileNames{modelIndex}(end-1));
                if t2BS ~= 0 
                    plot(betaMaxRange,betaMaxStorageMIT2(modelIndex,:),'-o','LineWidth',4, 'Color',colors_T2(t2BS,:))     
                    hold on 
                end
            end
            ax = gca;
            ax.FontSize = ticksFontSize; 
            xlabel('Maximum value of \textbf{$\gamma_1$}','fontweight', ...
            'bold','fontsize',fs,'Interpreter','latex') 
            ylabel('\textbf{Fidelity to $T_2$ (MI)}','fontweight', ...
                    'bold','fontsize', fs,'Interpreter','latex')
            set(gcf, 'Position', get(0, 'Screensize'));
            tempFigName = string(strcat(mainFolderName, modelsExtra, '/betaMaxT2MI',modelsExtra, figFormat));
            saveas(gcf,char(tempFigName),tempFigFormat)

        else
            % plotting the results for efficiency
            % we no longer calculate the newFidelity here, these are just lefover plots from 
            % the original efficiency plots (old fidelity vs beta1) and mean mRNA vs Beta1/2.
            % now the true new fidelity is calculated in the main loop by detecting the last 
            % iteration of the for loop
            if efficiencyBool && figsBool 

                figure(1)
                % each row is a model
                for modelIndex = 1:numModels
                    t1BS = str2double(shortFileNames{modelIndex}(end-2));
                    if t1BS ~= 0
                        % note that newFidelityT1 and efficiencyMatT1 (old fidelity) are different
                        % in this case. The new fidelity is not specific to T1 or T2 but in this case
                        % newFidelityT1 corresponds to the values of mRNA as we change beta1 
                        % values while newFidelityT2 as we change beta2 values.
                        if newFidelityBool
                            plot(betaRange1,newFidMatT1(modelIndex,:),'-o','LineWidth',4, 'Color',colors_T1(t1BS,:))               
                        else
                            plot(betaRange1,efficiencyMatT1(modelIndex,:),'-o', ...
                                'LineWidth',4, 'Color',colors_T1(t1BS,:))
                        end
                        hold on 
                    end
                end

                ax = gca;
                ax.FontSize = ticksFontSize; 
                xlabel('$T_1$ production rate \textbf{$\beta_1$}','fontweight', ...
                'bold','fontsize',fs,'Interpreter','latex') 
                if efficiencyBool && ~newFidelityBool
                    ylabel('\textbf{mRNA T1 correlation}','fontweight', ...
                        'bold','fontsize', fs -  25,'Interpreter','latex')
                else
                    ylabel('\textbf{Mean mRNA}','fontweight', ...
                        'bold','fontsize', fs,'Interpreter','latex')
                end
                set(gcf, 'Position', get(0, 'Screensize'));
                if efficiencyBool && ~newFidelityBool
                   tempFigName = string(strcat(mainFolderName, modelsExtra, '/efficiencyT1', figFormat));
                else
                    tempFigName = string(strcat(mainFolderName, modelsExtra, '/meanmRNABeta1', modelsExtra, figFormat));
                end 
                saveas(gcf,char(tempFigName),tempFigFormat)

                % making scale version where all  curves are normalized
                figure(1000)
                for modelIndex = 1:numModels
                    normalizedHolder = newFidMatT1(modelIndex,:)./max(newFidMatT1(modelIndex,:));
                    t1BS = str2double(shortFileNames{modelIndex}(end-2));
                    if t1BS ~= 0
                        plot(betaRange1,normalizedHolder,'-o','LineWidth',4, 'Color',colors_T1(t1BS,:))               
                        hold on 
                    end
                end    
                % this is only for new fidelity.
                ax = gca;
                ax.FontSize = ticksFontSize; 
                xlabel('$T_1$ production rate \textbf{$\beta_1$}','fontweight', ...
                'bold','fontsize',fs,'Interpreter','latex') 
                ylabel('\textbf{Normalized mRNA}','fontweight', ...
                    'bold','fontsize', fs-25,'Interpreter','latex')

                set(gcf, 'Position', get(0, 'Screensize'));
                tempFigName = string(strcat(mainFolderName, modelsExtra, '/newfidBeta1Normalized',modelsExtra, figFormat));
                saveas(gcf,char(tempFigName),tempFigFormat)


                figure(2)
                for modelIndex = 1:numModels
                    t2BS = str2double(shortFileNames{modelIndex}(end-1));
                    if t2BS ~= 0
                        if newFidelityBool           
                            plot(betaRange2,newFidMatT2(modelIndex,:),'-o','LineWidth',4, 'Color',colors_T2(t2BS,:))     
                        else
                            plot(betaRange2,efficiencyMatT2(modelIndex,:),'-o', ...
                                'LineWidth',4,'Color',colors_T2(t2BS,:))
                        end
                        hold on 
                    end
                end
                ax = gca;
                ax.FontSize = ticksFontSize; 
                xlabel('$T_2$ production rate \textbf{$\gamma_1$}','fontweight', ...
                'bold','fontsize',fs,'Interpreter','latex') 
                if efficiencyBool && ~newFidelityBool
                    ylabel('\textbf{mRNA T2 correlation}','fontweight', ...
                        'bold','fontsize', fs - 25,'Interpreter','latex')
                else
                    ylabel('\textbf{Mean mRNA}','fontweight', ...
                        'bold','fontsize', fs,'Interpreter','latex')
                end
                set(gcf, 'Position', get(0, 'Screensize'));
                if efficiencyBool && ~newFidelityBool
                   tempFigName = string(strcat(mainFolderName, modelsExtra, '/efficiencyT2', modelsExtra, figFormat));
                else
                    tempFigName = string(strcat(mainFolderName, modelsExtra, '/meanmRNAGamma1', modelsExtra, figFormat));
                end         

                saveas(gcf,char(tempFigName),tempFigFormat)

                % doing noise against T1
                figure(3)
                for modelIndex = 1:numModels
                    t1BS = str2double(shortFileNames{modelIndex}(end-1));
                    if t1BS ~= 0          
                        plot(betaRange1,noiseTfMatT1(modelIndex,:),'-o','LineWidth',4, 'Color',colors_T1(t1BS,:))               
                        hold on 
                    end
                end
                ax = gca;
                ax.FontSize = ticksFontSize; 
                xlabel('$T_1$ production rate \textbf{$\beta_1$}','fontweight', ...
                'bold','fontsize',fs,'Interpreter','latex') 
                ylabel('\textbf{Noise}','fontweight', ...
                        'bold','fontsize', fs - 25,'Interpreter','latex')
                set(gcf, 'Position', get(0, 'Screensize'));
                tempFigName = string(strcat(mainFolderName, modelsExtra, '/meanNoisebeta1', modelsExtra, figFormat));
                saveas(gcf,char(tempFigName),tempFigFormat)


                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if tfFidelityBool 
                    figure(200)
                    % each row is a model
                    for modelIndex = 1:numModels
                        t1BS = str2double(shortFileNames{modelIndex}(end-2));
                        if t1BS ~= 0
                            plot(betaRange1,tfFidMatT1(modelIndex,:),'-o', ...
                                    'LineWidth',4, 'Color',colors_T1(t1BS,:))
                            hold on 
                        end
                    end

                    ax = gca;
                    ax.FontSize = ticksFontSize; 
                    xlabel('$T_1$ production rate \textbf{$\beta_1$}','fontweight', ...
                    'bold','fontsize',fs,'Interpreter','latex') 
                    ylabel('\textbf{Mean $T_1$}','fontweight', ...
                            'bold','fontsize', fs -  25,'Interpreter','latex')
                    set(gcf, 'Position', get(0, 'Screensize'));
                    tempFigName = string(strcat(mainFolderName, modelsExtra, '/meanT1beta1', modelsExtra, figFormat));    
                    saveas(gcf,char(tempFigName),tempFigFormat) 

                    figure(201)
                    % each row is a model
                    for modelIndex = 1:numModels
                        t2BS = str2double(shortFileNames{modelIndex}(end-1));
                        if t2BS ~= 0
                            plot(betaRange2,tfFidMatT2(modelIndex,:),'-o', ...
                                    'LineWidth',4, 'Color',colors_T2(t2BS,:))
                            hold on 
                        end
                    end

                    ax = gca;
                    ax.FontSize = ticksFontSize; 
                    xlabel('$T_2$ production rate \textbf{$\gamma_1$}','fontweight', ...
                    'bold','fontsize',fs,'Interpreter','latex') 
                    ylabel('\textbf{Mean $T_2$}','fontweight', ...
                            'bold','fontsize', fs -  25,'Interpreter','latex')
                    set(gcf, 'Position', get(0, 'Screensize'));
                    tempFigName = string(strcat(mainFolderName, modelsExtra, '/meanT2gamma1', modelsExtra, figFormat));    
                    saveas(gcf,char(tempFigName),tempFigFormat) 

                    figure(202)
                    % normalized mean TF vs mean mRNA with all mean TF overlayed
                    for modelIndex = 1:numModels
                        normalizedHolderR = newFidMatT1(modelIndex,:)./max(newFidMatT1(modelIndex,:));
                        normalizedHolderTF = tfFidMatT1(modelIndex,:)./max(tfFidMatT1(modelIndex,:));
                        t1BS = str2double(shortFileNames{modelIndex}(end-2));
                        if t1BS ~= 0
                            plot(normalizedHolderTF, normalizedHolderR,'-o', ...
                                    'LineWidth',4, 'Color',colors_T1(t1BS,:))
                            hold on 
                            plot(normalizedHolderTF, normalizedHolderTF,'-x', ...
                                'MarkerSize',10, 'MarkerEdgeColor','black', 'MarkerFaceColor',[0 0 0], ...
                                'LineWidth',4, 'Color',colors_T1(t1BS,:))
                        end
                    end

                    ax = gca;
                    ax.FontSize = ticksFontSize; 
                    xlabel('Mean $T_1$ (Norm)','fontweight', ...
                    'bold','fontsize',fs,'Interpreter','latex') 
                    ylabel('\textbf{Mean mRNA} (Norm)','fontweight', ...
                            'bold','fontsize', fs -  25,'Interpreter','latex')
                    set(gcf, 'Position', get(0, 'Screensize'));
                    tempFigName = string(strcat(mainFolderName, modelsExtra, '/meanT1meanR',modelsExtra,  figFormat));    
                    saveas(gcf,char(tempFigName),tempFigFormat) 

                    figure(203)
                    % mean TF vs mean mRNA with all mean TF overlayed
                    for modelIndex = 1:numModels
                        normalizedHolderR = newFidMatT2(modelIndex,:)./max(newFidMatT2(modelIndex,:));
                        normalizedHolderTF = tfFidMatT2(modelIndex,:)./max(tfFidMatT2(modelIndex,:));
                        t2BS = str2double(shortFileNames{modelIndex}(end-2));
                        if t2BS ~= 0
                            plot(normalizedHolderTF, normalizedHolderR,'-o', ...
                                    'LineWidth', 4, 'Color',colors_T2(t2BS,:))
                            hold on 
                            plot(normalizedHolderTF, normalizedHolderTF,'-x', ...
                                'MarkerSize',10, 'MarkerEdgeColor','black', 'MarkerFaceColor',[0 0 0], ...
                                'LineWidth', 4, 'Color',colors_T2(t2BS,:))
                        end
                    end

                    ax = gca;
                    ax.FontSize = ticksFontSize; 
                    xlabel('Mean $T_2$ (Norm)','fontweight', ...
                    'bold','fontsize',fs,'Interpreter','latex') 
                    ylabel('\textbf{Mean mRNA} (Norm)','fontweight', ...
                            'bold','fontsize', fs -  25,'Interpreter','latex')
                    set(gcf, 'Position', get(0, 'Screensize'));
                    tempFigName = string(strcat(mainFolderName, modelsExtra, '/meanT2meanR',modelsExtra,  figFormat));    
                    saveas(gcf,char(tempFigName),tempFigFormat) 

                end

                if newFidelityBool
                    modelIndex = 4;
                    figure(10)
                    errorbar(betaRange1,newFidMatT1(modelIndex,:), ...
                        newFidStdMatT1(modelIndex,:),'-o','LineWidth',4, 'Color','k')    
                    ax = gca;
                    ax.FontSize = ticksFontSize; 
                    xlabel('$T_1$ production rate \textbf{$\beta_1$}','fontweight', ...
                    'bold','fontsize',fs,'Interpreter','latex') 
                    ylabel('\textbf{Mean mRNA}','fontweight', ...
                        'bold','fontsize', fs,'Interpreter','latex')
                    set(gcf, 'Position', get(0, 'Screensize'));
                    tempFigName = string(strcat(mainFolderName, modelsExtra, '/newfidBeta1ErrorBar', modelsExtra, figFormat));
                    saveas(gcf,char(tempFigName),tempFigFormat)


                    figure(11)
                    errorbar(betaRange2,newFidMatT2(modelIndex,:), ...
                        newFidStdMatT2(modelIndex,:),'-o','LineWidth',4, 'Color','k')   
                    ax = gca;
                    ax.FontSize = ticksFontSize; 
                    xlabel('$T_2$ production rate \textbf{$\gamma_1$}','fontweight', ...
                    'bold','fontsize',fs,'Interpreter','latex') 
                    ylabel('\textbf{Mean mRNA}','fontweight', ...
                        'bold','fontsize', fs,'Interpreter','latex')
                    set(gcf, 'Position', get(0, 'Screensize'));
                    tempFigName = string(strcat(mainFolderName, modelsExtra, '/newfidGamma1ErrorBar', modelsExtra, figFormat));
                    saveas(gcf,char(tempFigName),tempFigFormat)
                end


                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if efficiencyBool && ~newFidelityBool
                    figure(50)
                    for modelIndex = 1:numModels
                        t2BS = str2double(shortFileNames{modelIndex}(end-1));
                        if t2BS ~= 0
                            plot(betaRange1,efficiencyMatT2fixed(modelIndex,:),'-o', ...
                                'LineWidth',4,'Color',colors_T2(t2BS,:))
                            hold on
                        end 
                    end
                    ax = gca;
                    ax.FontSize = ticksFontSize; 
                    xlabel('\textbf{$\beta_1$}','fontweight', ...
                    'bold','fontsize',fs,'Interpreter','latex') 
                    ylabel('\textbf{Fidelity $T_2$}','fontweight', ...
                        'bold','fontsize', fs,'Interpreter','latex')
                    set(gcf, 'Position', get(0, 'Screensize'));
                    tempFigName = string(strcat(mainFolderName, modelsExtra, '/efficiencyFixedT2',modelsExtra,  figFormat));
                    saveas(gcf,char(tempFigName),tempFigFormat)


                    if topModelIndex < 4
                        figure(51)
                        for modelIndex = 1:numModels
                            t1BS = str2double(shortFileNames{modelIndex}(end-2));
                            if t1BS ~= 0
                                plot(betaRange2,efficiencyMatT1fixed(modelIndex,:),'-o', ...
                                    'LineWidth',4,'Color',colors_T1(t1BS,:))
                                hold on 
                            end
                        end
                        ax = gca;
                        ax.FontSize = ticksFontSize; 
                        xlabel('\textbf{$\gamma_1$}','fontweight', ...
                        'bold','fontsize',fs,'Interpreter','latex') 
                        ylabel('\textbf{Fidelity $T_1$}','fontweight', ...
                            'bold','fontsize', fs,'Interpreter','latex')
                        set(gcf, 'Position', get(0, 'Screensize'));
                        tempFigName = string(strcat(mainFolderName, modelsExtra, '/efficiencyFixedT1',modelsExtra,  figFormat));
                        saveas(gcf,char(tempFigName),tempFigFormat)        
                    end
                end
            end
        end
    end

% top model end
end

close all

toc     