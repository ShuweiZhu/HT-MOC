function DataName=InputData(DD)

switch DD
    %% Sythentic data sets  %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1.Small data sets
    case 1
        DataName='Flame';
    case 2
        DataName='Aggregation';
    case 3
        DataName='D31';
    case 4
        DataName = 'T4.8k';
    case 5
        DataName = 'T7.1k';        
    % 2.Large data sets
    case 6
        DataName = 'TB20K';
    case 7
        DataName = 'SF50K';
    case 8
        DataName = 'CC100K';
    case 9
        DataName = 'CG200K';        
    % 3.Huge data sets
    case 10
        DataName = 'GaussianL1';
    case 11
        DataName = 'GaussianL2';
    case 12
        DataName = 'GaussianL3';
        
    %% Real data sets  %%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1.Small data sets
    case 16
        DataName='Cancer';
    case 17
        DataName='Wine';
    case 18
        DataName='Control';
    case 19
        DataName='Opticaldigit';
    case 20
        DataName ='PenDigits';              
    % 2.Large data sets
    case 21
        DataName = 'USPS';
    case 22
        DataName = 'ISOLET';
    case 23
        DataName='C_Cube';
    case 24
        DataName='Magic';        
    case 25
        DataName = 'Shuttle';
    case 26
        DataName = 'Covertype';
        
end

