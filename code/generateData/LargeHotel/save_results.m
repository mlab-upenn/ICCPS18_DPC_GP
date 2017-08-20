function save_results(fileExt, fileType, savePython, saveMatlab)

if saveMatlab
    RESULTS = extract_results(['Output/' fileType 'Building.csv']);
    save(['../../data/' fileExt '-LargeHotel'], '-struct', 'RESULTS');
    
end

if savePython
    RESULTS = extract_results(['Output/' fileType 'Building.csv']);
    n = size(RESULTS.DOW,1);
    p = numel(fieldnames(RESULTS));
    allData = zeros(n,p);
    varNames = fieldnames(RESULTS);
    for idp = 1:p
        allData(:,idp) = RESULTS.(varNames{idp});
    end
    csvwrite_with_headers(['../../data/' fileExt '-LargeHotel.csv'],allData,varNames);

end