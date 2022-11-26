function nCircum = setCircum(nAge, sSex, sBody, sFile)
    %SETCIRCUM Fetches the circumference of the sBody from the sFile.
    %
    %   nCircum = setCircum(nAge, sSex, sBody, sFile)
    %
    %   INPUT =============================================================
    %
    %   nAge (double)
    %   Age of the subject. 
    %   Example: 23
    %
    %   sSex (string)
    %   Subject's sex. 
    %   Example: 'M' (options: "M", "F")
    %
    %   sBody (string)
    %   The name of the body of iterest. 
    %   Example: 'calf' (options: "waist", "thigh", "calf")
    %   
    %   sFile (string)
    %   Name of the xlsx file in which the circumference data are stored. 
    %   Example: "metaCircumference"
    %
    %   OUTPUT ============================================================
    %
    %   nCircum (double)
    %   Circumference of the body of interest chosen with respect to 
    %       subject's age and sex. 
    %
    %   AUTHOR ============================================================
    %   
    %   S.Bahdasariants, NEL, WVU, sb0220@mix.wvu.edu
    %
    %   See also MAIN, SOLVEDYNAMICS, SIMSWING, EXTRACTMETA, SAVESIM,...
    %   GETSWING, RUNSIM, GETKIN, SCALEANTHRO, SETMECHANICS, PARAXT,...
    %   FRUSTUMINERT
    %
    %   LITERATURE ========================================================
    %
    %   1. McDowell, M. A., Fryar, C. D., Hirsch, R. & Ogden, C. L. 
    %   Anthropometric reference data for children and adults: U.S. 
    %   population, 1999-2002. Adv Data 1–5 (2005).
    %
    %   2. Fryar, C. D., Carroll, M. D., Gu, Q., Afful, J. & Ogden,
    %   C. L. Anthropometric reference data for children and adults: 
    %   United States, 2015-2018. (2021).
    %
    %   ===================================================================
    
    %% LOAD DATA. Load the file with the circumerences specified for 
    %   different bodies, ages, and sexes. 
    
    metaCircum = readtable(sFile); % load the xlsx file
    
    %% FIND CIRCUMFERENCE. Find the row number that corresponds to the age,
    %   sex, and the name of the body of interest. Next, assing respective 
    %   circumference to the output variable. 
    
    % The body of interest is thigh
    if sBody == "thigh"     
        bThigh = 1;
        bCalf  = 0;
        bWaist = 0;

    % The body of interest is shank
    elseif sBody == "calf"  
        bThigh = 0;
        bCalf  = 1;
        bWaist = 0;
    
    % The body of interest is waist
    elseif sBody == "waist" 
        bThigh = 0;
        bCalf  = 0;
        bWaist = 1;
    end
    
    % Find the row index 
    idxCircum = (string(metaCircum.sSex) == sSex &...              % Sex
        metaCircum.nAgeMin<= nAge & metaCircum.nAgeMax>= nAge &... % Age
        metaCircum.bWaist == bWaist &...                           % Body
        metaCircum.bCalf  == bCalf &...
        metaCircum.bThigh == bThigh);
    
    % Fetch the circumference 
    nCircum = metaCircum.nCircumference(find(idxCircum)); 
    
end % function

