function htest = holmbon(nPList, nAlpha)
%HOLMBO Holm-Bonferroni method to control probability of false rejections.
%
%   htest = holmbon(nPList, nAlpha)
%
%   When considering several hypotheses, the problem of multiplicity ari-
%   ses: the more hypotheses are checked, the higher the probability of 
%   obtaining Type I errors (false positives). The Holm–Bonferroni method
%   is one of many approaches for controlling the FWER, i.e., the probabi-
%   lity that one or more Type I errors will occur, by adjusting the reje-
%   ction criteria for each of the individual hypotheses [1].
%
%   1. https://en.wikipedia.org/wiki/Holm%E2%80%93Bonferroni_method
%
%   INPUT =================================================================
%
%   nPList (numeric array)
%   List of p-values. 
%   Example: [0.001, 0.01, 0.2]
%
%   nAlpha (numeric)
%   Significance level.
%   Example: 0.05
%
%   OUTPUT ================================================================
%
%   htest (structure)
%   Contains original p-values and the corresponding decision to either
%   reject null hypothesis (1) or confirm it (0). The null hypothesis 
%   being: there is no significant difference between specified popula-
%   tions, any observed difference being due to sampling or experimental
%   error.
%
%   AUTHOR ================================================================
%   
%   S.Bahdasariants, NEL, WVU, sb0220@mix.wvu.edu
%
%   See also MULTICOMPARE RANKSUM TTEST
%
%   =======================================================================

% Suppose you have m p-values, sorted into order lowest-to-highest, and 
% their corresponding hypotheses（null hypotheses). You want the FWER to
% be no higher than a certain pre-specified significance level alpha.

% Number of p-values.
nM = length(nPList); 

% Sort p-values from lowest to highest. 
nPList = sort(nPList); 

%% Holm-Bonferroni Procedure.

% Loop through p-values. 
for iElement = 1:nM

    % Fetch currently analyzed p-value. 
    nP = nPList(iElement); 

    % Assign the p-value to the output structure. 
    htest.pvalue(iElement) = nP; 

    % Assume no significant difference. 
    htest.bHypothesis(iElement) = 0;

    % Apply Holm-Bonferroni correction. 
    if nP < (nAlpha/(nM + 1 - iElement))

        % Reject null hypothesis. The difference was found. 
        htest.bHypothesis(iElement) = 1;

    end % if

end % iElement

end % function