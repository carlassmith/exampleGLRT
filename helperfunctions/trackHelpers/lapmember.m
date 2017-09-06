function varargout = lapmember(A,B,flag1,flag2)
%ISMEMBER True for set member.
%   LIA = ISMEMBER(A,B) for arrays A and B returns an array of the same
%   size as A containing true where the elements of A are in B and false
%   otherwise.
%
%   LIA = ISMEMBER(A,B,'rows') for matrices A and B with the same number
%   of columns, returns a vector containing true where the rows of A are
%   also rows of B and false otherwise.
%
%   [LIA,LOCB] = ISMEMBER(A,B) also returns an array LOCB containing the
%   lowest absolute index in B for each element in A which is a member of
%   B and 0 if there is no such index.
%
%   [LIA,LOCB] = ISMEMBER(A,B,'rows') also returns a vector LOCB containing
%   the lowest absolute index in B for each row in A which is a member
%   of B and 0 if there is no such index.
%
%   The behavior of ISMEMBER has changed.  This includes:
%     -	occurrence of indices in LOCB switched from highest to lowest
%     -	tighter restrictions on combinations of classes
%
%   If this change in behavior has adversely affected your code, you may 
%   preserve the previous behavior with:
%
%      [LIA,LOCB] = ISMEMBER(A,B,'legacy')
%      [LIA,LOCB] = ISMEMBER(A,B,'rows','legacy')
%
%   Examples:
%
%      a = [9 9 8 8 7 7 7 6 6 6 5 5 4 4 2 1 1 1]
%      b = [1 1 1 3 3 3 3 3 4 4 4 4 4 9 9 9]
%
%      [lia1,locb1] = ismember(a,b)
%      % returns
%      lia1 = [1 1 0 0 0 0 0 0 0 0 0 0 1 1 0 1 1 1]
%      locb1 = [14 14 0 0 0 0 0 0 0 0 0 0 9 9 0 1 1 1]
%
%      [lia,locb] = ismember([1 NaN 2 3],[3 4 NaN 1])
%      % NaNs compare as not equal, so this returns
%      lia = [1 0 0 1], locb = [4 0 0 1]
%
%   Class support for inputs A and B, where A and B must be of the same
%   class unless stated otherwise:
%      - logical, char, all numeric classes (may combine with double arrays)
%      - cell arrays of strings (may combine with char arrays)
%      -- 'rows' option is not supported for cell arrays
%      - objects with methods SORT (SORTROWS for the 'rows' option), EQ and NE
%      -- including heterogeneous arrays derived from the same root class
%
%   See also UNIQUE, UNION, INTERSECT, SETDIFF, SETXOR, SORT, SORTROWS.

%   Copyright 1984-2013 The MathWorks, Inc.

if nargin == 2
    [varargout{1:max(1,nargout)}] = ismemberR2012a(A,B);
else
    % acceptable combinations, with optional inputs denoted in []
    % ismember(A,B, ['rows'], ['legacy'/'R2012a'])
    nflagvals = 3;
    flagvals = {'rows' 'legacy' 'R2012a'};
    % When a flag is found, note the index into varargin where it was found
    flaginds = zeros(1,nflagvals);
    for i = 3:nargin
        if i == 3
            flag = flag1;
        else
            flag = flag2;
        end
        foundflag = strcmpi(flag,flagvals);
        if ~any(foundflag)
            if ischar(flag)
                error(message('MATLAB:ISMEMBER:UnknownFlag',flag));
            else
                error(message('MATLAB:ISMEMBER:UnknownInput'));
            end
        end
        % Only 1 occurrence of each allowed flag value
        if flaginds(foundflag)
            error(message('MATLAB:ISMEMBER:RepeatedFlag',flag));
        end
        flaginds(foundflag) = i;
    end
        
    % Only 1 of each of the paired flags
    if flaginds(2) && flaginds(3)
        error(message('MATLAB:ISMEMBER:BehaviorConflict'))
    end
    % 'legacy' and 'R2012a' flags must be trailing
    if flaginds(2) && flaginds(2)~=nargin
        error(message('MATLAB:ISMEMBER:LegacyTrailing'))
    end
    if flaginds(3) && flaginds(3)~=nargin
        error(message('MATLAB:ISMEMBER:R2012aTrailing'))
    end
    
    if flaginds(3) % trailing 'R2012a' specified
        [varargout{1:max(1,nargout)}] = ismemberR2012a(A,B,logical(flaginds(1)));
    elseif flaginds(2) % trailing 'legacy' specified
        [varargout{1:max(1,nargout)}] = ismemberlegacy(A,B,logical(flaginds(1)));
    else % 'R2012a' (default behavior)
        [varargout{1:max(1,nargout)}] = ismemberR2012a(A,B,logical(flaginds(1)));
    end
end
end

function [tf,loc] = ismemberlegacy(a,s,isrows)
% 'legacy' flag implementation

if nargin == 3 && isrows
    flag = 'rows';
else
    flag = [];
end

numelA = numel(a);
numelS = numel(s);
nOut = nargout;

if ~(isa(a,'opaque') || isa(s,'opaque'))
    
    if isempty(flag)
        
        % Initialize types and sizes.
        
        tf = false(size(a));
        
        if nOut > 1
            loc = zeros(size(a));
        end
        
        % Handle empty arrays and scalars.
        
        if numelA == 0 || numelS <= 1
            if (numelA == 0 || numelS == 0)
                return
                % Scalar A handled below.
                % Scalar S: find which elements of A are equal to S.
            elseif numelS == 1
                tf = (a == s);
                if nOut > 1
                    % Use DOUBLE to convert logical "1" index to double "1" index.
                    loc = double(tf);
                end
                return
            end
        else
            % General handling.
            % Use FIND method for very small sizes of the input vector to avoid SORT.
            scalarcut = 5;
            if numelA <= scalarcut
                if nOut <= 1
                    for i=1:numelA
                        tf(i) = any(a(i)==s(:));   % ANY returns logical.
                    end
                else
                    for i=1:numelA
                        found = find(a(i)==s(:));  % FIND returns indices for LOC.
                        if ~isempty(found)
                            tf(i) = 1;
                            loc(i) = found(end);
                        end
                    end
                end
            else
                % Use method which sorts list, then performs binary search.
                % Convert to double for quicker sorting, to full to work in C helper.
                a = double(a);
                if issparse(a)
                    a = full(a);
                end
                
                s = double(s);
                if issparse(s)
                    s = full(s);
                end
                
                if (isreal(s))
                    % Find out whether list is presorted before sort
                    % If the list is short enough, SORT will be faster than ISSORTED
                    % If the list is longer, ISSORTED can potentially save time
                    checksortcut = 1000;
                    if numelS > checksortcut
                        sortedlist = issorted(s(:));
                    else
                        sortedlist = 0;
                    end
                    if nOut > 1
                        if ~sortedlist
                            [s,idx] = sort(s(:));
                        end
                    elseif ~sortedlist
                        s = sort(s(:));
                    end
                else
                    sortedlist = 0;
                    [~,idx] = sort(real(s(:)));
                    s = s(idx);
                end
                
                % Two C-Helper Functions are used in the code below:
                
                % ISMEMBC  - S must be sorted - Returns logical vector indicating which
                % elements of A occur in S
                % ISMEMBC2 - S must be sorted - Returns a vector of the locations of
                % the elements of A occurring in S.  If multiple instances occur,
                % the last occurrence is returned
                
                % Check for NaN values - NaN values will be at the end of S,
                % but may be anywhere in A.
                
                nana = isnan(a(:));
                
                if (any(nana) || isnan(s(numelS)))
                    % If NaNs detected, remove NaNs from the data before calling ISMEMBC.
                    ida = (nana == 0);
                    ids = (isnan(s(:)) == 0);
                    if nOut <= 1
                        ainfn = ismembc(a(ida),s(ids));
                        tf(ida) = ainfn;
                    else
                        loc1 = ismembc2(a(ida),s(ids));
                        tf(ida) = (loc1 > 0);
                        loc(ida) = loc1;
                        loc(~ida) = 0;
                    end
                else
                    % No NaN values, call ISMEMBC directly.
                    if nOut <= 1
                        tf = ismembc(a,s);
                    else
                        loc = ismembc2(a,s);
                        tf = (loc > 0);
                    end
                end
                
                if nOut > 1 && ~sortedlist
                    % Re-reference loc to original list if it was unsorted
                    loc(tf) = idx(loc(tf));
                end
            end
        end
        
    else    % 'rows' case
        
        rowsA = size(a,1);
        colsA = size(a,2);
        rowsS = size(s,1);
        colsS = size(s,2);
        
        % Automatically pad strings with spaces
        if ischar(a) && ischar(s),
            if colsA > colsS
                s = [s repmat(' ',rowsS,colsA-colsS)];
            elseif colsA < colsS
                a = [a repmat(' ',rowsA,colsS-colsA)];
            end
        elseif colsA ~= colsS && ~isempty(a) && ~isempty(s)
            error(message('MATLAB:ISMEMBER:AandBColnumAgree'));
        end
        
        % Empty check for 'rows'.
        if rowsA == 0 || rowsS == 0
            if (isempty(a) || isempty(s))
                tf = false(rowsA,1);
                loc = zeros(rowsA,1);
                return
            end
        end
        
        % General handling for 'rows'.
        
        % Duplicates within the sets are eliminated
        if (rowsA == 1)
            au = repmat(a,rowsS,1);
            d = au(1:end,:)==s(1:end,:);
            d = all(d,2);
            tf = any(d);
            if nOut > 1
                if tf
                    loc = find(d, 1, 'last');
                else
                    loc = 0;
                end
            end
            return;
        else
            [au,~,an] = unique(a,'rows','legacy');
        end
        if nOut <= 1
            su = unique(s,'rows','legacy');
        else
            [su,sm] = unique(s,'rows','legacy');
        end
        
        % Sort the unique elements of A and S, duplicate entries are adjacent
        [c,ndx] = sortrows([au;su]);
        
        % Find matching entries
        d = c(1:end-1,:)==c(2:end,:);     % d indicates matching entries in 2-D
        d = find(all(d,2));               % Finds the index of matching entries
        ndx1 = ndx(d);                    % NDX1 are locations of repeats in C
        
        if nOut <= 1
            tf = ismember(an,ndx1,'legacy');         % Find repeats among original list
        else
            szau = size(au,1);
            [tf,loc] = ismember(an,ndx1,'legacy');   % Find loc by using given indices
            newd = d(loc(tf));              % NEWD is D for non-unique A
            where = sm(ndx(newd+1)-szau);  % Index values of SU through UNIQUE
            loc(tf) = where;                % Return last occurrence of A within S
        end
    end
else
    % Handle objects that cannot be converted to doubles
    if isempty(flag)
        
        % Handle empty arrays and scalars.
        
        if numelA == 0 || numelS <= 1
            if (numelA == 0 || numelS == 0)
                tf = false(size(a));
                loc = zeros(size(a));
                return
                
                % Scalar A handled below.
                % Scalar S: find which elements of A are equal to S.
            elseif numelS == 1
                tf = (a == s);
                if nOut > 1
                    % Use DOUBLE to convert logical "1" index to double "1" index.
                    loc = double(tf);
                end
                return
            end
        else
            % General handling.
            % Use FIND method for very small sizes of the input vector to avoid SORT.
            scalarcut = 5;
            if numelA <= scalarcut
                tf = false(size(a));
                loc = zeros(size(a));
                if nOut <= 1
                    for i=1:numelA
                        tf(i) = any(a(i)==s);   % ANY returns logical.
                    end
                else
                    for i=1:numelA
                        found = find(a(i)==s);  % FIND returns indices for LOC.
                        if ~isempty(found)
                            tf(i) = 1;
                            loc(i) = found(end);
                        end
                    end
                end
            else
                
                % Duplicates within the sets are eliminated
                [au,~,an] = unique(a(:),'legacy');
                if nOut <= 1
                    su = unique(s(:),'legacy');
                else
                    [su,sm] = unique(s(:),'legacy');
                end
                
                % Sort the unique elements of A and S, duplicate entries are adjacent
                [c,ndx] = sort([au;su]);
                
                % Find matching entries
                d = c(1:end-1)==c(2:end);         % d indicates matching entries in 2-D
                d = find(d);                      % Finds the index of matching entries
                ndx1 = ndx(d);                    % NDX1 are locations of repeats in C
                
                if nOut <= 1
                    tf = ismember(an,ndx1,'legacy');         % Find repeats among original list
                else
                    szau = size(au,1);
                    [tf,loc] = ismember(an,ndx1,'legacy');   % Find loc by using given indices
                    newd = d(loc(tf));              % NEWD is D for non-unique A
                    where = sm(ndx(newd+1)-szau);   % Index values of SU through UNIQUE
                    loc(tf) = where;                % Return last occurrence of A within S
                end
            end
            tf = reshape(tf,size(a));
            if nOut > 1
                loc = reshape(loc,size(a));
            end
        end
        
    else    % 'rows' case
        
        rowsA = size(a,1);
        colsA = size(a,2);
        rowsS = size(s,1);
        colsS = size(s,2);
        
        % Automatically pad strings with spaces
        if ischar(a) && ischar(s),
            if colsA > colsS
                s = [s repmat(' ',rowsS,colsA-colsS)];
            elseif colsA < colsS
                a = [a repmat(' ',rowsA,colsS-colsA)];
            end
        elseif size(a,2)~=size(s,2) && ~isempty(a) && ~isempty(s)
            error(message('MATLAB:ISMEMBER:AandBColnumAgree'));
        end
        
        % Empty check for 'rows'.
        if rowsA == 0 || rowsS == 0
            if (isempty(a) || isempty(s))
                tf = false(rowsA,1);
                loc = zeros(rowsA,1);
                return
            end
        end
        
        % Duplicates within the sets are eliminated
        [au,~,an] = unique(a,'rows','legacy');
        if nOut <= 1
            su = unique(s,'rows','legacy');
        else
            [su,sm] = unique(s,'rows','legacy');
        end
        
        % Sort the unique elements of A and S, duplicate entries are adjacent
        [c,ndx] = sortrows([au;su]);
        
        % Find matching entries
        d = c(1:end-1,:)==c(2:end,:);     % d indicates matching entries in 2-D
        d = find(all(d,2));               % Finds the index of matching entries
        ndx1 = ndx(d);                    % NDX1 are locations of repeats in C
        
        if nOut <= 1
            tf = ismember(an,ndx1,'legacy');         % Find repeats among original list
        else
            szau = size(au,1);
            [tf,loc] = ismember(an,ndx1,'legacy');   % Find loc by using given indices
            newd = d(loc(tf));              % NEWD is D for non-unique A
            where = sm(ndx(newd+1)-szau);   % Index values of SU through UNIQUE
            loc(tf) = where;                % Return last occurrence of A within S
        end
    end
end
end

function [lia,locb] = ismemberR2012a(a,b,options)
% 'R2012a' flag implementation

% Error check flag
if nargin == 2
    byrow = false;
else
    byrow = options > 0;
end

classFlag = true;
% Check that one of A and B is double if A and B are non-homogeneous. Do a
% separate check if A is a heterogeneous object and only allow a B
% that is of the same root class.
if ~(isa(a,'handle.handle') || isa(b,'handle.handle'))
    if ~strcmpi(class(a),class(b))
        if isa(a,'matlab.mixin.Heterogeneous') && isa(b,'matlab.mixin.Heterogeneous')
            rootClassA = meta.internal.findHeterogeneousRootClass(a);
            if isempty(rootClassA) || ~isa(b,rootClassA.Name)
                error(message('MATLAB:ISMEMBER:InvalidInputsDataType',class(a),class(b)));
            end
        elseif ~(strcmpi(class(a),'double') || strcmpi(class(b),'double'))
            error(message('MATLAB:ISMEMBER:InvalidInputsDataType',class(a),class(b)));
        end
        classFlag = false;
    end
end

numelA = numel(a);

if ~byrow
    if ~(isa(a,'opaque') || isa(b,'opaque')) && (classFlag || numelA <= 5)     
        % Initialize types and sizes.
        if nargout > 1
            [lia,locb] = ismemberBuiltinTypes(a,b);
        else
            lia = ismemberBuiltinTypes(a,b);
        end

    else %Handle objects and allowed mixed data types
        if nargout > 1
            [lia,locb] = ismemberClassTypes(a,b);
        else
            lia = ismemberClassTypes(a,b);
        end
    end
else    % 'rows' case
    if ~(ismatrix(a) && ismatrix(b))
        error(message('MATLAB:ISMEMBER:NotAMatrix'));
    end
    
    [rowsA,colsA] = size(a);
    [rowsB,colsB] = size(b);
    
    % Automatically pad strings with spaces
    if ischar(a) && ischar(b),
        b = [b repmat(' ',rowsB,colsA-colsB)];
        a = [a repmat(' ',rowsA,colsB-colsA)];
    elseif colsA ~= colsB
        error(message('MATLAB:ISMEMBER:AandBColnumAgree'));
    end
    
    % Empty check for 'rows'.
    if rowsA == 0 || rowsB == 0
        lia = false(rowsA,1);
        locb = zeros(rowsA,1);
        return
    end
    
    % General handling for 'rows'.
    
    
    % overhead removed by PKR, data is already sorted, yeah!
%     % Duplicates within the sets are eliminated
%     if (rowsA == 1)
%         uA = repmat(a,rowsB,1);
%         d = uA(1:end,:)==b(1:end,:);
%         d = all(d,2);
%         lia = any(d);
%         if nargout > 1
%             if lia
%                 locb = find(d, 1, 'first');
%             else
%                 locb = 0;
%             end
%         end
%         return;
%     else
%         [uA,~,icA] = unique(a,'rows','sorted');
%     end
%     if nargout <= 1
%         uB = unique(b,'rows','sorted');
%     else
%         [uB,ib] = unique(b,'rows','sorted');
%     end
    
    % Sort the unique elements of A and B, duplicate entries are adjacent
    %[sortuAuB,IndSortuAuB] = sortrows([uA;uB]);
    
    % modified by PKR, don't need to carry over dummy variables, no changes
    % here.
    [sortuAuB,IndSortuAuB] = sortrows([a;b]);
    icA = (1:size(a,1))';
    
    % Find matching entries
    d = sortuAuB(1:end-1,:)==sortuAuB(2:end,:);     % d indicates matching entries
    d = all(d,2);                                   % Finds the index of matching entries
    ndx1 = IndSortuAuB(d);                          % NDX1 are locations of repeats in C
    
    if nargout <= 1
        lia = ismemberBuiltinTypes(icA,ndx1);           % Find repeats among original list
    else
        szuA = size(uA,1);
        [lia,locb] = ismemberBuiltinTypes(icA,ndx1);    % Find locb by using given indices
        d = find(d);
        newd = d(locb(lia));                    % NEWD is D for non-unique A
        where = ib(IndSortuAuB(newd+1)-szuA);   % Index values of uB through UNIQUE
        locb(lia) = where;                      % Return first or last occurrence of A within B
    end
end
end

function [lia,locb] = ismemberBuiltinTypes(a,b)
% General handling.
% Use FIND method for very small sizes of the input vector to avoid SORT.
lia = false(size(a));
if nargout > 1
    locb = zeros(size(a));
end
% Handle empty arrays and scalars.  
numelA = numel(a);
numelB = numel(b);
if numelA == 0 || numelB <= 1
    if (numelA == 0 || numelB == 0)
        return
        % Scalar A handled below.
        % Scalar B: find which elements of A are equal to B.
    elseif numelB == 1
        lia = (a == b);
        if nargout > 1
            % Use DOUBLE to convert logical "1" index to double "1" index.
            locb = double(lia);
        end
        return
    end
end

lia = false(size(a));
if nargout > 1
    locb = zeros(size(a));
end
scalarcut = 5;
numelA = numel(a);
numelB = numel(b);
if numelA <= scalarcut
    if nargout <= 1
        for i=1:numelA
            lia(i) = any(a(i)==b(:));   % ANY returns logical.
        end
    else
        for i=1:numelA
            found = a(i)==b(:);  % FIND returns indices for LOCB.
            if any(found)
                lia(i) = true;
                found = find(found);
                locb(i) = found(1);
            end
        end
    end
else
    % Use method which sorts list, then performs binary search.
    % Convert to full to work in C helper.
    if issparse(a)
        a = full(a);
    end
    if issparse(b)
        b = full(b);
    end
    
    if (isreal(b))
        % Find out whether list is presorted before sort
        % If the list is short enough, SORT will be faster than ISSORTED
        % If the list is longer, ISSORTED can potentially save time
        checksortcut = 1000;
        if numelB > checksortcut
            sortedlist = issorted(b(:));
        else
            sortedlist = 0;
        end
        if nargout > 1
            if ~sortedlist
                [b,idx] = sort(b(:));
            end
        elseif ~sortedlist
            b = sort(b(:));
        end
    else
        sortedlist = 0;
        [~,idx] = sort(real(b(:)));
        b = b(idx);
    end
    
    % C++ - Helper Function are used in the code below:
    
    % ISMEMBERONEOUTPUT(A,B)  - B must be sorted - Returns logical vector indicating which
    % elements of A occur in B
    % ISMEMBERLAST(A,B) - B must be sorted - Returns a vector of the locations of
    % the elements of A occurring in B.  If multiple instances occur,
    % the last occurrence is returned (Not currently being used)
    % ISMEMBERFIRST(A,B) - B must be sorted - Returns a vector of the
    % locations of the elements of A occurring in B.  If multiple
    % instances occur, the first occurence is returned.
    
    % Check for NaN values - NaN values will be at the end of B,
    % but may be anywhere in A.
    
    nana = isnan(a(:));
    if builtin('isnumeric',a) && builtin('isnumeric',b) && strcmp(class(a),class(b))
        if (any(nana) || isnan(b(numelB)))
            % If NaNs detected, remove NaNs from the data
            ida = ~nana;
            idb = ~isnan(b(:));
            if nargout <= 1
                ainfn = builtin('_ismemberoneoutput',a(ida),b(idb));
                lia(ida) = ainfn;
            else
                loc1 = builtin('_ismemberfirst',a(ida),b(idb));
                lia(ida) = (loc1 > 0);
                locb(ida) = loc1;
                locb(~ida) = 0;
            end
        else 
            if nargout <= 1
                lia = builtin('_ismemberoneoutput',a,b);
            else
                locb = builtin('_ismemberfirst',a,b);
                lia = (locb > 0);
            end
        end
    else %(a,b, are some other class like gpuArray, syb object)
        if nargout <= 1
            for i=1:numelA
                lia(i) = any(a(i)==b(:));   % ANY returns logical.
            end
        else
            for i=1:numelA
                found = a(i)==b(:); % FIND returns indices for LOCB.
                if any(found)
                    lia(i) = true;
                    found = find(found);
                    locb(i) = found(1);
                end
            end
        end
    end
    if nargout > 1 && ~sortedlist
        % Re-reference locb to original list if it was unsorted
        locb(lia) = idx(locb(lia));
    end
end
end
function [lia,locb] = ismemberClassTypes(a,b)
if issparse(a)
    a = full(a);
end
if issparse(b)
    b = full(b);
end
% Duplicates within the sets are eliminated
if isscalar(a) || isscalar(b)
    ab = [a(:);b(:)];
    sort(ab(1));
    numa = numel(a);
    lia = ab(1:numa)==ab(1+numa:end);
    if ~any(lia)
        lia  = false(size(a));
        locb = zeros(size(a));
        return
    end
    if ~isscalar(b)
        locb = find(lia);
        locb = locb(1);
        lia = any(lia);
    else
        locb = double(lia);
    end
else
    % Duplicates within the sets are eliminated
    [uA,~,icA] = unique(a(:),'sorted');
    if nargout <= 1
        uB = unique(b(:),'sorted');
    else
        [uB,ib] = unique(b(:),'sorted');
    end
    
    % Sort the unique elements of A and B, duplicate entries are adjacent
    [sortuAuB,IndSortuAuB] = sort([uA;uB]);
    
    % Find matching entries
    d = sortuAuB(1:end-1)==sortuAuB(2:end);         % d indicates the indices matching entries
    ndx1 = IndSortuAuB(d);                          % NDX1 are locations of repeats in C
    
    if nargout <= 1
        lia = ismemberBuiltinTypes(icA,ndx1);       % Find repeats among original list
    else
        szuA = size(uA,1);
        d = find(d);
        [lia,locb] = ismemberBuiltinTypes(icA,ndx1);% Find locb by using given indices
        newd = d(locb(lia));                        % NEWD is D for non-unique A
        where = ib(IndSortuAuB(newd+1)-szuA);       % Index values of uB through UNIQUE
        locb(lia) = where;                          % Return first or last occurrence of A within B
    end
end
lia = reshape(lia,size(a));
if nargout > 1
    locb = reshape(locb,size(a));
end
end