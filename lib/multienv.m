function varargout = multienv(tensorList,envList,legLinks,sequence)
% multienv v1.02 (C) R. N. C. Pfeifer 2013-2015
% ==============
% function varargout = multienv(tensorList,envList,legLinks,sequence)
% Performs simultaneous computation of multiple environments with efficient re-use of intermediate results.

% v1.01:
% Fixed implicit outer product of more than two subnetworks.
% Improved error message when performing an outer product using zero notation and tensors are connected by a non-trivial index.
% Fixed handling of trailing traces after an outer product on final two tensors.
% Added extra error-checking to interpretation of zeros-in-sequence notation.

% v1.02: (Beyond published version)
% If number of output arguments requested exceeds the largest environment number requested, return empty cells for the extra arguments.
% Multienv unofficially supports replacing one tensor with a dummy entry ([]) if computing only the environment of this tensor. 
%   However, checkinputs would sometimes incorrectly parse the dimensions of legs connecting to this tensor if it was not the final 
%   tensor in the list. This has been fixed.

    % Sanity-check inputs
    % ===================
    if ~exist('sequence','var')
        % Sequence not specified - use default: Positive indices in ascending order
        [altposindices,legLinks,legDims,dummyTensorPos,orderLegs] = checkinputs(tensorList,envList,legLinks);
        sequence = altposindices;
    else
        [altposindices,legLinks,legDims,dummyTensorPos,orderLegs] = checkinputs(tensorList,envList,legLinks,sequence);
    end
    
    % Perform any initial traces
    % ==========================
    traced = [];
    for a=1:numel(tensorList)
        % Determine which sequence entries correspond to traces, tensor by tensor
        tracinglegs = sort(legLinks{a});
        tracinglegs = tracinglegs(find(tracinglegs(1:end-1)==tracinglegs(2:end))); %#ok<FNDSB>
        if ~isempty(tracinglegs)
            % Perform these traces, and remove from sequence
            tensorList{a} = dotrace(tensorList{a},legLinks{a},tracinglegs);
            traced = [traced tracinglegs]; %#ok<AGROW>
            for b=1:numel(tracinglegs)
                legLinks{a}(legLinks{a}==tracinglegs(b)) = [];
            end
        end
    end
    if ~isempty(traced)
        if ~isequal(sort(traced),sort(sequence(1:numel(traced))))
            warning('multienv:singletensortracemisplaced','Found one or more single-tensor traces not appearing at the beginning of the contraction sequence -- performed them first anyway.')
        end
    end
    for a=1:numel(traced)
        sequence(sequence==traced(a)) = [];
    end
    
    
    % Construct contraction data tree
    % ===============================
    [datatree,finaltraces] = constructtree(legLinks,sequence,envList,legDims,dummyTensorPos,orderLegs);
    
    % Warn if doing some final traces early
    % =====================================
    dummyLegs = logical(sparse([],[],[],1,altposindices(end),numel(altposindices)));
    if dummyTensorPos~=0
        dummyLegs(legLinks{dummyTensorPos}) = 1;
    end
    if any((legDims(orderLegs(finaltraces))>1) & ~full(logical(dummyLegs(finaltraces))))
        warning('multienv:finaltraces','Contraction sequence ends with traces over non-trivial indices. This is always suboptimal, and contractions over these indices will instead be performed when the tensors bearing these indices are contracted together.')
    end
    
    % Evaluate environments
    % =====================
    varargout = cell(1,max([max(envList) nargout]));
    if numel(envList)==1
        varargout{1} = 1;
    else
        for a=reshape(find(envList~=0),1,[])
            [newEnv,datatree] = calculateenvironment(tensorList,legLinks,datatree,a);
            if isempty(varargout{envList(a)})
                varargout{envList(a)} = newEnv;
            else
                varargout{envList(a)} = varargout{envList(a)} + newEnv;
            end
        end
    end
end

function B = dotrace(A,leglabels,tindices) % Trace over all indices listed in tindices, each of which should occur twice on tensor A
    sz = size(A);
    sz = [sz ones(1,numel(leglabels)-numel(sz))];
    tpos = [];
    % Find positions of tracing indices
    for a=1:numel(tindices)
        tpos = [tpos find(leglabels==tindices(a))]; %#ok<AGROW>
    end
    % Reorder list of tracing indices so that they occur in two equivalent blocks
    tpos = [tpos(1:2:end) tpos(2:2:end)];
    % Identify non-tracing indices
    ind = 1:numel(leglabels);
    ind(tpos) = [];
    % Collect non-tracing and tracing indices
    A = reshape(permute(A,[ind tpos]),prod(sz(ind)),prod(sz(tpos))); % Separate indices to be traced and not to be traced
    B = 0;
    % Perform trace
    szcount = sqrt(prod(sz(tpos)));
    for a=1:szcount
        B = B + A(:,a+(a-1)*szcount); % Perform trace
    end
    B = reshape(B,[sz(ind) 1 1]);
end

function [altposindices,legLinks,sizes,dummytensorpos,orderLegs] = checkinputs(tensorList,envList,legLinks,sequence)
    % Check data sizes
    if size(tensorList,1)~=1 || size(tensorList,2)~=numel(tensorList)
        error('Array of tensors has incorrect dimension - should be 1xn')
    end
    for a=1:numel(legLinks)
        if size(legLinks{a},1)~=1 || size(legLinks{a},2)~=numel(legLinks{a})
            if isempty(legLinks{a})
                legLinks{a} = zeros(1,0);
            else
                error(['Leg link entry ' num2str(a) ' has wrong dimension - should be 1xn']);
            end
        end
    end
    if ~isequal(size(legLinks),size(tensorList))
        error('Array of links should be the same size as the array of tensors')
    end
    % Check all tensors are numeric
    for a=1:numel(tensorList)
        if ~isnumeric(tensorList{a})
            error('Tensor list must be a 1xn cell array of numerical objects')
        end
    end
    % Is there a dummy tensor? One tensor may be [] provided the only environment to be calculated belongs to that tensor.
    numdummytensors = 0;
    dummytensorpos = 0;
    for a=1:numel(tensorList)
        if isequal(tensorList{a},[])
            numdummytensors = numdummytensors+1;
            dummytensorpos = a;
        end
    end
    if numdummytensors > 1
        error('Only one dummy tensor is permitted')
    end
    if numdummytensors == 1
        envListCopy = envList;
        envListCopy(dummytensorpos) = [];
        if any(envListCopy~=0)
            error('When using a dummy tensor, only the environment of that tensor may be computed')
        end
    end
    % Get list of positive indices
    allindices = cell2mat(legLinks);
    if any(allindices==0)
        error('Zero entry in index list')
    elseif any(imag(allindices)~=0)
        error('Complex entry in index list')
    elseif any(allindices<0)
        error('Negative entry in index list')
    end
    [posindices,ix] = sort(allindices,'ascend');
    % Test all positive indices occur exactly twice
    if mod(numel(posindices),2)~=0
        maxposindex = posindices(end);
        posindices = posindices(1:end-1);
    end
    flags = (posindices(1:2:numel(posindices))-posindices(2:2:numel(posindices)))~=0;
    if any(flags)
        errorpos = 2*find(flags~=0,1,'first')-1;
        if errorpos>1 && posindices(errorpos-1)==posindices(errorpos)
            error(['Error in index list: Index ' num2str(posindices(errorpos)) ' appears more than twice']);
        else
            error(['Error in index list: Index ' num2str(posindices(errorpos)) ' only appears once']);
        end
    end
    if exist('maxposindex','var')
        if isempty(posindices)
            error(['Error in index list: Index ' num2str(maxposindex) ' only appears once']);
        end
        if posindices(end)==maxposindex
            error(['Error in index list: Index ' num2str(maxposindex) ' appears more than twice']);
        else
            error(['Error in index list: Index ' num2str(maxposindex) ' only appears once']);
        end
    end
    altposindices = posindices(1:2:numel(posindices));
    flags = altposindices(1:end-1)==altposindices(2:end);
    if any(flags)
        errorpos = find(flags,1,'first');
        error(['Error in index list: Index ' num2str(altposindices(errorpos)) ' appears more than twice']);
    end
    % Check sizes match
    sizes = ones(size(allindices));
    ptr = 1;
    for a=1:numel(tensorList)
        sz = size(tensorList{a});
        if numel(legLinks{a})==1 % Is a vector (1D)
            sz = max(sz);
        end
        sizes(ptr:ptr+numel(sz)-1) = sz;
        ptr = ptr + numel(legLinks{a});
    end
    sizes = sizes(ix); % Sort in ascending positive leg sequence
    flags = sizes(1:2:end)~=sizes(2:2:end);
    if numdummytensors==1
        for a=legLinks{dummytensorpos}
            flags(altposindices==a) = false; % Prevent size mismatch errors for dummy tensor
        end
    end
    if any(flags)
        errorpos = find(flags,1,'first');
        error(['Leg size mismatch on index ' num2str(altposindices(errorpos))]);
    end
    sizes = sizes(1:2:end);
    if numdummytensors==1
        % Ensure correct sizes recorded for dummy tensor legs
        realtensors = 1:numel(tensorList);
        realtensors(dummytensorpos) = [];
        for a=realtensors
             found = ismember(legLinks{a},legLinks{dummytensorpos});
             for b=find(found)
                 sizes(legLinks{a}(found)) = size(tensorList{a},b);
             end
        end
    end
    orderLegs = sparse([],[],[],1,altposindices(end),numel(altposindices));
    orderLegs(altposindices)=1:numel(altposindices);
    for a=1:numel(legLinks)
        if ~isempty(tensorList{a})
            if prod(sizes(orderLegs(legLinks{a})))~=numel(tensorList{a})
                error(['Size of tensor ' num2str(a) ' does not match the product of its labelled legs - is there an entry missing from the leg list?']);
            end
        end
    end
    % Shape 1D tensors as column vectors
    for a=1:numel(legLinks)
        if numel(legLinks{a})==1
            if numel(tensorList{a})~=size(tensorList{a},1)
                warning('multienv:rowvec',['Tensor ' num2str(a) ' is a 1D object and was supplied as a row vector. Converted to column vector.']);
                tensorList{a} = reshape(tensorList{a},[],1);
            end
        end
    end
    % Check envList is sane - numerical, correct size, and contains at least one non-zero entry
    if ~isnumeric(envList)
        error('Environment list must be a 1xn numeric array');
    end
    if ~isequal(size(envList),size(tensorList))
        error('Size of requested environment list does not match size of tensor list')
    end
    if ~isreal(envList) || any(envList<0) || any(uint32(envList)~=envList)
        error('The environment list must contain only non-negative integers and zeros')
    end
    if all(envList==0)
        error('The environment list must contain at least one entry which is not zero')
    end
    envListCopy = envList;
    while any(envListCopy~=0)
        newlist = find(envListCopy==envListCopy(find(envListCopy~=0,1)));
        tgtsize = size(tensorList{newlist(1)});
        for a=2:numel(newlist)
            if ~isequal(size(tensorList{newlist(a)}),tgtsize)
                error(['All tensors carrying environment label ' num2str(envListCopy(newlist(1))) ' must be the same size']);
            end
        end
        envListCopy(newlist)=0;
    end
    if exist('sequence','var')
        % Check sequence is a row vector of positive real integers, each occurring only once, and zeros.
        % Check they match the positive leg labels.
        if any(uint32(sequence)~=sequence)
            error('All entries in contraction sequence must be real positive integers or zero');
        end
        if numel(altposindices)~=sum(sequence>0)
            error('Each positive index must appear once and only once in the contraction sequence');
        end
        if any(altposindices~=sort(sequence(sequence>0)))
            error('Each positive index must appear once and only once in the contraction sequence');
        end
    end
end

function [datatree,finaltraces] = constructtree(legLinks,sequence,envList,legDims,dummyTensorPos,orderLegs)
    finaltraces = [];
    dummyTensorFlag = false(1,numel(legLinks));
    if dummyTensorPos~=0
        dummyTensorFlag(dummyTensorPos) = true;
    end
    
    vertex.objects = cell(1,3);
    vertex.objectLegLinks = cell(1,3);
    vertex.neededfor = cell(1,3); % What each output object of this vertex is needed for
    vertex.connect = zeros(3,2); % Connections to next vertex: Index 1 - leg of this vertex. Index 2 - which target vertex, which leg.
    
    objectID = [-(1:numel(legLinks)); zeros(1,numel(legLinks))].'; % First index is object ID, second index is "output" leg of that object
    
    datatree = {};
    
    % Create initial datatree
    if isempty(sequence) && numel(legLinks)>1
        sequence = zeros(1,numel(legLinks)-1); % If no contraction sequence, do outer products
    end
    while ~isempty(sequence)
        if sequence(1)==0
            % Outer product - Identify pair of tensors to contract
            if all(sequence==0) % Final outer product of all remaining objects - ensure enough zeros are present in the sequence
                sequence = zeros(1,numel(legLinks)-1);
            end
            % - Determine number of outer products pending
            numOPs = 1;
            while sequence(numOPs)==0 && numOPs < numel(sequence)
                numOPs = numOPs+1;
            end
            if sequence(numOPs)~=0
                numOPs = numOPs - 1;
            end
            % - Determine list of tensors on which OP is to be performed
            if numOPs == numel(legLinks)-1
                % OP of all remaining tensors
                OPlist = 1:numel(legLinks);
                if numel(sequence) > numOPs
                    finaltraces = sequence(numOPs+1:end); % Note any contractions to be performed as traces on the final object after all contractions are complete.
                                                          % In practice these are carried out as part of the pairwise fusion sequence as this is more
                                                          % efficient, and a warning is generated if their dimension is greater than 1.
                    sequence = sequence(1:numOPs); % Remove these final traces from the sequence
                end
            else
                % For OP of n tensors (n=numOPs+1) when more than n tensors remain, proceed past the zeros in the sequence and read nonzero indices until 
                % n+1 tensors accounted for, or the dummy tensor is encountered. Failure to achieve either of these conditions implies an invalid sequence.
                flags = false(1,numel(legLinks));
                ptr = numOPs+1;
                while sum(flags) < numOPs+2
                    % Flag tensors on which leg given by sequence(ptr) appears
                    count = 0;
                    for a=1:numel(legLinks)
                        if any(legLinks{a}==sequence(ptr))
                            flags(a) = true;
                            count = count + 1;
                        end
                    end
                    if count~=2
                        error(['Reading indices after an outer product, found an index ' num2str(sequence(ptr)) ' which appears on ' num2str(count) ' tensor(s). This violates the rules of zeros-in-sequence notation.']);
                    end
                    ptr = ptr + 1;
                end
                % Identify which of these tensors is _not_ participating in the OP (but is instead contracted with the result of the OP), and unflag it.
                % If one of them is the dummy tensor, that's it. Otherwise:
                % - Identify the two tensors on which the first nonzero index appears
                % - Examine consecutive nonzero indices until one matches only one of the two tensors. This is the tensor to unflag.
                if any(flags & dummyTensorFlag)
                    flags = flags & (~dummyTensorFlag);
                else
                    firsttensors = [0 0];
                    ptr = numOPs+1;
                    for a=1:numel(legLinks)
                        if any(legLinks{a}==sequence(ptr))
                            if firsttensors(1)==0
                                firsttensors(1) = a;
                            else
                                firsttensors(2) = a;
                                break;
                            end
                        end
                    end
                    done = false;
                    while ~done
                        nexttensors = [0 0];
                        ptr = ptr + 1;
                        for a=1:numel(legLinks)
                            if any(legLinks{a}==sequence(ptr))
                                if nexttensors(1)==0
                                    nexttensors(1) = a;
                                else
                                    nexttensors(2) = a;
                                    break;
                                end
                            end
                        end
                        if ~isequal(firsttensors,nexttensors)
                            done = true;
                        end
                    end
                    if any(firsttensors == nexttensors(1))
                        flags(nexttensors(1)) = false;
                    else
                        flags(nexttensors(2)) = false;
                    end
                end
                OPlist = find(flags);
            end
            % - Identify the first two tensors from this list on which an OP is to be performed (the two of smallest total dimension).
            OPsizes = zeros(1,numel(OPlist));
            for a=1:numel(OPlist)
                OPsizes(a) = prod(legDims(orderLegs(legLinks{OPlist(a)})));
            end
            [~, ix] = sort(OPsizes,'ascend');
            tensors = sort(OPlist(ix(1:2)));
        else
            % Find the two tensors connected by sequence(1)
            ptr = 1;
            while ~any(legLinks{ptr}==sequence(1))
                ptr = ptr + 1;
            end
            if sum(legLinks{ptr}==sequence(1))==2
                error('Inappropriate trace encountered - when contracting two tensors, ensure that all indices connecting these tensors appear consecutively in the sequence.')
                % This error condition should never be reached - it implies that some common indices were left uncontracted when two tensors were combined.
            end
            % Not a trace - find both tensors
            tensors = [ptr 0];
            ptr = ptr + 1;
            while ~any(legLinks{ptr}==sequence(1))
                ptr = ptr + 1;
            end
            tensors(2) = ptr;
            tensors = sort(tensors);
        end
        % Common to outer-product and index-linked contractions:
        % Identify indices linking these two tensors ("linking" indices) and indices not linking these two tensors ("free" indices)
        linkingLegs = legLinks{tensors(1)};
        freeLegs1 = legLinks{tensors(1)};
        freeLegs2 = legLinks{tensors(2)};
        for a=1:numel(legLinks{tensors(1)})
            if any(legLinks{tensors(2)}==legLinks{tensors(1)}(a))
                freeLegs1(freeLegs1==legLinks{tensors(1)}(a)) = [];
                freeLegs2(freeLegs2==legLinks{tensors(1)}(a)) = [];
            else
                linkingLegs(linkingLegs==legLinks{tensors(1)}(a)) = [];
            end
        end
        if ~isequal(sort(linkingLegs),sort(sequence(1:min([numel(linkingLegs) numel(sequence)]))))
            fromseq = sequence(1:min([numel(linkingLegs) numel(sequence)]));
            for a=numel(fromseq):-1:1
                if ~any(linkingLegs==fromseq(a))
                    fromseq(a) = [];
                end
            end
            for a=numel(linkingLegs):-1:1
                if any(fromseq==linkingLegs(a))
                    linkingLegs(a) = [];
                end
            end
            if any(legDims(orderLegs(linkingLegs))~=1) && ~any(dummyTensorFlag(tensors))
                % Warn for suboptimal sequences. Note that sequences which defer contraction over indices of trivial dimension are
                % acceptable and don't generate a warning (although in practice these indices will still be contracted at the earliest opportunity).
                if ~isempty(fromseq)
                    t = 'Sequence suboptimal: When contracting ind';
                    if numel(fromseq)==1
                        t = [t 'ex ' num2str(fromseq) ' please also contract ind']; %#ok<AGROW>
                    else
                        t = [t 'ices ' num2str(fromseq) ' please also contract ind']; %#ok<AGROW>
                    end
                    if numel(linkingLegs)==1
                        t = [t 'ex ' num2str(linkingLegs) ' as these indices connect the same two tensors.']; %#ok<AGROW>
                    else
                        t = [t 'ices ' num2str(linkingLegs) ' as these indices connect the same two tensors.']; %#ok<AGROW>
                    end
                    warning('multienv:suboptimalsequence',[t ' This has been corrected automatically.']); % 
                else
                    t = 'Sequence suboptimal: Instead of performing an outer product, please contract ind';
                    if numel(linkingLegs)==1
                        t = [t 'ex ' num2str(linkingLegs) '. This index connects the same two tensors and is non-trivial.']; %#ok<AGROW>
                    else
                        t = [t 'ices ' num2str(linkingLegs) '. These indices connect the same two tensors and are non-trivial.']; %#ok<AGROW>
                    end
                    warning('multienv:suboptimalsequence',[t ' This has been corrected automatically.']); % 
                end
            end
            linkingLegs = [fromseq linkingLegs]; %#ok<AGROW>
        end
        % Create tree node
        datatree{end+1} = vertex; %#ok<AGROW>
        datatree{end}.connect(1,:) = objectID(tensors(1),:);
        datatree{end}.connect(2,:) = objectID(tensors(2),:);
        if objectID(tensors(1),1)>0
            % Link from previous node
            datatree{objectID(tensors(1),1)}.connect(objectID(tensors(1),2),:) = [numel(datatree) 1]; %#ok<AGROW>
        end
        if objectID(tensors(2),1)>0
            % Link from previous node
            datatree{objectID(tensors(2),1)}.connect(objectID(tensors(2),2),:) = [numel(datatree) 2]; %#ok<AGROW>
        end
        % Update legLinks - remove objects just combined
        legLinks{tensors(1)} = [freeLegs1 freeLegs2];
        legLinks(tensors(2)) = [];
        % Update dummy tensor flag
        dummyTensorFlag(tensors(2)) = []; % It's neither of the tensors being combined, so just delete the one removed from the list to keep the position of the dummy tensor flag correct
        % Update objectID
        objectID(tensors(1),:) = [numel(datatree) 3];
        objectID = objectID([1:tensors(2)-1 tensors(2)+1:end],:);
        % Move to next uncontracted entry of sequence
        if sequence(1)==0 % If OP, it's done.
            sequence(1) = [];
        end
        for a=1:numel(linkingLegs) % If linked by shared indices (whether or not contraction is designated as an OP), these also have now been done.
            sequence(sequence==linkingLegs(a))=[];
        end
        if isempty(sequence)
            % All specified contractions completed - is the network still disjoint?
            if numel(legLinks)>1
                sequence = zeros(1,numel(legLinks)-1);
            end
        end
    end
    if ~isempty(datatree)
        % Tidy up end of tree
        % datatree{end-1}.connect(3,:) = datatree{end}.connect(3-datatree{end-1}.connect(3,2),:);
        if datatree{end}.connect(1,1) > 0
            datatree{datatree{end}.connect(1,1)}.connect(datatree{end}.connect(1,2),:) = datatree{end}.connect(2,:);
        end
        if datatree{end}.connect(2,1) > 0
            datatree{datatree{end}.connect(2,1)}.connect(datatree{end}.connect(2,2),:) = datatree{end}.connect(1,:);
        end
        datatree = datatree(1:end-1);
    end
    
    % Construct list of what each object is needed for
    for a=1:numel(datatree)
        for b=reshape(find(datatree{a}.connect(:,1)<0),1,[])
            if envList(-datatree{a}.connect(b,1))>0
                datatree = createneedslist(datatree,a,b,datatree{a}.connect(b,1));
            end
        end
    end
    
end

function datatree = createneedslist(datatree,node,link,flagas)
    % List of legs to follow in update process
    update = 1:3;
    update(link) = [];
    targetnodes = datatree{node}.connect(update,1);
    targetlinks = datatree{node}.connect(update,2);
    targetlinks(targetnodes<0) = [];
    targetnodes(targetnodes<0) = []; % Don't try to update neededfor entries on the "nodes" associated with original tensors!
    for a=1:numel(targetnodes)
        if ~any(datatree{targetnodes(a)}.neededfor{targetlinks(a)}==flagas) % If not previously done...
            % Leg "targetlinks(a)" on node "targetnodes(a)" enters node "node" on leg "link". Flag as required.
            datatree{targetnodes(a)}.neededfor{targetlinks(a)}(end+1) = flagas;
            % Flag as required the necessary inputs for the object on leg "targetlinks(a)" of node "targetnode(a)".
            datatree = createneedslist(datatree,targetnodes(a),targetlinks(a),targetnodes(a)+targetlinks(a)/10);
        end
    end
end

function [newEnv,datatree] = calculateenvironment(tensorList,legLinks,datatree,which)
    if numel(datatree)==0
        % Only two tensors - no tree. Any traces have already been performed, so environment of tensorList{which} is given by tensorList{3-which}, up to a possible permute.
        newEnv = tensorList{3-which};
        % Reorder legs of newEnv to match up with order of tensorList{which}
        seq = zeros(size(legLinks{3-which}));
        if numel(seq)>1
            for a=1:numel(seq)
                seq(a) = find(legLinks{3-which}==legLinks{which}(a));
            end
            newEnv = permute(newEnv,seq);
        end
    else
        % Find (in the tree) the object whose environment is to be calculated
        for a=1:numel(datatree)
            link = datatree{a}.connect(:,1)==-which;
            if any(link)
                node = a;
                break;
            end
        end
        % Request the objects required to compute the new environment
        required = 1:3;
        required(link) = [];
        node1 = datatree{node}.connect(required(1),1);
        link1 = datatree{node}.connect(required(1),2);
        [tensor1,legs1,datatree] = requestobject(datatree,tensorList,legLinks,node1,link1);
        node2 = datatree{node}.connect(required(2),1);
        link2 = datatree{node}.connect(required(2),2);
        [tensor2,legs2,datatree] = requestobject(datatree,tensorList,legLinks,node2,link2);
        % Compute the environment
        % - Set up output legs
        for a=1:numel(legLinks{which})
            legs1(legs1==legLinks{which}(a)) = -a;
            legs2(legs2==legLinks{which}(a)) = -a;
        end
        % - Check that there are no unshared positive index labels
        for a=legs1(legs1>0)
            if ~any(legs2==a)
                error('Error in implementation of multienv: Unshared positive index labels in final contraction!')
            end
        end
        % - Call tcontract to contract the two tensors
        % disp(['Env of ' num2str(which) '. Cost: X^' num2str(numel(legs1)+sum(legs2<0))]);
        newEnv = tcontract(tensor1,tensor2,legs1,legs2);
        % Update neededfor
        if node1>0
            datatree{node1}.neededfor{link1} = datatree{node1}.neededfor{link1}(datatree{node1}.neededfor{link1}~=-which);
            if isempty(datatree{node1}.neededfor{link1})
                datatree{node1}.objects{link1} = [];
                datatree{node1}.objectLegLinks{link1} = [];
            end
        end
        if node2>0
            datatree{node2}.neededfor{link2} = datatree{node2}.neededfor{link2}(datatree{node2}.neededfor{link2}~=-which);
            if isempty(datatree{node2}.neededfor{link2})
                datatree{node2}.objects{link2} = [];
                datatree{node2}.objectLegLinks{link2} = [];
            end
        end
    end
end

function [tensor,legs,datatree] = requestobject(datatree,tensorList,legLinks,node,link)
    if node<0
        % Return object from tensorList
        tensor = tensorList{-node};
        legs = legLinks{-node};
    else
        if ~isempty(datatree{node}.objects{link})
            % Return object from datatree
            tensor = datatree{node}.objects{link};
            legs = datatree{node}.objectLegLinks{link};
        else
            % Construct and return object from datatree
            % - Get prerequisites
            required = 1:3;
            required(link) = [];
            node1 = datatree{node}.connect(required(1),1);
            link1 = datatree{node}.connect(required(1),2);
            [tensor1,legs1,datatree] = requestobject(datatree,tensorList,legLinks,node1,link1);
            node2 = datatree{node}.connect(required(2),1);
            link2 = datatree{node}.connect(required(2),2);
            [tensor2,legs2,datatree] = requestobject(datatree,tensorList,legLinks,node2,link2);
            % - Compute, store, and return the object
            % -- Identify common legs (to be used as contraction sequence) and set up output legs
            linkingLegs = legs1;
            freeLegs1 = legs1;
            freeLegs2 = legs2;
            for a=1:numel(legs1)
                if any(legs2==legs1(a))
                    freeLegs1(freeLegs1==legs1(a)) = [];
                    freeLegs2(freeLegs2==legs1(a)) = [];
                else
                    linkingLegs(linkingLegs==legs1(a)) = [];
                end
            end
            for a=1:numel(freeLegs1)
                legs1(legs1==freeLegs1(a)) = -a;
            end
            for a=1:numel(freeLegs2)
                legs2(legs2==freeLegs2(a)) = -a-numel(freeLegs1);
            end
            % -- Contract tensors
            % disp(['Object ' num2str(node) '.' num2str(link) '. Cost: X^' num2str(numel(legs1)+sum(legs2<0))]);
            tensor = tcontract(tensor1,tensor2,legs1,legs2);
            legs = [freeLegs1 freeLegs2];
            % -- Store result - if no longer needed after this calculation, will be cleared by calling function
            datatree{node}.objects{link} = tensor;
            datatree{node}.objectLegLinks{link} = legs;
            % - Update neededfor
            flag = node + link/10;
            if node1>0
                datatree{node1}.neededfor{link1} = datatree{node1}.neededfor{link1}(datatree{node1}.neededfor{link1}~=flag);
                if isempty(datatree{node1}.neededfor{link1})
                    datatree{node1}.objects{link1} = [];
                    datatree{node1}.objectLegLinks{link1} = [];
                end
            end
            if node2>0
                datatree{node2}.neededfor{link2} = datatree{node2}.neededfor{link2}(datatree{node2}.neededfor{link2}~=flag);
                if isempty(datatree{node2}.neededfor{link2})
                    datatree{node2}.objects{link2} = [];
                    datatree{node2}.objectLegLinks{link2} = [];
                end
            end
        end
    end
end

function tensor = tcontract(T1,T2,legs1,legs2)
    % If either tensor is a number (no legs), add a trivial leg.
    if numel(legs1)==0
        legs1 = max(abs(legs2))+1;
        legs2 = [legs2 legs1];
    else
        if numel(legs2)==0
            legs2 = max(abs(legs1))+1;
            legs1 = [legs1 legs2];
        end
    end
    
    freelegs1 = find(legs1<0);
    [~, ix] = sort(legs1(freelegs1),'descend');
    freelegs1 = freelegs1(ix); % Sort free legs list in descending order (increasing magnitude)
    freelegs2 = find(legs2<0);
    [~, ix] = sort(legs2(freelegs2),'descend');
    freelegs2 = freelegs2(ix); % Sort free legs list in descending order (increasing magnitude)
    linklegs1 = find(legs1>0);
    linklegs2 = zeros(size(linklegs1));
    for a=1:numel(linklegs1)
        linklegs2(a) = find(legs2==legs1(linklegs1(a)),1);
    end
    sz1 = [size(T1) ones(1,numel(legs1)-ndims(T1))];
    sz2 = [size(T2) ones(1,numel(legs2)-ndims(T2))];
    if numel(legs1)>1
        T1 = permute(T1,[freelegs1 linklegs1]);
    end
    if numel(legs2)>1
        T2 = permute(T2,[linklegs2 freelegs2]);
    end
    linksize = prod(sz1(linklegs1)); % NB prod([]) = 1 if no linked legs
    T1 = reshape(T1,prod(sz1(freelegs1)),linksize);
    T2 = reshape(T2,linksize,prod(sz2(freelegs2)));
    tensor = T1 * T2;
    tensor = reshape(tensor,[sz1(freelegs1) sz2(freelegs2) 1 1]);
    if (numel(freelegs1)+numel(freelegs2)) > 1
        % Because of sorting the free leg list order earlier, there is a 
        % chance we may not need the final permute:
        leglist = -[legs1(freelegs1) legs2(freelegs2)];
        if ~isequal(leglist,1:numel(leglist))
            tensor = ipermute(tensor,leglist);
        end
    end
end
