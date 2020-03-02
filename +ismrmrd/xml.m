function header = xml(fname, schema)
% Read the flexible header of an ismrmrd file into a structure.
%
% FORMAT header = ismrmrd.xml(in,[schema])
%
% The input can be the path to an H5 or XML file, an XML string or a 
% DOM object.
% An XML schema can be provided to validate the XML. By default, th ISMRMRD
% schema is used.
% The output is a structure.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% -------------------------------------------------------------------------
% Read schema
if nargin < 2
    schema = '';
end
if ~isstruct(schema)
    schema = ismrmrd.schema(schema);
end

% -------------------------------------------------------------------------
% Read header
if ischar(fname)
    [~,~,ext] = fileparts(fname);
    switch ext
        case '.h5'
            if ~exist(fname,'file')
                error('File %s does not exist',fname);
            end
            xmlstring = h5read(fname, '/dataset/xml');
            xmlstring = xmlstring{1};
            DOMnode   = xmlreadstring(xmlstring);
        case '.xml'
            if ~exist(fname,'file')
                error('File %s does not exist',fname);
            end
            DOMnode   = xmlread(fname);
        otherwise 
            DOMnode   = xmlreadstring(fname);
    end
elseif isa(fname, 'org.apache.xerces.dom.DeferredDocumentImpl')
    DOMnode = fname;
else
    error('The input must be a h5 file, xml file, xml string or DOM object');
end

% -------------------------------------------------------------------------
% Parse header
header = parseNode(DOMnode, schema);

% =========================================================================
function obj = parseNode(DOMnode, schema)
% Recurse over node children.

nodeName   = char(DOMnode.getNodeName);

if isfield(schema, 'parse') && ~isempty(schema.parse)
    
    % ---------------------------------------------------------------------
    % Parsable value (= leaf)
    
    obj = [];
    
    childNodes    = DOMnode.getChildNodes;
    numChildNodes = childNodes.getLength;
    for c=1:numChildNodes
        theChild  = childNodes.item(c-1);
        childName = char(theChild.getNodeName);
        if childName(1) == '#'
            if ~isempty(char(theChild.getData))
                switch childName(2:end)
                    case 'text'
                        value = schema.parse(char(theChild.getData));
                        if isempty(obj)
                            obj = value;
                        else
                            obj(end+1) = value;
                        end
                    otherwise
                        warning('Unknown type %s',childName)
                end
            end
            continue
        end
    end
    if isfield(schema, 'minOccurs')
        if isfield(schema, 'default')
            value = schema.parse(schema.default);
            for i=1:(schema.minOccurs-numel(obj))
                if isempty(obj)
                    obj = value;
                else
                    obj(end+1) = value;
                end
            end
        end
        if numel(obj) < schema.minOccurs
            warning('Not enough values');
        end
    end
    if isfield(schema, 'maxOccurs')
        if numel(obj) > schema.maxOccurs
            warning('Too many values for %s (%d/%d)', nodeName, numel(obj), schema.maxOccurs);
            obj = obj(1:schema.maxOccurs);
        end
    end
    if numel(obj) == 1 && iscell(obj)
        obj = obj{1};
    end
    
else
    
    % ---------------------------------------------------------------------
    % Non-parsable value (= node)
    
    obj = struct;

    protected = {'type' 'parse' 'default' 'maxOccurs' 'minOccurs' ...
                 'enumeration' 'fractionDigits' 'length' ...
                 'maxExclusive' 'maxInclusive' 'maxLength' ...
                 'minExclusive' 'minInclusive' 'minLength' ...
                 'pattern' 'totalDigits' 'extension'};
             
    fields        = fieldnames(schema);
    childNodes    = DOMnode.getChildNodes;
    numChildNodes = childNodes.getLength;
    for c=1:numChildNodes
        theChild  = childNodes.item(c-1);
        childName = char(theChild.getNodeName);
        childName = strsplit(childName, ':');
        childName = childName{end};
        
        if any(strcmpi(childName, fields)) && ...
           ~any(strcmpi(childName, protected))
       
            child = parseNode(theChild, schema.(childName));
       
            if ~isfield(obj, childName)
                obj.(childName) = child;
            else
                obj.(childName)(end+1) = child;
            end
       
        end
        
    end
    
end