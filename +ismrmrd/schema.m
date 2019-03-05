function schema = schema(fname)
% Read the ISMRMRD XML schema into a matlab structure.
%
% FORMAT obj = ismrmrd.schema(in)
%
% The input can be the path to an XSD file, an XSD string or a DOM object.
% The output is a structure.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% -------------------------------------------------------------------------
% Read schema
if nargin < 1 || isempty(fname)
    dir   = fileparts(mfilename('fullpath'));
    fname = fullfile(dir, 'schema.xsd');
end
if ischar(fname)
    [~,~,ext] = fileparts(fname);
    switch ext
        case '.xsd'
            DOMnode   = xmlread(fname);
        otherwise 
            DOMnode   = xmlreadstring(fname);
    end
elseif isa(fname, 'org.apache.xerces.dom.DeferredDocumentImpl')
    DOMnode = fname;
end

% -------------------------------------------------------------------------
% Parse schema -> to struct
schema = parsNode(DOMnode);

% -------------------------------------------------------------------------
% Link named types
schema = linkTypes(schema);


% =========================================================================
function obj = linkTypes(obj, obj0)

if ~isstruct(obj)
    return
end
if nargin < 2
    obj0 = obj;
end

if (isfield(obj, 'parse') && isempty(obj.parse)) || isfield(obj, 'extension')
    if isfield(obj, 'type') && isfield(obj0, 'simpleType') && isfield(obj0.simpleType, obj.type)
        simpletype = obj.type;
        obj0.simpleType.(simpletype) = linkTypes(obj0.simpleType.(simpletype), obj0);
        typefields = fieldnames(obj0.simpleType.(simpletype));
        for cc=1:numel(typefields)
            obj.(typefields{cc}) = obj0.simpleType.(simpletype).(typefields{cc});
%                     obj.(typefields{cc}) = linkTypes(obj.(typefields{cc}), obj0);
        end
    end
    if isfield(obj, 'type') && isfield(obj0, 'complexType') && isfield(obj0.complexType, obj.type)
        complextype = obj.type;
        obj0.complexType.(complextype) = linkTypes(obj0.complexType.(complextype), obj0);
        typefields = fieldnames(obj0.complexType.(complextype));
        for cc=1:numel(typefields)
            obj.(typefields{cc}) = obj0.complexType.(complextype).(typefields{cc});
%                     obj.(typefields{cc}) = linkTypes(obj.(typefields{cc}), obj0);
        end
        if isfield(obj, 'parse')
            obj = rmfield(obj, 'parse');
        end
    end
    if isfield(obj, 'extension') && isfield(obj0, 'complexType') && isfield(obj0.complexType, obj.extension)
        complextype = obj.extension;
        obj0.complexType.(complextype) = linkTypes(obj0.complexType.(complextype), obj0);
        typefields = fieldnames(obj0.complexType.(complextype));
        for cc=1:numel(typefields)
            if ~isfield(obj, typefields{cc})
                obj.(typefields{cc}) = obj0.complexType.(complextype).(typefields{cc});
%                         obj.(typefields{cc}) = linkTypes(obj.(typefields{cc}), obj0);
            end
        end
        if isfield(obj, 'parse')
            obj = rmfield(obj, 'parse');
        end
    end
end

fields = fieldnames(obj);
for c=1:numel(fields)
    field = fields{c};
    if ~any(strcmpi(field, {'type' 'simpleType' 'complextype' 'extension' 'parse'}))
        obj.(field) = linkTypes(obj.(field), obj0);
    end
end

if isfield(obj, 'simpleType'),  obj = rmfield(obj, 'simpleType');  end
if isfield(obj, 'complexType'), obj = rmfield(obj, 'complexType'); end


% =========================================================================
function obj = parsNode(DOMnode, obj)
% Recurse over node children.

if nargin < 2
    obj = struct;
end

% -------------------------------------------------------------------------
% Parse name
nodeName   = char(DOMnode.getNodeName);
if ~strcmpi(nodeName(1:3), 'xs:')
    childNodes    = DOMnode.getChildNodes;
    numChildNodes = childNodes.getLength;
    for c=1:numChildNodes
        DOMnode = childNodes.item(c-1);
        obj     = parsNode(DOMnode, obj);
    end
    return
end
nodeName   = nodeName(4:end);

% -------------------------------------------------------------------------
% Parse attributes
DOMattrs   = DOMnode.getAttributes;
numAttr    = DOMattrs.getLength;
attributes = struct;
for c=1:numAttr
    DOMattr   = DOMattrs.item(c-1);
    attrName  = char(DOMattr.getName);
    attrValue = char(DOMattr.getValue);
    attrName  = strrep(attrName, ':', '_');
    attributes.(attrName) = attrValue;
end

% -------------------------------------------------------------------------
% Format
switch nodeName
    case 'schema'
        % This is the root node, let's skip it.
        childNodes    = DOMnode.getChildNodes;
        numChildNodes = childNodes.getLength;
        for c=1:numChildNodes
            DOMnode = childNodes.item(c-1);
            obj = parsNode(DOMnode, obj);
        end
        
    case 'element'
        % This is usually a leaf in the tree. However, the DOM might have
        % a child when the element type is not basic. In this case, the
        % child is a 'simpleType' or 'complexType'.
        name       = attributes.name;
        attributes = rmfield(attributes, 'name');
        obj.(name) = attributes;
        fields     = fieldnames(obj.(name));
        for c=1:numel(fields)
            field = fields{c};
            switch(field)
                case {'minOccurs' 'maxOccurs'}
                    switch obj.(name).(field)
                        case 'unbounded'
                            obj.(name).(field) = inf;
                        otherwise
                            obj.(name).(field) = str2double(obj.(name).(field));
                    end
                case 'type'
                    [parse, type] = parsefunc(obj.(name).(field));
                    obj.(name).(field) = type;
                    obj.(name).parse   = parse;
            end
        end
        childNodes    = DOMnode.getChildNodes;
        numChildNodes = childNodes.getLength;
        for c=1:numChildNodes
            DOMnode = childNodes.item(c-1);
            obj.(name) = parsNode(DOMnode, obj.(name));
        end
        
    case 'complexType'
        % A complex element is an element that contains other elements or
        % attributes. Therefore, the type needs to describe these
        % sub-elements.
        % Usually, its unique child will be one of 'all', 'sequence' or
        % 'complexContent'. However, for empty elements that contain only
        % attributes, the 'complexElement' part might be skipped, in which
        % case children will be of type 'attribute'.
        childNodes    = DOMnode.getChildNodes;
        numChildNodes = childNodes.getLength;
        tmp = struct;
        for c=1:numChildNodes
            DOMnode = childNodes.item(c-1);
            tmp = parsNode(DOMnode, tmp);
        end
        if isfield(attributes, 'name')
            name = attributes.name;
            obj.complexType.(name) = tmp;
        else
            fields = fieldnames(tmp);
            for c = 1:numel(fields)
               obj.(fields{c}) = tmp.(fields{c});
            end
        end
        
    case 'simpleType'
        % A simple element is an element that does not contains other 
        % elements or attributes. However, it may define non-trivial
        % validation rules. A 'simpleType' will usually contain a unique
        % 'restriction' child.
        childNodes    = DOMnode.getChildNodes;
        numChildNodes = childNodes.getLength;
        tmp = struct;
        for c=1:numChildNodes
            DOMnode = childNodes.item(c-1);
            tmp = parsNode(DOMnode, tmp);
        end
        if isfield(attributes, 'name')
            name = attributes.name;
            obj.simpleType.(name) = tmp;
        else
            fields = fieldnames(tmp);
            for c = 1:numel(fields)
               obj.(fields{c}) = tmp.(fields{c});
            end
        end
        
    case {'sequence' 'all'}
        % A 'sequence' contains an ordered list of elements, whereas 'all'
        % contains an unordered list of elements. This is related to
        % validation, which we do not do thoroughly here. We thus consider
        % these two node types as equivalent.
        childNodes    = DOMnode.getChildNodes;
        numChildNodes = childNodes.getLength;
        for c=1:numChildNodes
            DOMnode = childNodes.item(c-1);
            obj = parsNode(DOMnode, obj);
        end
        
    case 'restriction'
        % A restriction is a validation rule, that restricts the range of
        % possible values.
        if isfield(attributes, 'base')
            type = attributes.base;
        else
            type = '';
        end
        [obj.parse, obj.type] = parsefunc(type);
        childNodes    = DOMnode.getChildNodes;
        numChildNodes = childNodes.getLength;
        for c=1:numChildNodes
            DOMnode = childNodes.item(c-1);
            obj = parsNode(DOMnode, obj);
        end
        
    case 'enumeration'
        % An enumeration is restriction that lists all possible values.
        if ~isfield(obj, 'enumeration')
            obj.enumeration = {};
        end
        obj.enumeration{end+1} = obj.parse(attributes.value);
        
    case {'fractionDigits' 'length' ...
          'maxExclusive' 'maxInclusive' 'maxLength' ...
          'minExclusive' 'minInclusive' 'minLength' ...
          'pattern' 'totalDigits'}
        % These are all basic restrictions which acts on the limits or size
        % of the value. Some work esclusively with numbers or strings, some
        % work with both. Pattern contains a regex pattern to validate the
        % structure of a string.
        obj.(nodeName) = obj.parse(attributes.value);
        
    case 'whiteSpace'
        % If 'collapse', this enables a post-processing step where
        % redondant or trailing white spaces are removed.
        if strcmpi(attributes.value, 'collapse')
            % That's not enough. I should also convert all white-space
            % characters to spaces and collapse multiple spaces.
            obj.parse = @(X) obj.parse(deblank(X));
        end
        
    case {'complexContent' 'simpleContent'}
        % Content of a complexType.
        % Children:
        % ? annotation
        % 1 (restriction | extension)
        childNodes    = DOMnode.getChildNodes;
        numChildNodes = childNodes.getLength;
        for c=1:numChildNodes
            DOMnode = childNodes.item(c-1);
            obj = parsNode(DOMnode, obj);
        end
        
    case 'extension'
        % Specifies that a type is an extension of another type and
        % therefore inherits its members.
        % Attributes:
        % ? base = QName
        % Children:
        % ? annotation
        % * (attribute | attributeGroup)
        extension     = strsplit(attributes.base, ':');
        obj.extension = extension{end};
        childNodes    = DOMnode.getChildNodes;
        numChildNodes = childNodes.getLength;
        for c=1:numChildNodes
            DOMnode = childNodes.item(c-1);
            obj = parsNode(DOMnode, obj);
        end
        
    case 'attribute'
        % Allows to specify requirements for the attributes of XML nodes
        % Attributes:
        % ? default = string
        % ? fixed   = string
        % ? form    = (qualified | unqualified)
        % ? id      = ID
        % ? name    = NCName
        % ? ref     = QName
        % ? type    = QName
        % ? use     = (optional | prohibited | required) : optional
        % Children:
        % ? annotation
        % ? simpleType
        name = strrep(attributes.name, '-', '_');
        obj.attributes.(name) = attributes;
        if ~isfield(obj.attributes.(name), 'use')
            obj.attributes.(name).use = 'optional';
        end
        if isfield(obj.attributes.(name), 'type')
            [parse, type] = parsefunc(obj.attributes.(name).type);
            obj.attributes.(name).type  = type;
            obj.attributes.(name).parse = parse;
        end
        childNodes    = DOMnode.getChildNodes;
        numChildNodes = childNodes.getLength;
        for c=1:numChildNodes
            DOMnode = childNodes.item(c-1);
            obj.attributes.(name) = parsNode(DOMnode, obj.attributes.(name));
        end
        
    case 'annotation'
        % Provide human-readable information.
        % We don't really care here.
end

% =========================================================================
% Parse function associates conversion functions with basic types
% I have only implemented types that are used in the current version of the
% ismrmrd schema, but this could be easily extended to all types supported
% by the xml standard.
function [parse, type] = parsefunc(type)
    if numel(type) >= 3 && strcmpi(type(1:3), 'xs:')
        type = type(4:end);
        switch type
            case 'unsignedShort', parse = @(X) uint16(str2double(char(X)));
            case 'long',          parse = @(X) int64(str2double(char(X)));
            case 'float',         parse = @(X) single(str2double(char(X)));
            case 'double',        parse = @str2double;
            case 'string',        parse = @char;
            case 'date',          parse = @(X) datetime(char(X), 'InputFormat', 'yyyy-MM-dd');
            case 'time',          parse = @(X) datetime(char(X), 'InputFormat', 'HH:mm:ss');
            case 'base64Binary',  parse = @base64decode;
            otherwise,            parse = @char;
        end
    else
        type  = strsplit(type, ':');
        type  = type{end};
        parse = [];
    end