function [obj,s] = ismrmrd_xml(fname)
% FORMAT obj = ismrmrd_xml(fname)
%
% Read the flexible header of an ismrmrd file into a structure

xmlstring = h5read(fname, '/dataset/xml');
xmlstring = xmlstring{1};

s = xmlreadstring(xmlstring);

obj = parseNode(s);

function node = parseNode(theNode)
% Recurse over node children.

node = struct;
if theNode.hasChildNodes
    childNodes    = theNode.getChildNodes;
    numChildNodes = childNodes.getLength;
    for c=1:numChildNodes
        theChild  = childNodes.item(c-1);
        childName = char(theChild.getNodeName);
        if childName(1) == '#'
            if ~isempty(char(theChild.getData))
                switch childName(2:end)
                    case 'text'
                        node.get_data = char(theChild.getData);
                    otherwise
                        warning('Unknown type %s',childName)
                end
            end
            continue
        end
        if ~isfield(node, childName)
            idx = 1;
        else
            idx = numel(node.(childName)) + 1;
        end
        childStruct = parseNode(theChild);
        fields = fieldnames(childStruct);
        for f=1:numel(fields)
            field = fields{f};
            node.(childName)(idx).(field) = childStruct.(field);
        end
    end
end
node.get_attributes = parseAttributes(theNode);


function attributes = parseAttributes(theNode)
% Create attributes structure.

attributes = containers.Map;
if theNode.hasAttributes
    attributeNodes = theNode.getAttributes;
    numAttributes  = attributeNodes.getLength;
    for c=1:numAttributes
        theAttribute   = attributeNodes.item(c-1);
        attributeName  = char(theAttribute.getName);
        attributeValue = char(theAttribute.getValue);
        attributes(attributeName) = attributeValue;
    end
end