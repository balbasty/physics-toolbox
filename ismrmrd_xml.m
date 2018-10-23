function s = ismrmrd_xml(fname)
% FORMAT s = ismrmrd_xml(fname)
%
% Read the flexible header of an ismrmrd file into a structure

xmlstring = h5read(fname, '/dataset/xml');
xmlstring = xmlstring{1};

s = xmlreadstring(xmlstring);

