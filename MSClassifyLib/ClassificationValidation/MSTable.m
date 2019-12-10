classdef MSTable < matlab.mixin.Copyable
    % Container class for processing Tables
    % Verifications are performed during the classification study
    % Properties
    %   string: cell array of string identifiers
    %   object: processing object (can be a list for map modifiers and %
    %           preprocessings)
    %
    properties
        string; %cell array of string identifiers
        object; %processing object
    end
    methods
        function add(obj, string, object)
            % function adding a new processing element to the list
            % INPUT
            %   string: process identifier string
            %   object: processing object
            obj.string{end+1}=string;
            obj.object{end+1}=object;
        end        
    end
end