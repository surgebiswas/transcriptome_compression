classdef (Abstract) GeneSearcher < handle
    %
    
    properties
        
    end
    
    methods
        function obj = GeneSearcher(obj)
        end
        
    end
    
    methods (Abstract)
        function gene_idx = find_next_gene(obj)
        end
    end
    
end

