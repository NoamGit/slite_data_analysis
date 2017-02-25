classdef Iterator < handle
    % iterator for cellarray and array
    
    properties
        ii % bookmark
        iterable
        numToGo % number of elemnets
        type
        exc = struct('message','not iterable data!','identifier','Iterator:notIterable');

    end
    
    methods
        function obj = Iterator(elem_list)
            if iscell(elem_list) || ismatrix(elem_list)
                elem_list = reshape(elem_list,1,[]);
                obj.numToGo = size(elem_list,2);
                if obj.numToGo > 0 && size(elem_list,1) == 1
                    if iscell(elem_list);obj.type = 'cell';end;
                    if ismatrix(elem_list);obj.type = 'array';end;
                    obj.ii = 1;
                    obj.iterable = elem_list;
                else
                    error(obj.exc);
                end
            else 
                 error(obj.exc);
            end
                
        end
        
        function bool = hasNext(obj)
            bool = obj.numToGo > 0;
        end
        
        function out = getNext(obj)
            if obj.hasNext()
                out = obj.get(obj.ii);
                obj.ii = obj.ii + 1;
                obj.numToGo = obj.numToGo - 1;
            else
                out = [];
            end
        end
    end
    
    methods(Access = private)
        
        function element = get(obj, index)
            switch obj.type
                case 'cell'
                    element = obj.iterable{index};
                case 'array'
                    element = obj.iterable(index);
            end
        end
    end      
end
    

