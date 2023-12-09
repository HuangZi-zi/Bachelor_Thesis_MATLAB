classdef MinHeap
    properties
        heapArray
        heapSize
    end
    
    methods
        function obj = MinHeap(array)
            if nargin == 1
                obj.heapArray = array;
                obj.heapSize = numel(array);
                obj.buildHeap();
            else
                obj.heapArray = [];
                obj.heapSize = 0;
            end
        end
        
        function buildHeap(obj)
            for i = floor(obj.heapSize / 2):-1:1
                obj.minHeapify(i);
            end
        end
        
        function minHeapify(obj, index)
            left = 2 * index;
            right = 2 * index + 1;
            smallest = index;
            
            if left <= obj.heapSize && obj.heapArray(left) < obj.heapArray(smallest)
                smallest = left;
            end
            
            if right <= obj.heapSize && obj.heapArray(right) < obj.heapArray(smallest)
                smallest = right;
            end
            
            if smallest ~= index
                % Swap values
                temp = obj.heapArray(index);
                obj.heapArray(index) = obj.heapArray(smallest);
                obj.heapArray(smallest) = temp;
                
                % Recursively heapify the affected sub-tree
                obj.minHeapify(smallest);
            end
        end
        
        function insert(obj, value)
            obj.heapSize = obj.heapSize + 1;
            obj.heapArray(obj.heapSize) = value;
            
            currentIndex = obj.heapSize;
            
            % Maintain the heap property
            while currentIndex > 1 && obj.heapArray(floor(currentIndex / 2)) > obj.heapArray(currentIndex)
                % Swap values
                temp = obj.heapArray(floor(currentIndex / 2));
                obj.heapArray(floor(currentIndex / 2)) = obj.heapArray(currentIndex);
                obj.heapArray(currentIndex) = temp;
                
                currentIndex = floor(currentIndex / 2);
            end
        end
        
        function minElement = extractMin(obj)
            if obj.heapSize < 1
                error('Heap underflow');
            end
            
            minElement = obj.heapArray(1);
            obj.heapArray(1) = obj.heapArray(obj.heapSize);
            obj.heapSize = obj.heapSize - 1;
            
            obj.minHeapify(1);
        end
    end
end
