function [ path ] = getPath( )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

 if isdir('C:\Users\noambox\Documents\Sync')
     path = 'C:\Users\noambox\Documents\Sync';
 elseif isdir('C:\Users\Noam\Documents\Sync')
     path = ('C:\Users\Noam\Documents\Sync');
 else
     error('getPath error :: cumputer path is not recognized...');
 end
end

