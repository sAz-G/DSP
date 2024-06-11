clc
% Define the array to be sorted
lgth = 10000000;
t = linspace(0,1,lgth);
inpArr = rand(1,lgth); %[1,0,2,0,3,0,4,0,5,0,6,0,7,0,8,0,9,0,10,0];
inpArr_comp = rand(1,lgth);

final_input = zeros(1,2*lgth);
final_input(1:2:end-1) = inpArr;
final_input(2:2:end) = inpArr_comp;

inputArray = libpointer('singlePtr',final_input);
% Preallocate the output array
outputArray =  libpointer('singlePtr',final_input*0);

% Load the library
if ~libisloaded('myfft')
    loadlibrary('myfft', 'fft10.h');
end

% Call the C function
tic
calllib('myfft', 'local_fft_radix10', 0,1, 0, lgth, inputArray, outputArray);
toc 

% Display the sorted array
disp('Sorted Array:');
outputArray.Value;
% Unload the library
unloadlibrary('myfft');

tic
fft_res = fft(inpArr + 1j*inpArr_comp);
toc


final_output = outputArray.Value(1:2:end-1) + 1j*outputArray.Value(2:2:end);
max(abs(final_output) - abs(fft_res))