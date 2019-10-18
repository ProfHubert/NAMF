close all; 
clear, clc;

% load('img256/cameraman.mat');
% I = cameraman;

I = imread('data/img512/Barbara.bmp');

degree = 0.1:0.1:0.9;

JJ = I;

for i = 1:9
    
    I = imnoise(JJ, 'salt & pepper', degree(i));
    
    O = NAMF(I, 2, 20, 0.8);
    
    figure;
    
    imshow([JJ, I, O], []); 

    title(['density=',num2str(degree(i))])
    
end

