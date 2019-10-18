function DenoisedImg= NAMF(x, Ds, ds, B)

% x is the image to be processed, Ds and ds should be set as 2 and 20, and B should be set as 0.8,
% respectively, DenoisedImg is the restored image.

% This code is for our paper "NAMF: A NON-LOCAL ADAPTIVE MEAN FILTER FOR SALT-AND-PEPPER NOISE REMOVAL".
% You can download our paper on https://arxiv.org/abs/1910.07787.

% Author: Houwang Zhang, School of Automation, 
% China University of Geoscience, Wuhan, China. 
% Released Date: 2019.10.18
% If you have found any bugs, have any suggestions or problems, please contact me at
% Email: zhanghw@cug.edu.cn

    x = double(x);
    M = x;

    N = zeros(size(x));
    N(x~=0 & x~=255) = 1; 
    f = zeros(size(x));
    f(x==0 | x==255) = 1; 

    smax = 3;

    [xlen, ylen] =size(x);


    for i = 1:xlen 
        for j = 1:ylen

            if N(i,j) == 0 
                s = 1;
                g = N(max(i-s,1):min(i+s,xlen), max(j-s,1):min(j+s,ylen)); 
                num = sum(g(:));
                
                while num == 0 && s<smax
                    s = s+1;
                    g = N(max(i-s,1):min(i+s,xlen), max(j-s,1):min(j+s,ylen));
                    num = sum(g(:));
                end

                if s<=smax && sum(g(:)>0) 
                    
                   tmp = x(max(i-s,1):min(i+s,xlen), max(j-s,1):min(j+s,ylen));                  
                   Ws = 1 - f(max(i-s,1):min(i+s,xlen), max(j-s,1):min(j+s,ylen));
                   M(i,j) = sum(sum(tmp.*Ws)) / sum(sum(Ws));

                else
                   
                   tmp = x(max(i-s,1):min(i+s,xlen), max(j-s,1):min(j+s,ylen));
                   pd = length(find(tmp == x(i,j))) / length(tmp(:));
                   
                   if pd > B
                       f(i, j) = 0;
                       continue
                   end
                    
                   up = max(i - 1, 1); left = max(j - 1, 1);
                   M(i,j) = median([M(up, left), M(up, j), M(i, left)]);
                   
                end

            end
        end
    end

    y = M;


    noise = sum(sum(f))/(xlen*ylen);
    hs = 4.5595 + 6.0314*noise + 2.2186*(noise^2);


    I = y;
    [m,n]=size(I);
    PaddedImg = padarray(I,[Ds+ds+1,Ds+ds+1],'symmetric','both');
    PaddedV = padarray(I,[Ds,Ds],'symmetric','both');
    average=zeros(m,n);
    wmax=average;
    sweight=average;
    h2=hs*hs;
    d=(2*ds+1)^2;
    for t1=-Ds:Ds
        for t2=-Ds:Ds
            if(t1==0&&t2==0)
                continue;
            end
            Sd=integralImgSqDiff(PaddedImg,Ds,t1,t2);
            SqDist2=Sd(2*ds+2:end-1,2*ds+2:end-1)+Sd(1:end-2*ds-2,1:end-2*ds-2)...
                -Sd(2*ds+2:end-1,1:end-2*ds-2)-Sd(1:end-2*ds-2,2*ds+2:end-1);
            SqDist2=SqDist2/d;
            w=exp(-SqDist2/h2);
            v = PaddedV(1+Ds+t1:end-Ds+t1,1+Ds+t2:end-Ds+t2);
            average=average+w.*v;
            wmax=max(wmax,w);
            sweight=sweight+w;
        end
    end
    average=average+0.*I;
    average=average./(0+sweight);
    DenoisedImg = average;

    for i = 1:m
        for j = 1:n
            if f(i, j) == 0
                DenoisedImg(i, j) = y(i, j);
            end
        end
    end
    
end

function Sd = integralImgSqDiff(PaddedImg,Ds,t1,t2)
Dist2=(PaddedImg(1+Ds:end-Ds,1+Ds:end-Ds)-PaddedImg(1+Ds+t1:end-Ds+t1,1+Ds+t2:end-Ds+t2)).^2;
Sd = cumsum(Dist2,1);
Sd = cumsum(Sd,2);
end
