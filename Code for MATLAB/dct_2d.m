%  https://youtu.be/ensCM2hIFH0

%% DCT transform with 4 for loops
function DCT = dct_2d(image)


image = double(image);
N= length(image);
DCT=zeros(N);
% cu = sqrt(2/N);
% cv = sqrt(2/N);
cu = sqrt(1/2);
cv = sqrt(1/2);
A = ones(1,N);A(1) = 1/sqrt(2);

for u=0:N-1
    for v=0:N-1
        mysum = 0;
        for x=0:N-1
            for y=0:N-1
                mysum = mysum + A(x+1)*A(y+1) * cos(pi*u*(2*x+1)/(2*N))*cos(pi*v*(2*y+1)/(2*N))*image(x+1,y+1);

             end
        end
        DCT(u+1,v+1) = cu*cv*mysum;
    end
    % presenting the recent number of step from the whole steps
    disp([num2str(u)  ' of ' num2str(N-1)])
end
end



