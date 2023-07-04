function [total_time, min_d, motion_vectors] = EBMA(frame_1, frame_2, rows, cols, blockSize, searchRange, error_fcn)
    start_time = cputime;    
    min_d = zeros(searchRange,1); 
    b=1; 
    G = imgradientxy(frame_1);
    %G = Gx+Gy;
    for row = 1:blockSize:rows 
        for col = 1:blockSize:cols
            error = zeros(blockSize,1);
            count = 0;
            for block_row = -searchRange:searchRange 
                for block_col = -searchRange:searchRange
                    if(row+block_row<1 || row+block_row+blockSize-1>rows || col+block_col<1 || col+block_col+blockSize-1>cols) % Boundary condition
                        continue;
                    end
                    count = count +1;
                    if count == 1
                        block_start_m = row+block_row;
                        block_start_n = col+block_col;
                    end          
                    if strcmp(error_fcn, 'dfd')
                        block1 = frame_1(row:row+blockSize-1, col:col+blockSize-1);
                        block2 = frame_2(row+block_row:row+blockSize-1+block_row, col+block_col:col+blockSize-1+block_col);  
                        error(row+block_row, col+block_col) = sum(sum(abs(int16(block1) - int16(block2))));
                    elseif strcmp(error_fcn, 'of')
                        grd = G(row+block_row:row+blockSize-1+block_row, col+block_col:col+blockSize-1+block_col);
                        %display(G(row:row+blockSize-1, col:col+blockSize-1))
                        %display([0:blockSize-1; 0:blockSize-1])
                        % grd = G(row:row+blockSize-1, col:col+blockSize-1)*[1:blockSize; 1:blockSize;]';
                        block1 = frame_1(row:row+blockSize-1, col:col+blockSize-1);
                        block2 = frame_2(row:row+blockSize-1, col:col+blockSize-1);
                        error(row+block_row, col+block_col) = sum(sum(((int16(grd) + int16(block1) - int16(block2)).^2)));
                    end
                end
            end
            ab1 = find(any(error,2),1,'last');
            ab2 = find(any(error,1),1,'last');
            error(~any(error,2), :) = [];
            error(:, ~any(error,1)) = [];
            b1 = find(any(error,2),1,'last');
            b2 = find(any(error,1),1,'last');
            [~, tmp2] = min(error(:));
            if col+searchRange >= cols-searchRange+1
               new_centre_x = b1 - searchRange;
               new_centre_y = b2;
            elseif row+searchRange >= rows-searchRange+1
               new_centre_x = b1;
               new_centre_y = b2 - searchRange;
            else
               new_centre_x = b1 - searchRange;
               new_centre_y = b2 - searchRange;
            end
            min_d(b,1) = block_start_m;
            min_d(b,2) = block_start_n;
            [min_d(b,3), min_d(b,4)] = ind2sub(size(error()), tmp2);
            min_d(b,3) = min_d(b,3) - new_centre_x;
            min_d(b,4)  = min_d(b,4) - new_centre_y;
            min_d(b,5) = ab1;
            min_d(b,6) = ab2;
            b = b+1;
        end
    end
    total_time = cputime - start_time;
    tmp = 1;
    for i = 2:length(min_d)
        motion_vectors(tmp,1:4) = min_d(i,1:4);
        motion_vectors(tmp+1,1) = min_d(i,1) + 4;
        motion_vectors(tmp+1,2) = min_d(i,2);
        motion_vectors(tmp+2,1) = min_d(i,1);
        motion_vectors(tmp+2,2) = min_d(i,2) + 4;
        motion_vectors(tmp+3,1:2) = min_d(i,1:2) + 4;
        motion_vectors(tmp+1:tmp+3,3:4) = repmat(min_d(i,3:4),3,1);
        tmp = tmp+4;
    end
%     fig = figure();
%     set(gcf, 'renderer', 'zbuffer')
%     imshow(imrotate(frame_1, 90), 'InitialMagnification',500,'Interpolation',"bilinear");
%     hold on;
%     quiver(motion_vectors(:,1), motion_vectors(:,2), motion_vectors(:,3), motion_vectors(:,4),0, 'r');
%     camroll(-90);
%     shading interp;
%     if strcmp(error_fcn, 'dfd')
%         title('Motion vector - DFD')
%     elseif strcmp(error_fcn, 'of')
%         title('Motion vector - OF')
%     end
end
