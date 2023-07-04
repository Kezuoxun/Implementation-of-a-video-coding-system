function [mp_error, rec_img] = prediction_error(min_d, frame_1, frame_2, D)
    [rows, cols] = size(frame_1);
    for m = 1:length(min_d) 
        for row = min_d(m,1):min_d(m,5)
            for col = min_d(m,2):min_d(m,6)
                if(row+min_d(m,3) > 0 && col+min_d(m,4) > 0 && row+min_d(m,3) < rows && col+min_d(m,4) < cols)
                    mp_error(row,col) = (frame_2(row,col)) - (frame_1(row+min_d(m,3), col+min_d(m,4))); 
                    rec_img(row, col) = frame_1(row+min_d(m,3), col+min_d(m,4));
                else
                    mp_error(row,col) = 255;
                    rec_img(row, col) = 255;
                end
            end
        end
    end
    mp_error(min_d(end,5)-1:rows,min_d(end,6)-1:cols) = 0; 
    rec_img(min_d(end,5):rows,min_d(end,6):cols) = frame_2(min_d(end,5)-D:rows-D,min_d(end,6)-D:cols-D);  
end