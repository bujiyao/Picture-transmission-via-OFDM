function interleaver = Matrix_interleaver(info_len)
    block = fix(sqrt(info_len)) + 1;
    temp = [1:info_len -ones(1,block^2-info_len)]; 
    temp = reshape(temp,[],block)';
    temp2 = [temp;temp];

    interleaver = zeros(1,numel(temp));
    i1 = 1;
    i2 = 1;
    i_itv = 1;
    for i = 1:max(size(temp))
        for j = 0:min(size(temp))-1
            interleaver(i_itv) = temp2(i1+j,i2+j);
            i_itv = i_itv + 1;
        end
        i1 = i1+1;
    end

    temp_index = interleaver == -1;
    interleaver(temp_index) = [];
end