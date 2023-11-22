function interleaver = HELICAL_interleaver(info_len)

    block = fix(sqrt(info_len)) + 1;
    temp = [1:info_len -ones(1,block*(block+1)-info_len)]; %不能用方正矩陣
    temp = reshape(temp,[],block)';
    
    interleaver = zeros(1,numel(temp));
    i1 = size(temp,1);
    i2 = 1;
    for i = 1:numel(temp)
        interleaver(i) = temp(i1,i2);
        if i1 - 1 >= 1 && i2 + 1 <= size(temp,2) %往右上
            i1 = i1 - 1;
            i2 = i2 + 1;
    
        elseif i1 - 1 >= 1 && i2 - size(temp,2) + 1 >= 1 %往最左、往上
            i1 = i1 -1;
            i2 = i2 - size(temp,2) + 1;
        else %往右、往最下
            i1 = i1 + size(temp,1) - 1;
            i2 = i2 + 1;
        end
    end

    temp_index = interleaver == -1;
    interleaver(temp_index) = [];

end