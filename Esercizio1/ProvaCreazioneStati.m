temp = combnk(1:4,1);
for i = 1 : 4
    X1(i) = Stato;
    X1(i).statiProcessati = temp(i, :);
    X1(i).ultimoJob = i;
end
temp = combnk(1:4,2);
for i = 1 : length(temp(:, :))
    for j = 1 : 2
        index = 2*(i-1)+j;
        X2(index) = Stato;
        X2(index).statiProcessati = temp(i, :);
        X2(index).ultimoJob = temp(i, j);
    end
end
temp = combnk(1:4,3);
for i = 1 : length(temp(:, :))
    for j = 1 : 3
        index = 3*(i-1)+j;
        X3(index) = Stato;
        X3(index).statiProcessati = temp(i, :);
        X3(index).ultimoJob = temp(i, j);
    end
end
temp = combnk(1:4,4);
for j = 1 : 4
        X4(j) = Stato;
        X4(j).statiProcessati = temp(1, :);
        X4(j).ultimoJob = temp(1, j);
end