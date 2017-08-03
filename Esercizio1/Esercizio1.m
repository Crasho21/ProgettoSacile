clear;
clc;

%Definizione numero job e costi
J = [1 2 3 4];
cleanup = [2 4 3 4];
setup = [10000 4 5 3;
         2 10000 4 6;
         3 5 10000 3;
         3 2 1 10000];
startup = [2 2 5 6];

%Definizione combinazioni degli stati
for k = 1 : length(J)
    temp = combnk(1 : length(J), k);
    for i = 1 : length(temp(:, 1))
        for j = 1 : k
            index = k * (i - 1) + j;
            X{k}(index) = Stato;
            X{k}(index).statiProcessati = temp(i, :);
            X{k}(index).ultimoJob = temp(i, j);
        end
    end
end

%Calcolo numero di stati per ogni stadio
for k = length(J) : -1 : 1
    stati(k) = length(X{k});
end

%Inizializzazione matrice dei percorsi per ogni stadio
for k = 1 : length(J)
    for i = 1 : max(stati)
        percorsi{k, i} = zeros(1, length(J));
    end
end

%Passo N
%Inizializzazione costi ottimi con costi di cleanup
Go{length(J)} = cleanup;

%Passo k = N - 1, ... , 1
%Passo iterativo per il calcolo del percorso con costo minore
for k = length(J) - 1 : -1 : 1
    %Inizializzazione matrice dei costi allo stadio k
    G{k} = 10000 * ones(stati(k),stati(k + 1));
    for i = 1 : stati(k)
        for j = 1 : stati(k + 1)
            %Se lo stato considerato dello stadio k è contenuto in quello
            %considerato dello stadio k + 1
            if(ismember(X{k}(i).statiProcessati, X{k + 1}(j).statiProcessati))
                %Calcolo il job di differenza tra i 2 stati considerati
                controllo(i, j) = setdiff(X{k + 1}(j).statiProcessati, X{k}(i).statiProcessati);
                %Se il job di differenza corrisponde con il job eseguito
                %allo stato dello stadio successivo
                if(isequal(controllo(i, j), X{k + 1}(j).ultimoJob))
                    %Stampa di debug
                    %disp('Transizione possibile');
                    %Aggiorno la matrice dei costi
                    G{k}(i, j) = Go{k + 1}(j) + setup(X{k}(i).ultimoJob, controllo(i, j));
                end
            end
        end
        %Calcolo il costo minimo ed il percorso ottimo dallo stadio k allo
        %stadio N
        [Go{k}(i), percorsi{k + 1, i}(1, k + 1)] = min(G{k}(i, :));
    end
end

%Passo 0
%Calcolo i costi allo stadio 0, considerando anche i costi di startup
for i = 1 : length(J)
        G0(i) = Go{1}(i) + startup(X{1}(i).ultimoJob);
end

%Calcolo il costo minimo ed il percorso ottimo dallo stadio 0 allo stadio N
[Go0, percorsi{1, 1}(1, 1)] = min(G0);

%Stampo il costo ottimo
disp(['Costo ottimo = ', num2str(Go0)]);

%Creo il testo contenente il percorso ottimo
index = 1;
text = 'Percorso ottimo = ';
for k = 1 : length(J)
    index = percorsi{k, index}(1, k);
    text = [text num2str(X{k}(index).ultimoJob) ' '];
end

%Stampo il percorso ottimo
disp(text);