

% Definición de los archivos de datos
files = ["Tabla 3 Canchacalla.csv", ...
         "Tabla 7 Santa EULAlia.csv", "Tabla 8 Santiago (DE TUNA).csv", "Tabla 9 Ñaña.csv"];

% Cargar y transformar cada tabla
for i = 1:length(files)
    % Leer cada archivo con las opciones de formato específico
    T = readtable(files(i), 'Delimiter', ';', 'VariableNamingRule', 'preserve');
    % Transformar la tabla
    transformedTables{i} = transform_table(T);
end

% Unir las tablas transformadas en una sola tabla 'all_data'
all_data = transformedTables{1};
for i = 2:length(transformedTables)
    all_data = join(all_data, transformedTables{i}, 'Keys', {'Year', 'Month'});
end


% Cambiar los nombres de las columnas para reflejar las estaciones
all_data.Properties.VariableNames = {'Year', 'Month', 'Canchacalla', 'SantaEula', 'SantTuna', 'Nana'};

% Cargar y transformar la tabla de Chosica
Pre_Chosica = readtable("PrecCHOSICA.csv", 'Delimiter', ';', 'VariableNamingRule', 'preserve');
Pre_Chosica = transform_table(Pre_Chosica);

% Analizar y eliminar datos atípicos
mean_chosica = mean(Pre_Chosica.Precipitation, 'omitnan');
std_chosica = std(Pre_Chosica.Precipitation, 'omitnan');
outlier_index = find(Pre_Chosica.Precipitation > mean_chosica + 2 * std_chosica | Pre_Chosica.Precipitation < mean_chosica - 2 * std_chosica);

% Unir la precipitación de Chosica con el resto de los datos
all_data = join(all_data, Pre_Chosica, 'Keys', {'Year', 'Month'});
all_data.Properties.VariableNames{end} = 'Chosica';

%Eliminar los datos atipicos de chosica en el all_data
all_data(outlier_index, :) = [];

% Calcular los promedios anuales para todas las estaciones, incluyendo Chosica
%annual_averages = calculate_annual_averages(all_data);
annual_averages = calculate_annual_averages(all_data);


%Definincion de mi matri
X = table2array(annual_averages(:, 1:end-2)); % Variables independientes (todas las estaciones excepto Chosica)
y = annual_averages.Chosica; % Variable dependiente (precipitación en Chosica)

% Añadir término de intercepción a X
X = [ones(size(X, 1), 1), X];
% Calcular los coeficientes de regresión
beta = (X' * X) \ (X' * y);


%evaular el modelo
% Calcular los valores predichos
y_pred = X * beta;

%Errores
% Calcular la varianza residual
n = length(y);
k = size(X, 2) - 1;
residuals = y - y_pred;
sigma2 = sum(residuals.^2) / (n - k - 1);

% Calcular la matriz de varianza-covarianza de los coeficientes
var_beta = sigma2 * inv(X' * X);

% Calcular los errores estándar
SE = sqrt(diag(var_beta));

% Calcular los estadísticos t
tStat = beta ./ SE;

% Calcular los valores p
pValue = 2 * (1 - tcdf(abs(tStat), n - k - 1));

% Calcular el R^2
SS_res = sum((y - y_pred).^2);
SS_tot = sum((y - mean(y)).^2);
R2 = 1 - (SS_res / SS_tot);

% Calcular el estadístico F
MS_reg = (SS_tot - SS_res) / k;
MS_res = SS_res / (n - k - 1);
F_stat = MS_reg / MS_res;
pValue_F = 1 - fcdf(F_stat, k, n - k - 1);

%la homegenidad de los resulatdos-> prueba estadística de cramer
% Calcular el estadístico de Cramer
%Cramer_stat = sqrt(R2 / (1 - R2))
% Mostrar los resultados con etiquetas


% Mostrar los resultados
disp('Coeficientes de regresión:');
disp(beta);
disp('Errores estándar (SE):');
disp(SE);
disp('Estadísticos t (tStat):');
disp(tStat);
disp('Valores p (pValue):');
disp(pValue);
disp(['R^2: ', num2str(R2)]);
disp(['Estadístico F: ', num2str(F_stat)]);
disp(['Valor p del estadístico F: ', num2str(pValue_F)]);



Ano = annual_averages.Year;
%grafica de estos datos, pero en el eje x que muestre el año en el eje x y la precipitación en el eje y
figure;
plot(y, 'b', 'DisplayName', 'Real');
hold on;
plot(y_pred, 'r', 'DisplayName', 'Predicho');
hold off;
%mostras mas años en el eje x
xticks(1:length(Ano));
xticklabels(Ano);

xlabel('Año');
ylabel('Precipitación media anual promedio (mm)');
title('Precipitación Real vs. Predicha en Chosica');
grid on;
legend;


% Función para transformar la tabla
function data = transform_table(T)
    % Convertir la tabla en un array para facilitar la manipulación
    data_array = table2array(T(:, 2:end-1));  % Ignorar la columna de Año y Total Anual
    years = T{:, 1};  % Obtener los años
    
    % Crear una matriz donde cada fila es un mes de un año específico
    months = ["Ene", "Feb", "Mar", "Abr", "May", "Jun", "Jul", "Ago", "Sep", "Oct", "Nov", "Dic"];
    num_years = size(data_array, 1);
    num_months = length(months);
    
    % Inicializar la tabla resultante
    data = table;
    
    for i = 1:num_years
        for j = 1:num_months
            new_row = table(years(i), months(j), data_array(i, j), 'VariableNames', {'Year', 'Month', 'Precipitation'});
            data = [data; new_row];
        end
    end
end

function annual_averages = calculate_annual_averages(data)
    % Extraer años únicos y estaciones (columnas desde la tercera hasta la última)
    unique_years = unique(data.Year);
    stations = data.Properties.VariableNames(3:end); % Incluye Chosica al final

    % Inicializar la tabla para guardar los promedios anuales
    annual_averages = array2table(zeros(length(unique_years), length(stations)), 'VariableNames', stations);
    annual_averages.Year = unique_years;

    % Calcular el promedio anual por estación y año
    for i = 1:length(unique_years)
        year_data = data(data.Year == unique_years(i), :);
        for j = 1:length(stations)
            % Suma total anual por estación y divide por el número de meses con datos
            valid_data = year_data{:, stations(j)};
            valid_data = valid_data(~isnan(valid_data)); % Eliminar NaNs
            if isempty(valid_data)
                annual_averages{i, j} = NaN; % Asignar NaN si no hay datos ese año
            else
                annual_averages{i, j} = mean(valid_data); % Guardar el promedio anual
            end
        end
    end
end

