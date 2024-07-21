
% Cargar las tablas
Pre_Arahuay = readtable("Tabla 1 - Arahuay.csv", 'Delimiter', ';', 'VariableNamingRule', 'preserve');
Pre_Authisha = readtable("Tabla 2 Autisha.csv", 'Delimiter', ';', 'VariableNamingRule', 'preserve');
Pre_Canchacalla = readtable("Tabla 3 Canchacalla.csv", 'Delimiter', ';', 'VariableNamingRule', 'preserve');
Pre_Carampoma = readtable("Tabla 4 Carampoma.csv", 'Delimiter', ';', 'VariableNamingRule', 'preserve');
Pre_Matucana = readtable("Tabla 5 Matucana.csv", 'Delimiter', ';', 'VariableNamingRule', 'preserve');
Pre_RioBlanco = readtable("Tabla 6 Rio Blanco.csv", 'Delimiter', ';', 'VariableNamingRule', 'preserve');
Pre_SantaEula = readtable("Tabla 7 Santa EULAlia.csv", 'Delimiter', ';', 'VariableNamingRule', 'preserve');
Pre_SantTuna = readtable("Tabla 8 Santiago (DE TUNA).csv", 'Delimiter', ';', 'VariableNamingRule', 'preserve');
Pre_Nana = readtable("Tabla 9 Ñaña.csv", 'Delimiter', ';', 'VariableNamingRule', 'preserve');
Pre_Chosica = readtable("PrecCHOSICA.csv", 'Delimiter', ';', 'VariableNamingRule', 'preserve');

% Transformar cada tabla
Pre_Arahuay = transform_table(Pre_Arahuay);
Pre_Authisha = transform_table(Pre_Authisha);
Pre_Canchacalla = transform_table(Pre_Canchacalla);
Pre_Carampoma = transform_table(Pre_Carampoma);
Pre_Matucana = transform_table(Pre_Matucana);
Pre_RioBlanco = transform_table(Pre_RioBlanco);
Pre_SantaEula = transform_table(Pre_SantaEula);
Pre_SantTuna = transform_table(Pre_SantTuna);
Pre_Nana = transform_table(Pre_Nana);
Pre_Chosica = transform_table(Pre_Chosica);

% Unir las tablas en una sola
all_data = join(Pre_Arahuay, Pre_Authisha, 'Keys', {'Year', 'Month'});
all_data = join(all_data, Pre_Canchacalla, 'Keys', {'Year', 'Month'});
all_data = join(all_data, Pre_Carampoma, 'Keys', {'Year', 'Month'});
all_data = join(all_data, Pre_Matucana, 'Keys', {'Year', 'Month'});
all_data = join(all_data, Pre_RioBlanco, 'Keys', {'Year', 'Month'});
all_data = join(all_data, Pre_SantaEula, 'Keys', {'Year', 'Month'});
all_data = join(all_data, Pre_SantTuna, 'Keys', {'Year', 'Month'});
all_data = join(all_data, Pre_Nana, 'Keys', {'Year', 'Month'});
all_data = join(all_data, Pre_Chosica, 'Keys', {'Year', 'Month'});

% Cambiar los nombres de las columnas para reflejar las estaciones
all_data.Properties.VariableNames = {'Year', 'Month', 'Arahuay', 'Authisha', 'Canchacalla', 'Carampoma', 'Matucana', 'RioBlanco', 'SantaEula', 'SantTuna', 'Nana', 'Chosica'};

% Definir las variables independientes (X) y la variable dependiente (y)
X = [ones(size(all_data, 1), 1), all_data{:, {'Arahuay', 'Authisha', 'Canchacalla', 'Carampoma', 'Matucana', 'RioBlanco', 'SantaEula', 'SantTuna', 'Nana'}}];
y = all_data.Chosica;

% Calcular los coeficientes de regresión antes de eliminar el outlier
beta = (X' * X) \ (X' * y);
y_pred = X * beta;

% Calcular el R^2 antes de eliminar el outlier
R2_before = calcular_R2(y, y_pred);

% Eliminar el outlier (febrero de 1976)
outlier_index = find(all_data.Year == 1976 & all_data.Month == "Feb");
all_data_cleaned = all_data;
all_data_cleaned(outlier_index, :) = [];

% Redefinir las variables independientes (X) y la variable dependiente (y) sin el outlier
X_cleaned = [ones(size(all_data_cleaned, 1), 1), all_data_cleaned{:, {'Arahuay', 'Authisha', 'Canchacalla', 'Carampoma', 'Matucana', 'RioBlanco', 'SantaEula', 'SantTuna', 'Nana'}}];
y_cleaned = all_data_cleaned.Chosica;

% Calcular los coeficientes de regresión después de eliminar el outlier
beta_cleaned = (X_cleaned' * X_cleaned) \ (X_cleaned' * y_cleaned);
y_pred_cleaned = X_cleaned * beta_cleaned;

% Calcular el R^2 después de eliminar el outlier
R2_after = calcular_R2(y_cleaned, y_pred_cleaned);

% Graficar los resultados antes y después de eliminar el outlier
figure;
subplot(2, 1, 1);
plot(y, 'b', 'DisplayName', 'Real');
hold on;
plot(y_pred, 'r', 'DisplayName', 'Predicho');
hold off;
set(gca, 'XTick', 1:12:length(y), 'XTickLabel', Ano(1:12:length(y)));
title(['Precipitación Real vs. Predicha en Chosica (Antes de eliminar el outlier), R^2 = ', num2str(R2_before)]);
xlabel('Índice de Tiempo');
ylabel('Precipitación');
grid on;
legend;

subplot(2, 1, 2);
plot(y_cleaned, 'b', 'DisplayName', 'Real');
hold on;
plot(y_pred_cleaned, 'r', 'DisplayName', 'Predicho');
hold off;
set(gca, 'XTick', 1:12:length(y), 'XTickLabel', Ano(1:12:length(y)));
title(['Precipitación Real vs. Predicha en Chosica (Después de eliminar el outlier), R^2 = ', num2str(R2_after)]);
xlabel('Índice de Tiempo');
ylabel('Precipitación');
grid on;
legend;

disp(['R^2 antes de eliminar el outlier: ', num2str(R2_before)]);
disp(['R^2 después de eliminar el outlier: ', num2str(R2_after)]);

% Función para calcular R^2
function R2 = calcular_R2(y_real, y_predicho)
SS_res = sum((y_real - y_predicho).^2);
SS_tot = sum((y_real - mean(y_real)).^2);
R2 = 1 - (SS_res / SS_tot);
end

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

