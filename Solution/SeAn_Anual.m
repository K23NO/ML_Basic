% Leer las tablas de datos
Pre_Arahuay = readtable("Tabla 1 - Arahuay.csv", 'Delimiter', ';', 'VariableNamingRule', 'preserve');
Pre_Authisha = readtable("Tabla 2 Autisha.csv", 'Delimiter', ';', 'VariableNamingRule', 'preserve');
Pre_Canchacalla = readtable("Tabla 3 Canchacalla.csv", 'Delimiter', ';', 'VariableNamingRule', 'preserve');
Pre_Carampoma = readtable("Tabla 4 Carampoma.csv", 'Delimiter', ';', 'VariableNamingRule', 'preserve');
Pre_Matucana = readtable("Tabla 5 Matucana.csv", 'Delimiter', ';', 'VariableNamingRule', 'preserve');
Pre_RioBlanco = readtable("Tabla 6 Rio Blanco.csv", 'Delimiter', ';', 'VariableNamingRule', 'preserve');
Pre_SantaEula = readtable("Tabla 7 Santa EULAlia.csv", 'Delimiter', ';', 'VariableNamingRule', 'preserve');
Pre_SantTuna = readtable("Tabla 8 Santiago (DE TUNA).csv", 'Delimiter', ';', 'VariableNamingRule', 'preserve');
Pre_Nana = readtable("Tabla 9 Ñaña.csv", 'Delimiter', ';', 'VariableNamingRule', 'preserve');

% Transformar cada tabla a datos anuales
Pre_Arahuay = transform_table_annual(Pre_Arahuay);
Pre_Authisha = transform_table_annual(Pre_Authisha);
Pre_Canchacalla = transform_table_annual(Pre_Canchacalla);
Pre_Carampoma = transform_table_annual(Pre_Carampoma);
Pre_Matucana = transform_table_annual(Pre_Matucana);
Pre_RioBlanco = transform_table_annual(Pre_RioBlanco);
Pre_SantaEula = transform_table_annual(Pre_SantaEula);
Pre_SantTuna = transform_table_annual(Pre_SantTuna);
Pre_Nana = transform_table_annual(Pre_Nana);

% Unir las tablas en una sola
all_data = join(Pre_Arahuay, Pre_Authisha, 'Keys', {'Year'});
all_data = join(all_data, Pre_Canchacalla, 'Keys', {'Year'});
all_data = join(all_data, Pre_Carampoma, 'Keys', {'Year'});
all_data = join(all_data, Pre_Matucana, 'Keys', {'Year'});
all_data = join(all_data, Pre_RioBlanco, 'Keys', {'Year'});
all_data = join(all_data, Pre_SantaEula, 'Keys', {'Year'});
all_data = join(all_data, Pre_SantTuna, 'Keys', {'Year'});
all_data = join(all_data, Pre_Nana, 'Keys', {'Year'});

% Cambiar los nombres de las columnas para reflejar las estaciones
all_data.Properties.VariableNames = {'Year', 'Arahuay', 'Authisha', 'Canchacalla', 'Carampoma', 'Matucana', 'RioBlanco', 'SantaEula', 'SantTuna', 'Nana'};

% Leer y transformar la tabla de Chosica
Pre_Chosica = readtable("PrecCHOSICA.csv", 'Delimiter', ';', 'VariableNamingRule', 'preserve');
Pre_Chosica = transform_table_annual(Pre_Chosica);
% Unir la precipitación de Chosica con el resto de los datos
all_data = join(all_data, Pre_Chosica, 'Keys', {'Year'});

% Cambiar el nombre de la columna de Chosica
all_data.Properties.VariableNames{end} = 'Chosica';

% Definir las variables independientes (X) y la variable dependiente (y)
X = [ones(size(all_data, 1), 1), all_data{:, {'Arahuay', 'Authisha', 'Canchacalla', 'Carampoma', 'Matucana', 'RioBlanco', 'SantaEula', 'SantTuna', 'Nana'}}];
y = all_data.Chosica;

% Calcular los coeficientes de regresión
beta = (X' * X) \ (X' * y);

% Calcular los valores predichos
y_pred = X * beta;

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

% Graficar los datos
figure;
plot(y, 'b', 'DisplayName', 'Real');
hold on;
plot(y_pred, 'r', 'DisplayName', 'Predicho');
hold off;
set(gca, 'XTick', 1:length(y), 'XTickLabel', all_data.Year);
xlabel('Año');
ylabel('Precipitación media anual (mm)');
title('Precipitación Real vs. Predicha en Chosica');
grid on;
legend;

% Función para transformar la tabla a datos anuales
function data = transform_table_annual(T)
    % Convertir la tabla en un array para facilitar la manipulación
    data_array = table2array(T(:, 2:end-1));  % Ignorar la columna de Año y Total Anual
    years = T{:, 1};  % Obtener los años
    
    % Sumar las precipitaciones mensuales para cada año
    annual_precipitation = sum(data_array, 2)/12;
    
    % Crear una nueva tabla con datos anuales
    data = table(years, annual_precipitation, 'VariableNames', {'Year', T.Properties.VariableNames{2}});
end