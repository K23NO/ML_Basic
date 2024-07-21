%LLamamos a la data en csv a matlab
Pre_Arahuay = readtable("Tabla 1 - Arahuay.csv", 'Delimiter', ';', 'VariableNamingRule', 'preserve');
Pre_Authisha = readtable("Tabla 2 Autisha.csv", 'Delimiter', ';', 'VariableNamingRule', 'preserve');
Pre_Canchacalla = readtable("Tabla 3 Canchacalla.csv", 'Delimiter', ';', 'VariableNamingRule', 'preserve');
Pre_Carampoma = readtable("Tabla 4 Carampoma.csv", 'Delimiter', ';', 'VariableNamingRule', 'preserve');
Pre_Matucana = readtable("Tabla 5 Matucana.csv", 'Delimiter', ';', 'VariableNamingRule', 'preserve');
Pre_RioBlanco = readtable("Tabla 6 Rio Blanco.csv", 'Delimiter', ';', 'VariableNamingRule', 'preserve');
Pre_SantaEula = readtable("Tabla 7 Santa EULAlia.csv", 'Delimiter', ';', 'VariableNamingRule', 'preserve');
Pre_SantTuna = readtable("Tabla 8 Santiago (DE TUNA).csv", 'Delimiter', ';', 'VariableNamingRule', 'preserve');
Pre_Nana = readtable("Tabla 9 Ñaña.csv", 'Delimiter', ';', 'VariableNamingRule', 'preserve');
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


% Unir las tablas en una sola
all_data = join(Pre_Arahuay, Pre_Authisha, 'Keys', {'Year', 'Month'});
all_data = join(all_data, Pre_Canchacalla, 'Keys', {'Year', 'Month'});
all_data = join(all_data, Pre_Carampoma, 'Keys', {'Year', 'Month'});
all_data = join(all_data, Pre_Matucana, 'Keys', {'Year', 'Month'});
all_data = join(all_data, Pre_RioBlanco, 'Keys', {'Year', 'Month'});
all_data = join(all_data, Pre_SantaEula, 'Keys', {'Year', 'Month'});
all_data = join(all_data, Pre_SantTuna, 'Keys', {'Year', 'Month'});
all_data = join(all_data, Pre_Nana, 'Keys', {'Year', 'Month'});


% Cambiar los nombres de las columnas para reflejar las estaciones
all_data.Properties.VariableNames = {'Year', 'Month', 'Arahuay', 'Authisha', 'Canchacalla', 'Carampoma', 'Matucana', 'RioBlanco', 'SantaEula', 'SantTuna', 'Nana'};

% Leer y transformar la tabla de Chosica
Pre_Chosica = readtable("PrecCHOSICA.csv", 'Delimiter', ';', 'VariableNamingRule', 'preserve');
Pre_Chosica = transform_table(Pre_Chosica);

% Unir la precipitación de Chosica con el resto de los datos
all_data = join(all_data, Pre_Chosica, 'Keys', {'Year', 'Month'});

% Cambiar el nombre de la columna de Chosica
all_data.Properties.VariableNames{end} = 'Chosica';

% Definir las variables independientes (X) y la variable dependiente (y)
X = [ones(size(all_data, 1), 1), all_data{:, {'Arahuay', 'Authisha', 'Canchacalla', 'Carampoma', 'Matucana', 'RioBlanco', 'SantaEula', 'SantTuna', 'Nana'}}];
y = all_data.Chosica;

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

% station_names = {'Intercepto', 'Arahuay', 'Authisha', 'Canchacalla', 'Carampoma', 'Matucana', 'RioBlanco', 'SantaEula', 'SantTuna', 'Nana'};
% disp('Coeficientes de regresión:');
% for i = 1:length(beta)
%     disp([station_names{i}, ': ', num2str(beta(i))]);
% end

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


%Grafica de los datos
% figure;
% plot(y, 'b', 'DisplayName', 'Real');
% hold on;
% plot(y_pred, 'r', 'DisplayName', 'Predicho');
% hold off;
% xlabel('Mes');
% ylabel('Precipitación');
% title('Precipitación Real vs. Predicha en Chosica');
% legend;

Ano = all_data.Year;
Mes = all_data.Month;
%grafica de estos datos, pero en el eje x que muestre el año en el eje x y la precipitación en el eje y
figure;
plot(y, 'b', 'DisplayName', 'Real');
hold on;
plot(y_pred, 'r', 'DisplayName', 'Predicho');
hold off;
set(gca, 'XTick', 1:12:length(y), 'XTickLabel', Ano(1:12:length(y)));
xlabel('Año');
ylabel('Precipitación total mensual (mm)');
title('Precipitación Real vs. Predicha en Chosica');
grid on;
legend;




%Grafica para mostar los datos, pero que ne el eje x, se muestre el año y el mes
% figure;
% plot(y, 'b', 'DisplayName', 'Real');
% hold on;
% plot(y_pred, 'r', 'DisplayName', 'Predicho');
% hold off;
% % solo mostrar los meses abajo, no el año
% set(gca, 'XTick', 1:12:length(y), 'XTickLabel', Mes(1:12:length(y)));
% xlabel('Mes');
% ylabel('Precipitación');
% title('Precipitación Real vs. Predicha en Chosica');
% legend;

% Analisis de errores MSE
MSE = sum((y - y_pred).^2) / n;
disp(['Error cuadrático medio (MSE): ', num2str(MSE)]);
%Normalizado: NMSE

%Error_cuadratico_medio = sum((y - y_pred).^2) / n;

% Graficar los datos reales vs. los datos predichos
figure;
scatter(1:length(y), y, 'filled', 'DisplayName', 'Datos Reales'); % Datos reales
hold on;
plot(1:length(y_pred), y_pred, 'r-', 'LineWidth', 2, 'DisplayName', 'Predicción del Modelo'); % Datos predichos
xlabel('Índice de Tiempo');
ylabel('Precipitación (mm)');
title('Comparación de la Precipitación Real y Predicha en Chosica');
legend('show');
grid on;



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
