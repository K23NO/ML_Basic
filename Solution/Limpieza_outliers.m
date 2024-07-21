Pre_Chosica = readtable("PrecCHOSICA.csv", 'Delimiter', ';', 'VariableNamingRule', 'preserve');
Pre_Chosica = transform_table(Pre_Chosica);
% Unir la precipitación de Chosica con el resto de los datos
% all_data = join(Pre_Chosica, 'Keys', {'Year', 'Month'});


%has un scatter plot de la precipitacion de chosica
scatter(Pre_Chosica.Year, Pre_Chosica.Precipitation, 'filled', 'MarkerFaceColor', 'b');
xlabel('Año');
ylabel('Precipitación (mm)');
title('Precipitación en Chosica');
grid on;

% halla la media y la desviacion estandar de la precipitacion de chosica
mean_chosica = mean(Pre_Chosica.Precipitation);
std_chosica = std(Pre_Chosica.Precipitation);

% Encuentra los índices de los datos atípicos
outlier_index = find(Pre_Chosica.Precipitation > mean_chosica + 2 * std_chosica | Pre_Chosica.Precipitation < mean_chosica - 2 * std_chosica);

%encuentra los datos atipicos solo con la desviacion estandar
%oulier_index_std = find(Pre_Chosica.Precipitation > mean_chosica + 2 * std_chosica);
%muestra estos datos atipicos en la grafica
hold on;
scatter(Pre_Chosica.Year(outlier_index), Pre_Chosica.Precipitation(outlier_index), 'r', 'filled', 'MarkerFaceColor', 'r');
hold off;
legend('Datos', 'Datos atípicos');
%calcula el porcentaje de datos atipicos
porcentaje_outliers = length(outlier_index) / length(Pre_Chosica.Precipitation) * 100;
fprintf('Porcentaje de datos atípicos: %.2f%%\n', porcentaje_outliers);


%suma de cuantos datos atipicos hay
num_outliers = length(outlier_index);
fprintf('Número de datos atípicos: %d\n', num_outliers);


% Eliminar los datos atípicos
Pre_Chosica(outlier_index, :) = [];

% Guardar los datos limpios en un archivo CSV
writetable(Pre_Chosica, 'PrecCHOSICA_limpio.csv', 'Delimiter', ';');



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
