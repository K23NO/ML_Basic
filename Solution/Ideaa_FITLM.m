%LLamamos a la data en csv a matlab
Pre_Arahuay = readtable("Tabla 1 - Arahuay.csv", 'Delimiter', ';', 'VariableNamingRule', 'preserve');
Pre_Authisha = readtable("Tabla 2 Autisha.csv", 'Delimiter', ';', 'VariableNamingRule', 'preserve');
Pre_Cancha = readtable("Tabla 3 Canchacalla.csv", 'Delimiter', ';', 'VariableNamingRule', 'preserve');
Pre_Carampoma = readtable("Tabla 4 Carampoma.csv", 'Delimiter', ';', 'VariableNamingRule', 'preserve');
Pre_Matuca = readtable("Tabla 5 Matucana.csv", 'Delimiter', ';', 'VariableNamingRule', 'preserve');
Pre_RioBlanco = readtable("Tabla 6 Rio Blanco.csv", 'Delimiter', ';', 'VariableNamingRule', 'preserve');
Pre_SantaEula = readtable("Tabla 7 Santa EULAlia.csv", 'Delimiter', ';', 'VariableNamingRule', 'preserve');
Pre_SantTuna = readtable("Tabla 8 Santiago (DE TUNA).csv", 'Delimiter', ';', 'VariableNamingRule', 'preserve');
Pre_Nana = readtable("Tabla 9 Ñaña.csv", 'Delimiter', ';', 'VariableNamingRule', 'preserve');
%Pre_sheque = readtable("Tabla 10 Sheque.csv", 'Delimiter', ';', 'VariableNamingRule', 'preserve');

% Transformar cada tabla
Pre_Arahuay = transform_table(Pre_Arahuay);
Pre_Authisha = transform_table(Pre_Authisha);
Pre_Cancha = transform_table(Pre_Cancha);
Pre_Carampoma = transform_table(Pre_Carampoma);
Pre_Matuca = transform_table(Pre_Matuca);
Pre_RioBlanco = transform_table(Pre_RioBlanco);
Pre_SantaEula = transform_table(Pre_SantaEula);
Pre_SantTuna = transform_table(Pre_SantTuna);
Pre_Nana = transform_table(Pre_Nana);



% Unir las tablas en una sola
all_data = join(Pre_Arahuay, Pre_Authisha, 'Keys', {'Year', 'Month'});
all_data = join(all_data, Pre_Cancha, 'Keys', {'Year', 'Month'});
all_data = join(all_data, Pre_Carampoma, 'Keys', {'Year', 'Month'});
all_data = join(all_data, Pre_Matuca, 'Keys', {'Year', 'Month'});
all_data = join(all_data, Pre_RioBlanco, 'Keys', {'Year', 'Month'});
all_data = join(all_data, Pre_SantaEula, 'Keys', {'Year', 'Month'});
all_data = join(all_data, Pre_SantTuna, 'Keys', {'Year', 'Month'});
all_data = join(all_data, Pre_Nana, 'Keys', {'Year', 'Month'});

% Cambiar los nombres de las columnas para reflejar las estaciones
all_data.Properties.VariableNames = {'Year', 'Month', 'Arahuay', 'Authisha', 'Cancha', 'Carampoma', 'Matuca', 'RioBlanco', 'SantaEula', 'SantTuna', 'Nana'};



% Leer y transformar la tabla de Chosica
Pre_Chosica = readtable("PrecCHOSICA.csv", 'Delimiter', ';', 'VariableNamingRule', 'preserve');
Pre_Chosica = transform_table(Pre_Chosica);

% Unir la precipitación de Chosica con el resto de los datos
all_data = join(all_data, Pre_Chosica, 'Keys', {'Year', 'Month'});

% Cambiar el nombre de la columna de Chosica
all_data.Properties.VariableNames{end} = 'Chosica';


% Definir las variables independientes (X) y dependiente (y) -> AQUI
% COMIENZA EL USO DE LA FUNCIÓN FILTM
X = all_data{:, {'Arahuay', 'Authisha', 'Cancha', 'Carampoma', 'Matuca', 'RioBlanco', 'SantaEula', 'SantTuna', 'Nana'}};
y = all_data.Chosica;

% Ajustar el modelo de regresión lineal múltiple
mdl = fitlm(X, y);

% Mostrar los resultados del modelo
disp(mdl);

% Evaluar el modelo
R2 = mdl.Rsquared.Ordinary;
disp(['R^2: ', num2str(R2)]);

%Muestras de la función mld
anova(mdl,'summary')

plot(mdl)



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
