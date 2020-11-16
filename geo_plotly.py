import plotly.graph_objs as go
import numpy as np
import math


def rotation_body(radius_array, z_array, nt = 50):
    theta = np.linspace(0, 2*np.pi, nt)
    fi, r = np.meshgrid(theta, radius_array)
    _,  z = np.meshgrid(theta, z_array)
    x = r * np.cos(fi)
    y = r * np.sin(fi)
    z = z
    return x, y, z    


def geodetic_line():
    theta = np.linspace(0, 2*np.pi, 50)
    r = 50
    z = np.linspace(0, 50, 50)

    x = r * np.cos(theta)
    y = r * np.sin(theta)
    z = z
    return x, y, z


def cartesian_from_polar(r, fi, h):
    x = r * np.cos(fi)
    y = r * np.sin(fi)
    z = h
    return x, y, z 


# Начальные геометрические условия
z_data = [0, 15, 80, 118, 258, 308]
diameters_data = [91, 91, 150, 170, 236, 236]
r_data = [diam/2 for diam in diameters_data]

# Зеркально отображаем данные
z_reverse = [z_data[-1] * 2 - x for x in reversed(z_data[:-1])]
z_data = z_data + z_reverse
r_data = r_data + list(reversed(r_data[:-1]))

# Линейная интерполяция для получения num_points точек
num_points = z_data[-1] - z_data[0] + 1
z_poits = np.linspace(z_data[0], z_data[-1], num_points)
r_points = np.interp(z_poits, z_data, r_data)

# Точки для построения поверхности вращения
x1, y1, z1 = rotation_body(r_points, z_poits)

# Дифференциируем массив радиусов
diff_r_points = np.diff(r_points)
diff_r_points = np.append(diff_r_points, diff_r_points[-1])

# Находим коэффициент альфа, зависящий от угла армирования
# Угол армирования psi_0 задаём для конкретного радиуса r_0
psi_0 = 87
r_0 = r_points[0]
alfa = 1 / (r_0 * math.sin(psi_0 * np.pi / 180) )

# Находим dz
dz = (z_data[-1] - z_data[0] + 1) / num_points

# Вычисляем подыинтегральное выражение (65) со стр. 426 книги Бухгольца
# для каждой точки
f = r_points
dfdz = diff_r_points
geo = dz / f * np.sqrt((1 + dfdz * dfdz) / (alfa * alfa * f * f - 1))

# Находим углы как интеграл в каждой точке
fi = np.cumsum(geo)

# Точки для построения геодезической линии
x2, y2, z2 = cartesian_from_polar(r_points, fi, z_poits)

#Опрелеляем угол армировния в большем сечении через формулу Клеро r*sin = const
width = 15 # пока с потолка
r_center = r_points[-1]
psi_center = math.asin(r_0 / r_center * math.sin(psi_0 * np.pi / 180))
width_center = width / math.sin(psi_center)
n = (2 * r_center * np.pi) / width_center

# Добавляем второй виток
fi2 = fi + fi[-1]
r_points = np.concatenate((r_points, r_points))
fi = np.concatenate((fi, fi2))
z_poits = np.concatenate((z_poits, np.flip(z_poits)))


x, y, z = cartesian_from_polar(r_points, fi, z_poits)
curve = go.Scatter3d(x=z, y=y, z=x, mode='lines', line={'width': 10})


'''curves = []
for i in range(19):
    x, y, z = cartesian_from_polar(r_points, fi + i * 2 * np.pi / 19, z_poits)
    curves.append(go.Scatter3d(x=z, y=y, z=x, mode='lines', line={'width': 10}))'''


# Рисуем графики
surface = go.Surface(x=z1, y=y1, z=x1, colorscale='Viridis', opacity=0.5)
#curve = go.Scatter3d(x=z2, y=y2, z=x2, mode='lines', line={'width': 10})
data = [surface, curve]
#data = [surface] + curves
fig = go.Figure(data=data)
fig.update_layout(
    title={
        'text': "Угол армирования " + str(psi_0),
        'y':0.9,
        'x':0.5,
        'xanchor': 'center',
        'yanchor': 'top'},
    #scene={'bgcolor': 'cyan'},
    plot_bgcolor='cyan')
fig.write_html('tmp.html', auto_open=True)