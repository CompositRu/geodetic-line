import plotly.graph_objs as go
import numpy as np
from math import sin, cos, asin, radians, degrees, tan, sqrt
import math


# -------------------------------------------------------------------------------------------
# ------- НАХОЖДЕНИЕ КООРДИНАТ ГЕОДЕЗИЧЕСКОЙ ЛИНИИ НА ПОВЕРХНОСТИ ТЕЛА ВРАЩЕНИЯ -------------
# -------------------------------------------------------------------------------------------


def get_shift_angle(np_fi, num_turns):
    central_angle = 2 * np.pi / num_turns
    delta_fi = (np_fi[-1] - np_fi[0]) % (2 * np.pi) 
    num_central_angles = int(delta_fi // central_angle) + 1
    # Для равномерного покрытия 
    while (math.gcd(num_turns, num_central_angles) != 1 or num_turns == num_central_angles):
        num_central_angles += 1
    shift_angle = central_angle * num_central_angles - delta_fi
    print('shift_angle rad:', shift_angle)
    print('shift_angle deg:', degrees(shift_angle))
    return shift_angle


def get_rare(ar):
    # Разряжаем массив опираясь на внутреннее чутьё
    def f(ar):
        idx1 = len(ar) // 1000
        idx2 = len(ar) // 100
        idx3 = len(ar) // 10
        ar1 = ar[:idx1]
        ar2 = ar[idx1:idx2:3]
        ar3 = ar[idx2:idx3:30]
        ar4 = ar[idx3::300]
        if type(ar) == 'numpy.ndarray':
            return np.concatenate((ar1, ar2, ar3, ar4))
        else:
            return ar1 + ar2 + ar3 + ar4
    def rev(ar):
        if type(ar) == 'numpy.ndarray':
            return np.flip(ar)
        else:
            return list(reversed(ar)) 
    idx_ar = list(range(len(ar)))        
    idx = f(idx_ar[:len(idx_ar)//2]) + rev(f(rev(idx_ar[len(idx_ar)//2:])))
    if type(ar) == 'numpy.ndarray':
        return numpy_ar[idx]
    else:
        res = []
        for i in idx:
            res.append(ar[i])
        return res


# Начальные геометрические условия
# start_z_data = [15, 70, 90, 118, 170, 244, 268, 308]
# start_diameters_data = [91, 141, 156, 172, 199, 230, 236, 236]

start_z_data = [15, 80, 118, 258, 308]
start_z_data = [x - 15 for x in start_z_data]
start_diameters_data = [91, 150, 170, 236, 236]
start_r_data = [diam / 2 for diam in start_diameters_data]
#z_data = [0, 15, 80, 118, 258, 308]
#diameters_data = [91, 91, 150, 170, 236, 236]
#z_data = [15, 80, 118, 258, 286, 290, 292, 294, 308]
#diameters_data = [91, 150, 170, 236, 106, 100, 100, 106, 236]

# Зеркально отображаем данные и дополняем стартовые массивы
z_reverse = [start_z_data[-1] * 2 - x for x in reversed(start_z_data[:-1])]
start_z_data = start_z_data + z_reverse
start_r_data = start_r_data + list(reversed(start_r_data[:-1]))

# Линейная интерполяция для получения num_points точек
num_points = (start_z_data[-1] - start_z_data[0] + 1) * 10
np_z = np.linspace(start_z_data[0], start_z_data[-1], num_points)
np_r = np.interp(np_z, start_z_data, start_r_data)

# Дифференциируем массив радиусов
diff_r_points = np.diff(np_r)
diff_r_points = np.append(diff_r_points, diff_r_points[-1])

# Находим коэффициент альфа, зависящий от угла армирования
# Угол армирования psi_0 задаём для конкретного радиуса r_0
psi_0 =89
r_0 = np_r[0]
alfa = 1 / (r_0 * math.sin(psi_0 * np.pi / 180) )

# Находим dz для интеграла
dz = (start_z_data[-1] - start_z_data[0]) / num_points

# Вычисляем подыинтегральное выражение (65) со стр. 426 книги Бухгольца
# для каждой точки
f = np_r
dfdz = diff_r_points
geo = dz / f * np.sqrt((1 + dfdz * dfdz) / (alfa * alfa * f * f - 1))

# Находим углы как интеграл в каждой точке
np_fi = np.cumsum(geo)

# Опрелеляем угол армировния в большем сечении через формулу Клеро 
# r * sin = const
r_center = max(np_r)
psi_center = asin(r_0 / r_center * sin(radians(psi_0)))

# Выбор количества витков    
num_turns = 50 #количество витков для полного покрытия наибольшего диаметра (взято с потолка)
width = 2 * np.pi * r_center / num_turns 
print('width:', width) 
width_center = width / cos(psi_center) # данное значение должно быть близко к реальному значению ширины жгута
print('width_center:', width_center)

# Находим смещение второго витка, чтобы все витки равномерно покрывали поверхность
shift_angle = get_shift_angle(np_fi, num_turns)

# Добавляем второй виток
np_r_1_turn = np_r
np_fi_1_turn = np_fi
np_z_1_turn = np_z

np_fi_2 = np_fi_1_turn + np_fi[-1] + shift_angle
np_r_2_turns = np.concatenate((np_r_1_turn, np_r_1_turn))
np_fi_2_turns = np.concatenate((np_fi_1_turn, np_fi_2))
np_z_2_turns = np.concatenate((np_z_1_turn, np.flip(np_z_1_turn)))

# Перед разряжение массивов данных создаём отдельно массив данных для будущего рассчёта векторов схода нити
# Объединяем r, fi, z в один массив данных и добавляем в каждую точку координаты следующей (чтобы эти данные не потерялись при расзряжении)
directions_1 = np.stack((np_r_2_turns, np_fi_2_turns, np_z_2_turns), axis=1) # либо axis=1, либо транспонируем полученную матрицу 
directions_2 = np.stack((np.roll(np_r_2_turns, -1), np.roll(np_fi_2_turns, -1), np.roll(np_z_2_turns, -1)), axis=1)
directions = np.concatenate((directions_1, directions_2), axis = 1)

# Для разряжения массива np_array подготовим список индексов
range_1 = list(range(len(np_r_1_turn))) 
range_2 = list(range(len(np_r_1_turn), 2 * len(np_r_1_turn))) 
range_1 = get_rare(range_1)
range_2 = get_rare(range_2)
idx_range = range_1 + range_2

# Делаем массив точек разреженным
print('old len(np_r_2_turn)', len(np_r_2_turns))
np_r_2_turns = np_r_2_turns[idx_range]
np_fi_2_turns = np_fi_2_turns[idx_range]
np_z_2_turns = np_z_2_turns[idx_range]
directions = directions[idx_range]
print('new len(np_r_2_turn)', len(np_r_2_turns))

# Создаём массивы из num_turns витков
np_r = np.array([])
np_fi = np.array([])
np_z = np.array([])

# for i in range(int(num_turns / 2)):
for i in range(1):
    np_fi_2 = np_fi_2_turns if i == 0 else np_fi_2_turns + np_fi[-1] + shift_angle 
    np_r = np.concatenate((np_r, np_r_2_turns))
    np_fi = np.concatenate((np_fi, np_fi_2))
    np_z = np.concatenate((np_z, np_z_2_turns))

print('len(np_fi)', len(np_fi))
print('len(np_r)', len(np_r))
print('len(np_z)', len(np_z))

# -------------------------------------------------------------------------------------------
# ------- ПОСТРОЕНИЕ ТРАЕКТОРИИ ДВИЖЕНИЯ РАСКЛАДОЧНОЙ ГОЛОВЫ --------------------------------
# -------------------------------------------------------------------------------------------

# Траектрия движения раздаточной головы. Условно считаем, что идеальная траектория - эквидестанта поверхности.
# Эквидестанту получаем простым параллельным переносом
# По краям эквидестанты добавляем две точки для увеличения длины траектории
z_plane = -175 # высота плоскости, в которой движется кольцо раздаточной головы
# shift = 0 # сдвиг
# increasing_distance = 150 # Увеличиваем вылет траектории по краям
# extreme_y = 0
# first_point = [(-increasing_distance, extreme_y)] 
# curve_trajectory = [(x, y + shift) for x, y in zip(start_z_data, start_r_data)]
# end_point = [(start_z_data[-1] + increasing_distance, extreme_y)]
# trajectory_list = first_point + curve_trajectory + end_point

# Итоговая траектория близка к эквидестанте и подобрана вручную, опираясь на внутреннее чутьё
trajectory_list = [(-150, 0), (0, 70), (65, 90.0), (103, 95.0), (243, 118.0), (293, 118.0), (343, 118.0), (483, 95.0), (521, 90.0), (586, 70), (736, 0)]
print('trajectory_list', trajectory_list)

# -------------------------------------------------------------------------------------------
# ------- ПОСТРОЕНИЕ 3D МОДЕЛИ --------------------------------------------------------------
# -------------------------------------------------------------------------------------------

def rotation_body(radius_array, z_array, nt = 50):
    theta = np.linspace(0, 2*np.pi, nt)
    np_fi, r = np.meshgrid(theta, radius_array)
    _,  z = np.meshgrid(theta, z_array)
    x = r * np.cos(np_fi)
    y = r * np.sin(np_fi)
    z = z
    return x, y, z    


def cartesian_from_polar(r, np_fi, h):
    x = r * np.cos(np_fi)
    y = r * np.sin(np_fi)
    z = h
    return x, y, z


def inside(x, interval):
    if x >= interval[0] and x <= interval[1]:
        return True
    return False


def intersection_plane(point, ray, z_plain):
    t = (z_plain - point[2]) / ray[2]
    return point + t * ray


def get_distance_on_horizontal_plane(a, b, point):
    np_a = np.array(a)
    np_b = np.array(b)
    np_p = np.array(point) if type(point) != 'numpy.ndarray' else point
    np_p = np_p[:2] if len(np_p) > 2 else np_p
    length = np.linalg.norm
    if np.dot(np_b - np_a, np_p - np_a) > 0:
        if np.dot(np_a - np_b, np_p - np_b) > 0:
            dx, dy = np_a - np_b
            c = np_a[0] * np_b[1] - np_b[0] * np_a[1]
            return abs(dy * np_p[0] - dx * np_p[1] + c) / length(np_b - np_a)
        else:
            return length(np_p - np_b)
    else:
        return length(np_p - np_a)


def get_distance(segments, point):
    distances  = []
    for i in range(len(segments) - 1):
        line = (segments[i], segments[i + 1])
        distances.append(get_distance_on_horizontal_plane(*line, point))
    return min(distances)


# Точки для построения поверхности вращения
x1, y1, z1 = rotation_body(np_r_1_turn, np_z_1_turn)
surface = go.Surface(x=z1, y=y1, z=x1, colorscale='Viridis', opacity=0.5)

# Точки для построения геодезической линии
x2, y2, z2 = cartesian_from_polar(np_r, np_fi, np_z)
curve = go.Scatter3d(x=z2, y=y2, z=x2, mode='lines', line={'width': 10})

# Точки для построения траектории движения раскладочной головы
np_trajectory = np.array(trajectory_list)
z_column = np.array([z_plane] * len(np_trajectory))
np_trajectory = np.column_stack((np_trajectory, z_column)).T
x, y, z = np_trajectory[0], np_trajectory[1], np_trajectory[2]
trajectory = go.Scatter3d(x=x, y=y, z=z, mode='lines', line={'width': 10})

# Точки для построения границ трактории
x1 = np.min(np_trajectory[0])
x = np.array([x1, x1])
y = np.array([-100, 150])
z = np.array([z_plane, z_plane])
border_1 = go.Scatter3d(x=x, y=y, z=z, mode='lines', line={'width': 10}, marker={'color':'red'})
x2 = np.max(np_trajectory[0])
x = np.array([x2, x2])
border_2 = go.Scatter3d(x=x, y=y, z=z, mode='lines', line={'width': 10}, marker={'color':'red'})
borders = [border_1, border_2]
borders_interval = [x1.tolist(), x2.tolist()]


def get_intersect_point_exit_point_angle(start_polar, end_polar):
    dfi = (end_polar - start_polar)[1]
    start_beta = radians(45)
    finish_beta = radians(170)
    steps = 200
    dbeta = (finish_beta - start_beta) / steps;
    min_distance = None
    point_on_plain = None
    start_point = None
    angle = start_beta
    for step in range(steps):
        beta = start_beta + step * dbeta;
        z1, y1, x1 = cartesian_from_polar(start_polar[0], beta, start_polar[2])
        z2, y2, x2 = cartesian_from_polar(end_polar[0], beta + dfi, end_polar[2])
        np_start = np.array((x1, y1, z1))
        np_end = np.array((x2, y2, z2))
        np_ray = np_end - np_start
        intersection_point = intersection_plane(np_start, np_ray, z_plane)
        d = get_distance(trajectory_list, intersection_point)
        if step == 0:
            min_distance = d
            point_on_plain = intersection_point
            start_point = np_start
        else:
            if d < min_distance:
                min_distance = d
                point_on_plain = intersection_point
                start_point = np_start
                angle = beta
    return point_on_plain, start_point, angle, min_distance


# Точки для построения векторов схода нити
vectors = [] # Массив для построения линий схода нити в 3D
thread_points = [] # Массив некоторых координат линий схода нити для последующего получения координат управлющих органов станка
max_distance = 0; # Параметр исключительно для проверки результатов
for i, item in enumerate(directions):
    start_polar = item[:3]
    end_polar = item[3:]
    point_on_plain, start_point, angle, min_distance = get_intersect_point_exit_point_angle(start_polar, end_polar)
    x, y, fi = point_on_plain[0].tolist(), point_on_plain[1].tolist(), angle
    thread_points.append([x, y, fi])
    v = np.stack((start_point, point_on_plain)).T
    x, y, z = v[0], v[1], v[2]
    color = 'cyan' if inside(point_on_plain[0].tolist(), borders_interval) else 'coral'
    vectors.append(go.Scatter3d(x=x, y=y, z=z, mode='lines', line={'width': 10}, opacity=0.5, marker={'color':color}))
    max_distance = min_distance if min_distance > max_distance else max_distance
print('max_distance', max_distance) # параметр исключительно для проверки

# Создание html файла 3D модели
data = [surface, curve, trajectory, *borders, *vectors]
fig = go.Figure(data=data)
fig.update_layout(
    title={
        'text': f'Углы армирования:<br>на краю - {psi_0}<br>в центре - {round(degrees(psi_center))}',
        'y':0.9,
        'x':0.5,
        'xanchor': 'center',
        'yanchor': 'top'})
fig.write_html('tmp.html', auto_open=True) # Для получения графиков убираем # с этой строки


# -------------------------------------------------------------------------------------------
# ------- ПОЛУЧАЕМ МАССИВ КООРДИНАТ УПРАВЛЯЮЩИХ ОРГАНОВ СТАНКА ------------------------------
# -------------------------------------------------------------------------------------------

